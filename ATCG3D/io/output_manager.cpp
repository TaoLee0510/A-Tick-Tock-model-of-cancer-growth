#include "io/output_manager.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string_view>

#include "io/preview_sampler.hpp"
#include "io/snapshot.hpp"
#include "io/vtkhdf_writer.hpp"
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
#include "io/checkpoint_hdf5.hpp"
#endif

namespace atcg3d {
namespace {

bool same_time(double lhs, double rhs) {
    return std::abs(lhs - rhs) <= 1e-10 * std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

std::string frame_name(std::size_t index) {
    std::ostringstream value;
    value << "frame_" << std::setw(8) << std::setfill('0') << index << ".vtkhdf";
    return value.str();
}

#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
std::string checkpoint_name(std::uint64_t completed_events) {
    std::ostringstream value;
    value << "checkpoint_" << std::setw(16) << std::setfill('0') << completed_events << ".h5";
    return value.str();
}
#endif

void require_unused_path(const std::filesystem::path& path) {
    if (std::filesystem::exists(path) ||
        std::filesystem::exists(std::filesystem::path(path.string() + ".tmp"))) {
        throw std::runtime_error("refusing to overwrite existing output: " + path.string());
    }
}

template <class Integer>
Integer parse_integer(std::string_view value,
                      const std::filesystem::path& path,
                      std::string_view field) {
    Integer result{};
    const char* begin = value.data();
    const char* end = begin + value.size();
    const auto parsed = std::from_chars(begin, end, result);
    if (value.empty() || parsed.ec != std::errc{} || parsed.ptr != end) {
        throw std::runtime_error("invalid lineage " + std::string(field) + " in " +
                                 path.string());
    }
    return result;
}

double parse_finite_double(std::string_view value,
                           const std::filesystem::path& path,
                           std::string_view field) {
    std::size_t consumed = 0;
    double result{};
    try {
        result = std::stod(std::string(value), &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid lineage " + std::string(field) + " in " +
                                 path.string());
    }
    if (consumed != value.size() || !std::isfinite(result)) {
        throw std::runtime_error("invalid lineage " + std::string(field) + " in " +
                                 path.string());
    }
    return result;
}

std::array<std::string_view, 5> split_lineage_row(const std::string& line,
                                                  const std::filesystem::path& path) {
    std::array<std::string_view, 5> fields;
    const std::string_view row(line);
    std::size_t begin = 0;
    for (std::size_t index = 0; index < fields.size(); ++index) {
        const std::size_t comma = row.find(',', begin);
        if ((index + 1 == fields.size()) != (comma == std::string_view::npos)) {
            throw std::runtime_error("invalid lineage row in " + path.string());
        }
        const std::size_t end = comma == std::string_view::npos ? row.size() : comma;
        fields[index] = row.substr(begin, end - begin);
        begin = end + 1;
    }
    return fields;
}

std::vector<LineageEdge> read_existing_lineage(const std::filesystem::path& path) {
    if (!std::filesystem::exists(path)) return {};
    if (!std::filesystem::is_regular_file(path)) {
        throw std::runtime_error("lineage path is not a regular file: " + path.string());
    }
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("unable to read lineage: " + path.string());
    std::string line;
    if (!std::getline(stream, line) ||
        line != "birth_time,child_uid,parent_uid,clone_id,type") {
        throw std::runtime_error("invalid lineage header: " + path.string());
    }
    std::vector<LineageEdge> edges;
    while (std::getline(stream, line)) {
        if (line.empty()) throw std::runtime_error("empty lineage row: " + path.string());
        const auto fields = split_lineage_row(line, path);
        LineageEdge edge;
        edge.birth_time = parse_finite_double(fields[0], path, "birth_time");
        if (edge.birth_time < 0.0) {
            throw std::runtime_error("negative lineage birth_time in " + path.string());
        }
        edge.child_uid = parse_integer<CellUid>(fields[1], path, "child_uid");
        edge.parent_uid = parse_integer<CellUid>(fields[2], path, "parent_uid");
        edge.clone_id = parse_integer<std::uint32_t>(fields[3], path, "clone_id");
        const unsigned int type = parse_integer<unsigned int>(fields[4], path, "type");
        if (type != static_cast<unsigned int>(CellType::r) &&
            type != static_cast<unsigned int>(CellType::K)) {
            throw std::runtime_error("invalid lineage cell type in " + path.string());
        }
        edge.type = static_cast<CellType>(type);
        edges.push_back(edge);
    }
    if (!stream.eof()) throw std::runtime_error("unable to finish lineage read: " + path.string());
    return edges;
}

template <class Writer>
void write_text_atomic(const std::filesystem::path& path, Writer writer) {
    std::filesystem::create_directories(path.parent_path());
    const std::filesystem::path temporary = path.string() + ".tmp";
    std::filesystem::remove(temporary);
    {
        std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
        if (!stream) throw std::runtime_error("unable to create " + temporary.string());
        writer(stream);
        stream.flush();
        if (!stream) throw std::runtime_error("unable to finish " + temporary.string());
    }
    std::filesystem::rename(temporary, path);
}

bool same_lineage_edge(const LineageEdge& lhs, const LineageEdge& rhs) noexcept {
    return lhs.birth_time == rhs.birth_time && lhs.child_uid == rhs.child_uid &&
           lhs.parent_uid == rhs.parent_uid && lhs.clone_id == rhs.clone_id &&
           lhs.type == rhs.type;
}

void write_lineage_atomic(const std::filesystem::path& path,
                          const std::vector<LineageEdge>& edges,
                          std::size_t begin,
                          std::size_t end) {
    if (begin > end || end > edges.size()) {
        throw std::logic_error("invalid lineage recovery range");
    }
    write_text_atomic(path, [&](std::ostream& out) {
        out << "birth_time,child_uid,parent_uid,clone_id,type\n" << std::setprecision(17);
        for (std::size_t index = begin; index < end; ++index) {
            const LineageEdge& edge = edges[index];
            out << edge.birth_time << ',' << edge.child_uid << ',' << edge.parent_uid << ','
                << edge.clone_id << ',' << static_cast<unsigned int>(edge.type) << '\n';
        }
    });
}

struct ParsedArtifactName {
    std::uint64_t index{};
    bool temporary{};
};

std::optional<ParsedArtifactName> parse_indexed_artifact(
    std::string_view name,
    std::string_view prefix,
    std::size_t digits,
    std::string_view extension) {
    bool temporary = false;
    constexpr std::string_view temporary_suffix = ".tmp";
    if (name.size() >= temporary_suffix.size() &&
        name.substr(name.size() - temporary_suffix.size()) == temporary_suffix) {
        temporary = true;
        name.remove_suffix(temporary_suffix.size());
    }
    if (name.size() != prefix.size() + digits + extension.size() ||
        name.substr(0, prefix.size()) != prefix ||
        name.substr(prefix.size() + digits) != extension) {
        return std::nullopt;
    }
    const std::string_view number = name.substr(prefix.size(), digits);
    std::uint64_t index{};
    const auto parsed = std::from_chars(number.data(), number.data() + number.size(), index);
    if (parsed.ec != std::errc{} || parsed.ptr != number.data() + number.size()) {
        return std::nullopt;
    }
    return ParsedArtifactName{index, temporary};
}

struct RecoveryArtifact {
    std::filesystem::path source;
    std::filesystem::path relative;
    bool safe_before_catalog_update{};
};

void collect_if_present(std::vector<RecoveryArtifact>& artifacts,
                        const std::filesystem::path& run_directory,
                        const std::filesystem::path& relative,
                        bool safe_before_catalog_update) {
    const std::filesystem::path source = run_directory / relative;
    if (std::filesystem::exists(source)) {
        artifacts.push_back({source, relative, safe_before_catalog_update});
    }
}

void collect_frame_artifacts(std::vector<RecoveryArtifact>& artifacts,
                             const std::filesystem::path& run_directory,
                             std::string_view relative_directory,
                             std::size_t kept_count) {
    const std::filesystem::path directory = run_directory / relative_directory;
    if (!std::filesystem::is_directory(directory)) return;
    for (const std::filesystem::directory_entry& entry :
         std::filesystem::directory_iterator(directory)) {
        const auto parsed = parse_indexed_artifact(
            entry.path().filename().string(), "frame_", 8, ".vtkhdf");
        if (!parsed || (!parsed->temporary && parsed->index < kept_count)) continue;
        artifacts.push_back({entry.path(), entry.path().lexically_relative(run_directory),
                             parsed->temporary});
    }
}

void collect_checkpoint_artifacts(std::vector<RecoveryArtifact>& artifacts,
                                  const std::filesystem::path& run_directory,
                                  std::uint64_t completed_events) {
    const std::filesystem::path directory = run_directory / "checkpoints";
    if (!std::filesystem::is_directory(directory)) return;
    for (const std::filesystem::directory_entry& entry :
         std::filesystem::directory_iterator(directory)) {
        const auto parsed = parse_indexed_artifact(
            entry.path().filename().string(), "checkpoint_", 16, ".h5");
        if (!parsed || (!parsed->temporary && parsed->index <= completed_events)) continue;
        artifacts.push_back({entry.path(), entry.path().lexically_relative(run_directory),
                             parsed->temporary});
    }
}

std::size_t checkpoint_prefix_size(const std::vector<SeriesEntry3D>& entries,
                                   double checkpoint_time) {
    std::size_t size = 0;
    while (size < entries.size() &&
           (entries[size].time_hours < checkpoint_time ||
            same_time(entries[size].time_hours, checkpoint_time))) {
        ++size;
    }
    return size;
}

std::filesystem::path create_recovery_directory(
    const std::filesystem::path& run_directory,
    std::uint64_t completed_events) {
    const std::filesystem::path parent = run_directory / "recovery";
    std::filesystem::create_directories(parent);
    std::ostringstream stem;
    stem << "checkpoint_" << std::setw(16) << std::setfill('0') << completed_events;
    for (std::size_t attempt = 0; attempt < 10000; ++attempt) {
        std::ostringstream name;
        name << stem.str();
        if (attempt != 0) name << '_' << std::setw(4) << std::setfill('0') << attempt;
        const std::filesystem::path candidate = parent / name.str();
        if (std::filesystem::create_directory(candidate)) return candidate;
    }
    throw std::runtime_error("unable to allocate a unique output recovery directory");
}

void quarantine_artifact(const RecoveryArtifact& artifact,
                         const std::filesystem::path& recovery_directory,
                         std::vector<std::string>& quarantined) {
    const std::filesystem::path destination = recovery_directory / artifact.relative;
    std::filesystem::create_directories(destination.parent_path());
    if (std::filesystem::exists(destination)) {
        throw std::runtime_error("refusing to overwrite recovery artifact: " +
                                 destination.string());
    }
    std::filesystem::rename(artifact.source, destination);
    quarantined.push_back(artifact.relative.generic_string());
}

}  // namespace

OutputManager3D::OutputManager3D(const Model3DConfig& config)
    : config_(config), run_directory_(config.output_directory) {
    if (!config_.output_enabled) return;
    if ((config_.preview_every_hours > 0.0 || config_.full_every_hours > 0.0) &&
        !vtkhdf_output_available()) {
        throw std::runtime_error(
            "3D visualization output was requested but VTK-HDF support is not built");
    }
#ifndef ATCG3D_HAS_HDF5_CHECKPOINT
    if (config_.checkpoint_every_hours > 0.0) {
        throw std::runtime_error(
            "checkpoint output was requested but HDF5 checkpoint support is not built");
    }
#endif
    resume_mode_ = config_.run_mode == "resume";
    if (resume_mode_) {
        ExistingRunOutput3D existing = load_existing_run_output(run_directory_, config_);
        preview_ = std::move(existing.preview);
        full_ = std::move(existing.full);
        vessels_ = std::move(existing.vessels);
        if (!preview_.empty()) last_preview_time_ = preview_.back().time_hours;
        if (!full_.empty()) last_full_time_ = full_.back().time_hours;
        existing_lineage_ =
            read_existing_lineage(run_directory_ / "lineage" / "edges.csv");
        lineage_written_ = existing_lineage_.size();
        return;
    }
    if (std::filesystem::exists(run_directory_) &&
        (!std::filesystem::is_directory(run_directory_) ||
         std::filesystem::directory_iterator(run_directory_) !=
             std::filesystem::directory_iterator{})) {
        throw std::runtime_error("refusing to overwrite existing run directory: " +
                                 run_directory_.string());
    }
    std::filesystem::create_directories(run_directory_ / "metrics");
    std::filesystem::create_directories(run_directory_ / "checkpoints");
    std::filesystem::create_directories(run_directory_ / "lineage");
    std::filesystem::create_directories(run_directory_ / "viz" / "preview");
    std::filesystem::create_directories(run_directory_ / "viz" / "full");
    std::filesystem::create_directories(run_directory_ / "viz" / "vessels");
    resume_validated_ = true;
    update_metadata();
}

bool OutputManager3D::due(double now, double next, double interval) const noexcept {
    return interval > 0.0 && (now > next || same_time(now, next));
}

double OutputManager3D::advance(double next, double now, double interval) noexcept {
    if (interval <= 0.0) return next;
    do {
        next += interval;
    } while (now > next || same_time(now, next));
    return next;
}

void OutputManager3D::observe(const Simulation3D& simulation) {
    if (!config_.output_enabled) return;
    if (!resume_validated_) validate_resume_state(simulation);
    const double now = simulation.clock().time_hours;
    const bool full_due = due(now, next_full_, config_.full_every_hours);
    const bool preview_due = due(now, next_preview_, config_.preview_every_hours);
    if ((preview_due || full_due) && !same_time(last_preview_time_, now)) {
        write_preview(simulation);
    }
    if (full_due && !same_time(last_full_time_, now)) {
        write_full(simulation);
    }
    if (due(now, next_checkpoint_, config_.checkpoint_every_hours)) {
        write_checkpoint(simulation);
    }
    if (preview_due) next_preview_ = advance(next_preview_, now, config_.preview_every_hours);
    if (full_due) next_full_ = advance(next_full_, now, config_.full_every_hours);
    if (due(now, next_checkpoint_, config_.checkpoint_every_hours)) {
        next_checkpoint_ = advance(next_checkpoint_, now, config_.checkpoint_every_hours);
    }
    append_lineage(simulation);
}

void OutputManager3D::write_preview(const Simulation3D& simulation) {
    const std::vector<Slot> slots = stable_preview_sample(
        simulation.cells(), static_cast<std::size_t>(config_.preview_max_cells),
        config_.preview_seed);
    const SimulationSnapshotView3D snapshot{
        simulation.cells(), slots, simulation.clock(), config_.display_radius};
    const std::string name = frame_name(preview_.size());
    const std::filesystem::path path = run_directory_ / "viz" / "preview" / name;
    require_unused_path(path);
    write_vtkhdf_points_atomic(path, snapshot);
    write_vessels(simulation, name);
    preview_.push_back({"viz/preview/" + name, simulation.clock().time_hours});
    last_preview_time_ = simulation.clock().time_hours;
    update_metadata();
}

void OutputManager3D::write_vessels(const Simulation3D& simulation,
                                    const std::string& frame) {
    const std::filesystem::path path = run_directory_ / "viz" / "vessels" / frame;
    require_unused_path(path);
    write_vtkhdf_vessels_atomic(path, simulation.vessel_nodes());
    vessels_.push_back({"viz/vessels/" + frame, simulation.clock().time_hours});
}

void OutputManager3D::write_full(const Simulation3D& simulation) {
    const std::vector<Slot> slots = simulation.cells().alive_slots();
    const SimulationSnapshotView3D snapshot{
        simulation.cells(), slots, simulation.clock(), config_.display_radius};
    const std::string name = frame_name(full_.size());
    const std::filesystem::path path = run_directory_ / "viz" / "full" / name;
    require_unused_path(path);
    write_vtkhdf_points_atomic(path, snapshot);
    full_.push_back({"viz/full/" + name, simulation.clock().time_hours});
    last_full_time_ = simulation.clock().time_hours;
    write_series_atomic(run_directory_ / "full.vtkhdf.series", full_);
    update_metadata();
}

void OutputManager3D::write_checkpoint(const Simulation3D& simulation) {
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
    const std::filesystem::path path = run_directory_ / "checkpoints" /
                                       checkpoint_name(simulation.clock().completed_events);
    require_unused_path(path);
    write_hdf5_checkpoint(path, simulation);
#else
    (void)simulation;
    throw std::runtime_error("HDF5 checkpoint support is not built");
#endif
}

void OutputManager3D::validate_resume_state(const Simulation3D& simulation) {
    const double now = simulation.clock().time_hours;
    if (!std::isfinite(now) || now < 0.0) {
        throw std::runtime_error("resume simulation has an invalid clock");
    }
    const auto& lineage = simulation.lineage();
    const std::size_t common_lineage = std::min(lineage.size(), existing_lineage_.size());
    for (std::size_t index = 0; index < common_lineage; ++index) {
        if (!same_lineage_edge(lineage[index], existing_lineage_[index])) {
            throw std::runtime_error(
                "existing lineage does not match the resume checkpoint prefix");
        }
    }

    const std::size_t preview_prefix = checkpoint_prefix_size(preview_, now);
    const std::size_t full_prefix = checkpoint_prefix_size(full_, now);
    const std::size_t vessel_prefix = checkpoint_prefix_size(vessels_, now);
    if (preview_prefix != vessel_prefix) {
        throw std::runtime_error(
            "preview and vessel output have different checkpoint prefixes");
    }

    const bool catalogs_trimmed = preview_prefix != preview_.size() ||
                                  full_prefix != full_.size() ||
                                  vessel_prefix != vessels_.size();
    const bool lineage_trimmed = existing_lineage_.size() > lineage.size();

    std::vector<RecoveryArtifact> artifacts;
    collect_if_present(artifacts, run_directory_, "preview.vtkhdf.series.tmp", true);
    collect_if_present(artifacts, run_directory_, "full.vtkhdf.series.tmp", true);
    collect_if_present(artifacts, run_directory_, "vessels.vtkhdf.series.tmp", true);
    collect_if_present(artifacts, run_directory_, "run.json.tmp", true);
    collect_if_present(artifacts, run_directory_, "lineage/edges.csv.tmp", true);
    collect_if_present(artifacts, run_directory_, "metrics/final.json.tmp", true);
    collect_frame_artifacts(artifacts, run_directory_, "viz/preview", preview_prefix);
    collect_frame_artifacts(artifacts, run_directory_, "viz/full", full_prefix);
    collect_frame_artifacts(artifacts, run_directory_, "viz/vessels", vessel_prefix);
    collect_checkpoint_artifacts(artifacts, run_directory_,
                                 simulation.clock().completed_events);

    const bool has_non_temporary_artifact = std::any_of(
        artifacts.begin(), artifacts.end(), [](const RecoveryArtifact& artifact) {
            return !artifact.safe_before_catalog_update;
        });
    const bool state_rollback = catalogs_trimmed || lineage_trimmed ||
                                has_non_temporary_artifact;
    if (state_rollback) {
        collect_if_present(artifacts, run_directory_, "metrics/final.json", true);
    }

    if (catalogs_trimmed || lineage_trimmed || !artifacts.empty()) {
        const std::filesystem::path recovery_directory = create_recovery_directory(
            run_directory_, simulation.clock().completed_events);
        std::vector<std::string> quarantined;

        // Temporary and final-metrics files are never referenced by a live
        // catalog, so preserve them before any atomic rewrite would reuse the
        // same .tmp name.
        for (const RecoveryArtifact& artifact : artifacts) {
            if (artifact.safe_before_catalog_update) {
                quarantine_artifact(artifact, recovery_directory, quarantined);
            }
        }

        if (lineage_trimmed) {
            write_lineage_atomic(recovery_directory / "lineage" / "edges_tail.csv",
                                 existing_lineage_, lineage.size(),
                                 existing_lineage_.size());
            quarantined.push_back("lineage/edges_tail.csv");
            write_lineage_atomic(run_directory_ / "lineage" / "edges.csv",
                                 existing_lineage_, 0, lineage.size());
            existing_lineage_.resize(lineage.size());
        }
        lineage_written_ = existing_lineage_.size();

        if (catalogs_trimmed) {
            preview_.resize(preview_prefix);
            full_.resize(full_prefix);
            vessels_.resize(vessel_prefix);
            update_metadata();
        }

        // Catalogs now describe only the verified checkpoint prefix. Moving
        // superseded final frames cannot leave a published catalog dangling.
        for (const RecoveryArtifact& artifact : artifacts) {
            if (!artifact.safe_before_catalog_update) {
                quarantine_artifact(artifact, recovery_directory, quarantined);
            }
        }

        write_text_atomic(recovery_directory / "recovery.json", [&](std::ostream& out) {
            out << std::setprecision(17)
                << "{\n"
                << "  \"schema\": \"atcg.output-recovery\",\n"
                << "  \"schema_version\": 1,\n"
                << "  \"checkpoint_time_hours\": " << now << ",\n"
                << "  \"checkpoint_completed_events\": "
                << simulation.clock().completed_events << ",\n"
                << "  \"kept_preview_frames\": " << preview_prefix << ",\n"
                << "  \"kept_full_frames\": " << full_prefix << ",\n"
                << "  \"kept_vessel_frames\": " << vessel_prefix << ",\n"
                << "  \"kept_lineage_edges\": " << lineage_written_ << ",\n"
                << "  \"quarantined\": [\n";
            for (std::size_t index = 0; index < quarantined.size(); ++index) {
                out << "    \"" << quarantined[index] << '"';
                if (index + 1 != quarantined.size()) out << ',';
                out << '\n';
            }
            out << "  ]\n}\n";
        });
    } else {
        lineage_written_ = existing_lineage_.size();
    }

    last_preview_time_ = preview_.empty() ? -1.0 : preview_.back().time_hours;
    last_full_time_ = full_.empty() ? -1.0 : full_.back().time_hours;

    // A resume never re-emits the restored instant. Continue each periodic
    // schedule at the first regular boundary strictly after the restored time.
    next_preview_ = advance(0.0, now, config_.preview_every_hours);
    next_full_ = advance(0.0, now, config_.full_every_hours);
    next_checkpoint_ = advance(0.0, now, config_.checkpoint_every_hours);
    resume_time_ = now;
    resume_validated_ = true;
}

void OutputManager3D::append_lineage(const Simulation3D& simulation) {
    const auto& lineage = simulation.lineage();
    if (lineage_written_ >= lineage.size()) return;
    const std::filesystem::path path = run_directory_ / "lineage" / "edges.csv";
    const bool write_header = !std::filesystem::exists(path);
    std::ofstream stream(path, std::ios::app);
    if (!stream) throw std::runtime_error("unable to append lineage: " + path.string());
    if (write_header) stream << "birth_time,child_uid,parent_uid,clone_id,type\n";
    stream << std::setprecision(17);
    for (; lineage_written_ < lineage.size(); ++lineage_written_) {
        const LineageEdge& edge = lineage[lineage_written_];
        stream << edge.birth_time << ',' << edge.child_uid << ',' << edge.parent_uid << ','
               << edge.clone_id << ',' << static_cast<unsigned int>(edge.type) << '\n';
    }
    stream.flush();
    if (!stream) throw std::runtime_error("unable to finish lineage: " + path.string());
}

void OutputManager3D::update_metadata() {
    // Publish dependent catalogs before preview, which is the viewer's
    // timeline. A live reader that discovers a new preview time can therefore
    // already resolve its exact vessel frame (and any exact full frame).
    write_series_atomic(run_directory_ / "vessels.vtkhdf.series", vessels_);
    write_series_atomic(run_directory_ / "full.vtkhdf.series", full_);
    write_series_atomic(run_directory_ / "preview.vtkhdf.series", preview_);
    write_run_manifest_atomic(run_directory_, config_, preview_, full_, vessels_,
                              config_.checkpoint_every_hours > 0.0);
}

void OutputManager3D::write_metrics(const Simulation3D& simulation) {
    write_text_atomic(run_directory_ / "metrics" / "final.json", [&](std::ostream& out) {
        out << std::setprecision(17)
            << "{\n"
            << "  \"time_hours\": " << simulation.clock().time_hours << ",\n"
            << "  \"completed_events\": " << simulation.clock().completed_events << ",\n"
            << "  \"alive_cells\": " << simulation.cells().alive_count() << ",\n"
            << "  \"vessel_nodes\": " << simulation.vessel_nodes().alive_count() << ",\n"
            << "  \"cell_store_bytes\": " << simulation.cells().allocated_bytes() << ",\n"
            << "  \"grid_bytes\": " << simulation.grid().allocated_bytes() << ",\n"
            << "  \"state_checksum\": " << simulation.state_checksum() << "\n"
            << "}\n";
    });
}

void OutputManager3D::finalize(const Simulation3D& simulation) {
    if (!config_.output_enabled || finalized_) return;
    observe(simulation);
    const bool advanced_since_resume =
        !resume_mode_ || (simulation.clock().time_hours > resume_time_ &&
                          !same_time(simulation.clock().time_hours, resume_time_));
    if (!same_time(last_preview_time_, simulation.clock().time_hours) &&
        advanced_since_resume &&
        (config_.preview_every_hours > 0.0 || config_.full_every_hours > 0.0)) {
        write_preview(simulation);
    }
    append_lineage(simulation);
    write_metrics(simulation);
    update_metadata();
    finalized_ = true;
}

}  // namespace atcg3d
