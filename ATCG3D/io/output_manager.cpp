#include "io/output_manager.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <charconv>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <optional>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string_view>

#include "io/preview_sampler.hpp"
#include "io/snapshot.hpp"
#include "io/vtkhdf_writer.hpp"
#include "engine/parallelism.hpp"
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

std::string checkpoint_name(const SimulationClock3D& clock) {
    std::ostringstream value;
    value << "checkpoint_" << std::setw(16) << std::setfill('0')
          << clock.completed_events << "_time_" << std::hex << std::setw(16)
          << std::setfill('0') << std::bit_cast<std::uint64_t>(clock.time_hours)
          << ".h5";
    return value.str();
}

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

struct ParsedCheckpointArtifactName {
    std::uint64_t completed_events{};
    std::optional<double> time_hours;
    bool temporary{};
};

std::optional<ParsedCheckpointArtifactName> parse_checkpoint_artifact(
    std::string_view name) {
    bool temporary = false;
    constexpr std::string_view temporary_suffix = ".tmp";
    if (name.size() >= temporary_suffix.size() &&
        name.substr(name.size() - temporary_suffix.size()) == temporary_suffix) {
        temporary = true;
        name.remove_suffix(temporary_suffix.size());
    }

    constexpr std::string_view prefix = "checkpoint_";
    constexpr std::string_view separator = "_time_";
    constexpr std::string_view extension = ".h5";
    constexpr std::size_t event_digits = 16;
    constexpr std::size_t time_hex_digits = 16;
    const auto parse_events = [&](std::string_view number)
        -> std::optional<std::uint64_t> {
        std::uint64_t value{};
        const auto parsed = std::from_chars(
            number.data(), number.data() + number.size(), value);
        if (parsed.ec != std::errc{} ||
            parsed.ptr != number.data() + number.size()) return std::nullopt;
        return value;
    };

    const std::size_t legacy_size = prefix.size() + event_digits + extension.size();
    if (name.size() == legacy_size && name.substr(0, prefix.size()) == prefix &&
        name.substr(prefix.size() + event_digits) == extension) {
        const auto events = parse_events(name.substr(prefix.size(), event_digits));
        if (!events) return std::nullopt;
        return ParsedCheckpointArtifactName{*events, std::nullopt, temporary};
    }

    const std::size_t timed_size = prefix.size() + event_digits + separator.size() +
                                   time_hex_digits + extension.size();
    if (name.size() != timed_size || name.substr(0, prefix.size()) != prefix ||
        name.substr(prefix.size() + event_digits, separator.size()) != separator ||
        name.substr(name.size() - extension.size()) != extension) {
        return std::nullopt;
    }
    const auto events = parse_events(name.substr(prefix.size(), event_digits));
    if (!events) return std::nullopt;
    const std::size_t time_begin = prefix.size() + event_digits + separator.size();
    const std::string_view time_hex = name.substr(time_begin, time_hex_digits);
    std::uint64_t time_bits{};
    const auto parsed_time = std::from_chars(
        time_hex.data(), time_hex.data() + time_hex.size(), time_bits, 16);
    if (parsed_time.ec != std::errc{} ||
        parsed_time.ptr != time_hex.data() + time_hex.size()) return std::nullopt;
    const double time_hours = std::bit_cast<double>(time_bits);
    if (!std::isfinite(time_hours) || time_hours < 0.0) return std::nullopt;
    return ParsedCheckpointArtifactName{*events, time_hours, temporary};
}

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
                                  std::uint64_t completed_events,
                                  double checkpoint_time) {
    const std::filesystem::path directory = run_directory / "checkpoints";
    if (!std::filesystem::is_directory(directory)) return;
    for (const std::filesystem::directory_entry& entry :
         std::filesystem::directory_iterator(directory)) {
        const auto parsed = parse_checkpoint_artifact(
            entry.path().filename().string());
        if (!parsed) continue;
        const bool within_checkpoint = parsed->time_hours
            ? (*parsed->time_hours < checkpoint_time ||
               same_time(*parsed->time_hours, checkpoint_time))
            : parsed->completed_events <= completed_events;
        if (!parsed->temporary && within_checkpoint) continue;
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
        last_scheduled_preview_time_ = last_preview_time_;
        last_scheduled_full_time_ = last_full_time_;
        start_worker();
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
    start_worker();
}

OutputManager3D::~OutputManager3D() noexcept {
    if (!worker_.joinable()) return;
    {
        std::lock_guard lock(queue_mutex_);
        worker_stopping_ = true;
    }
    queue_ready_.notify_all();
    worker_.join();
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

OutputManager3D::FrozenSimulationState3D OutputManager3D::freeze(
    const Simulation3D& simulation) const {
    FrozenSimulationState3D state;
    state.cell_slots = simulation.snapshot_cell_slots();
    state.cells.resize(state.cell_slots.size());
    state.lesion_ids.resize(state.cell_slots.size());
    deterministic_parallel_for(
        state.cell_slots.size(), std::max(1, config_.threads),
        [&](std::size_t index) {
            const Slot slot = state.cell_slots[index];
            state.cells[index] = simulation.cells().snapshot(slot);
            state.lesion_ids[index] = simulation.lesion_index()
                .lesion_for_anchor(state.cells[index].anchor)
                .value_or(kNoLesionId);
        });
    state.free_slots = simulation.cells().free_slots();
    state.next_uid = simulation.next_uid();
    state.clock = simulation.clock();
    state.stats = simulation.stats();
    state.lineage = simulation.lineage();
    state.vasculature = simulation.snapshot_vasculature();
    state.cell_slot_count = simulation.cells().slot_count();
    state.state_checksum = simulation.state_checksum();
    return state;
}

OutputManager3D::FrozenPreviewState3D OutputManager3D::freeze_preview(
    const Simulation3D& simulation) const {
    FrozenPreviewState3D state;
    state.cell_slots = stable_preview_sample(
        simulation.cells(),
        static_cast<std::size_t>(config_.preview_max_cells),
        config_.preview_seed);
    state.cells.resize(state.cell_slots.size());
    state.lesion_ids.resize(state.cell_slots.size());
    deterministic_parallel_for(
        state.cell_slots.size(), std::max(1, config_.threads),
        [&](std::size_t index) {
            const Slot slot = state.cell_slots[index];
            state.cells[index] = simulation.cells().snapshot(slot);
            state.lesion_ids[index] = simulation.lesion_index()
                .lesion_for_anchor(state.cells[index].anchor)
                .value_or(kNoLesionId);
        });
    state.clock = simulation.clock();
    state.stats = simulation.stats();
    state.vasculature = simulation.snapshot_vasculature();
    state.cell_slot_count = simulation.cells().slot_count();
    return state;
}

bool OutputManager3D::checkpoint_base_due(double now) const noexcept {
    return config_.storage_mode == "self_contained_v1" ||
           checkpoint_parent_path_.empty() ||
           checkpoint_delta_chain_ >=
               config_.checkpoint_max_delta_chain ||
           checkpoint_base_time_ < 0.0 ||
           now - checkpoint_base_time_ + 1e-10 >=
               config_.checkpoint_base_every_hours;
}

OutputManager3D::FrozenCheckpointJournal3D
OutputManager3D::freeze_checkpoint_journal(
    const Simulation3D& simulation,
    const std::filesystem::path& path) {
    if (checkpoint_parent_path_.empty() ||
        checkpoint_parent_time_ < 0.0) {
        throw std::logic_error(
            "checkpoint journal has no scheduled parent");
    }
    const auto& lineage = simulation.lineage();
    if (checkpoint_lineage_scheduled_ > lineage.size()) {
        throw std::logic_error(
            "checkpoint lineage prefix exceeds current lineage");
    }
    FrozenCheckpointJournal3D state;
    state.cells = simulation.cells().take_checkpoint_journal();
    state.next_uid = simulation.next_uid();
    state.clock = simulation.clock();
    state.stats = simulation.stats();
    state.lineage_prefix_count = checkpoint_lineage_scheduled_;
    state.lineage_tail.assign(
        lineage.begin() +
            static_cast<std::ptrdiff_t>(checkpoint_lineage_scheduled_),
        lineage.end());
    state.vasculature = simulation.snapshot_vasculature();
    state.state_checksum = simulation.state_checksum();
    state.path = path;
    state.parent_path = checkpoint_parent_path_;
    state.parent_state_checksum =
        checkpoint_parent_state_checksum_;
    state.parent_time_hours = checkpoint_parent_time_;
    state.chain_length = checkpoint_delta_chain_ + 1;
    return state;
}

void OutputManager3D::start_worker() {
    if (!config_.output_enabled || !config_.output_async_enabled ||
        worker_.joinable()) return;
    worker_ = std::thread([this] { worker_loop(); });
}

void OutputManager3D::throw_worker_error() {
    if (!worker_failed_.load(std::memory_order_acquire)) return;
    std::exception_ptr error;
    {
        std::lock_guard lock(queue_mutex_);
        error = worker_error_;
    }
    if (error) std::rethrow_exception(error);
}

std::uint64_t OutputManager3D::estimate_job_bytes(
    const OutputJob3D& job) noexcept {
    std::uint64_t bytes = sizeof(OutputJob3D);
    if (job.preview_state) {
        bytes += job.preview_state->cells.capacity() * sizeof(CellInit);
        bytes += job.preview_state->cell_slots.capacity() * sizeof(Slot);
        bytes +=
            job.preview_state->lesion_ids.capacity() * sizeof(LesionId);
    }
    if (job.full_state) {
        bytes += job.full_state->cells.capacity() * sizeof(CellInit);
        bytes += job.full_state->cell_slots.capacity() * sizeof(Slot);
        bytes += job.full_state->free_slots.capacity() * sizeof(Slot);
        bytes += job.full_state->lesion_ids.capacity() * sizeof(LesionId);
        bytes +=
            job.full_state->lineage.capacity() * sizeof(LineageEdge);
    }
    if (job.checkpoint_journal) {
        bytes += job.checkpoint_journal->cells.mutations.capacity() *
                 sizeof(CheckpointCellMutation3D);
        bytes +=
            job.checkpoint_journal->cells.free_list_mutations.capacity() *
            sizeof(FreeListMutation3D);
        bytes +=
            job.checkpoint_journal->lineage_tail.capacity() *
            sizeof(LineageEdge);
    }
    return bytes;
}

void OutputManager3D::enqueue(OutputJob3D job) {
    const std::uint64_t job_bytes = estimate_job_bytes(job);
    std::unique_lock lock(queue_mutex_);
    if (job.live_preview &&
        config_.preview_overflow_policy == "coalesce_latest") {
        const auto replaceable = std::find_if(
            jobs_.rbegin(), jobs_.rend(), [](const OutputJob3D& queued) {
                return queued.live_preview && !queued.preview &&
                       !queued.full && !queued.checkpoint_base &&
                       !queued.checkpoint_journal;
            });
        if (replaceable != jobs_.rend()) {
            const std::uint64_t old_bytes =
                estimate_job_bytes(*replaceable);
            *replaceable = std::move(job);
            queued_snapshot_bytes_ =
                old_bytes > queued_snapshot_bytes_
                ? job_bytes
                : queued_snapshot_bytes_ - old_bytes + job_bytes;
            lock.unlock();
            queue_ready_.notify_one();
            return;
        }
    }
    queue_space_.wait(lock, [this, job_bytes] {
        const bool within_depth =
            jobs_.size() < config_.output_async_queue_depth;
        const bool within_bytes =
            queued_snapshot_bytes_ <=
            config_.output_async_max_pending_bytes -
                std::min(job_bytes,
                         config_.output_async_max_pending_bytes);
        // A single large full snapshot is allowed through an otherwise-empty
        // queue; the byte cap prevents several such snapshots accumulating.
        return (within_depth && (within_bytes || jobs_.empty())) ||
               worker_error_ || worker_stopping_;
    });
    if (worker_error_) {
        const std::exception_ptr error = worker_error_;
        lock.unlock();
        std::rethrow_exception(error);
    }
    if (worker_stopping_) {
        throw std::runtime_error("asynchronous output worker is stopping");
    }
    jobs_.push_back(std::move(job));
    queued_snapshot_bytes_ += job_bytes;
    lock.unlock();
    queue_ready_.notify_one();
}

void OutputManager3D::worker_loop() noexcept {
    try {
        for (;;) {
            OutputJob3D job;
            {
                std::unique_lock lock(queue_mutex_);
                queue_ready_.wait(lock, [this] {
                    return !jobs_.empty() || worker_stopping_;
                });
                if (jobs_.empty() && worker_stopping_) break;
                const std::uint64_t job_bytes =
                    estimate_job_bytes(jobs_.front());
                job = std::move(jobs_.front());
                jobs_.pop_front();
                queued_snapshot_bytes_ =
                    job_bytes > queued_snapshot_bytes_
                    ? 0
                    : queued_snapshot_bytes_ - job_bytes;
            }
            queue_space_.notify_one();

            if (job.live_preview) {
                if (!job.preview_state) {
                    throw std::logic_error(
                        "live preview job has no sampled snapshot");
                }
                write_live_preview(*job.preview_state);
            }
            if (job.preview) {
                if (job.preview_state) {
                    write_preview(*job.preview_state);
                } else if (job.full_state) {
                    write_preview(*job.full_state);
                } else {
                    throw std::logic_error(
                        "preview output job has no snapshot");
                }
            }
            if (job.full) {
                if (!job.full_state) {
                    throw std::logic_error(
                        "full output job has no full snapshot");
                }
                write_full(*job.full_state);
            }
            if (job.checkpoint_base) {
                if (!job.full_state) {
                    throw std::logic_error(
                        "base checkpoint job has no full snapshot");
                }
                write_checkpoint(*job.full_state);
            } else if (job.checkpoint_journal) {
                write_checkpoint(*job.checkpoint_journal);
            }
        }
    } catch (...) {
        {
            std::lock_guard lock(queue_mutex_);
            worker_error_ = std::current_exception();
            worker_stopping_ = true;
            jobs_.clear();
            queued_snapshot_bytes_ = 0;
        }
        worker_failed_.store(true, std::memory_order_release);
        queue_space_.notify_all();
        queue_ready_.notify_all();
    }
}

void OutputManager3D::finish_worker() {
    if (!worker_.joinable()) return;
    {
        std::lock_guard lock(queue_mutex_);
        worker_stopping_ = true;
    }
    queue_ready_.notify_all();
    worker_.join();
    throw_worker_error();
}

bool OutputManager3D::viewer_attached() const {
    if (!config_.live_preview_when_attached) return false;
    const std::filesystem::path marker =
        run_directory_ / "control" / "viewer.attached";
    std::error_code error;
    const auto written =
        std::filesystem::last_write_time(marker, error);
    if (error) return false;
    return std::filesystem::file_time_type::clock::now() - written <=
           std::chrono::seconds(5);
}

void OutputManager3D::observe(const Simulation3D& simulation) {
    if (!config_.output_enabled) return;
    if (!resume_validated_) validate_resume_state(simulation);
    throw_worker_error();
    const double now = simulation.clock().time_hours;
    const bool full_due = due(now, next_full_, config_.full_every_hours);
    const bool preview_due = due(now, next_preview_, config_.preview_every_hours);
    const bool write_preview_now =
        (preview_due || full_due) &&
        !same_time(config_.output_async_enabled
                       ? last_scheduled_preview_time_ : last_preview_time_, now);
    const bool write_full_now =
        full_due &&
        !same_time(config_.output_async_enabled
                       ? last_scheduled_full_time_ : last_full_time_, now);
    const bool write_checkpoint_now =
        due(now, next_checkpoint_, config_.checkpoint_every_hours);
    if (write_preview_now || write_full_now || write_checkpoint_now) {
        OutputJob3D job;
        job.preview = write_preview_now;
        job.full = write_full_now;
        const bool base_checkpoint =
            write_checkpoint_now && checkpoint_base_due(now);
        job.checkpoint_base = base_checkpoint;
        if (write_full_now || base_checkpoint) {
            job.full_state = freeze(simulation);
        }
        if (write_preview_now && !job.full_state) {
            job.preview_state = freeze_preview(simulation);
        }
        if (write_checkpoint_now) {
            const std::filesystem::path checkpoint_path =
                run_directory_ / "checkpoints" /
                checkpoint_name(simulation.clock());
            if (base_checkpoint) {
                simulation.cells().reset_checkpoint_journal();
                checkpoint_base_time_ = now;
                checkpoint_delta_chain_ = 0;
                checkpoint_parent_path_ = checkpoint_path;
                checkpoint_parent_state_checksum_ =
                    job.full_state->state_checksum;
                checkpoint_parent_time_ = now;
                checkpoint_lineage_scheduled_ =
                    simulation.lineage().size();
            } else {
                job.checkpoint_journal =
                    freeze_checkpoint_journal(simulation,
                                              checkpoint_path);
                checkpoint_parent_path_ = checkpoint_path;
                checkpoint_parent_state_checksum_ =
                    job.checkpoint_journal->state_checksum;
                checkpoint_parent_time_ = now;
                checkpoint_delta_chain_ =
                    job.checkpoint_journal->chain_length;
                checkpoint_lineage_scheduled_ =
                    simulation.lineage().size();
            }
        }
        last_scheduled_state_time_ = now;
        if (write_preview_now) last_scheduled_preview_time_ = now;
        if (write_full_now) last_scheduled_full_time_ = now;
        if (config_.output_async_enabled) {
            enqueue(std::move(job));
        } else {
            if (job.preview) {
                if (job.preview_state) {
                    write_preview(*job.preview_state);
                } else {
                    write_preview(*job.full_state);
                }
            }
            if (job.full) write_full(*job.full_state);
            if (job.checkpoint_base) {
                write_checkpoint(*job.full_state);
            } else if (job.checkpoint_journal) {
                write_checkpoint(*job.checkpoint_journal);
            }
        }
    }
    const auto wall_now = std::chrono::steady_clock::now();
    const bool live_wall_due =
        last_live_preview_wall_.time_since_epoch().count() == 0 ||
        std::chrono::duration<double>(
            wall_now - last_live_preview_wall_).count() >=
            config_.live_preview_wall_interval_seconds;
    if (live_wall_due && viewer_attached() &&
        !write_preview_now && !write_full_now) {
        OutputJob3D live_job;
        live_job.live_preview = true;
        live_job.preview_state = freeze_preview(simulation);
        if (config_.output_async_enabled) {
            enqueue(std::move(live_job));
        } else {
            write_live_preview(*live_job.preview_state);
        }
        last_live_preview_wall_ = wall_now;
    }
    append_lineage(simulation);
    if (preview_due) next_preview_ = advance(next_preview_, now, config_.preview_every_hours);
    if (full_due) next_full_ = advance(next_full_, now, config_.full_every_hours);
    if (write_checkpoint_now) {
        next_checkpoint_ = advance(next_checkpoint_, now, config_.checkpoint_every_hours);
    }
}

void OutputManager3D::write_preview(const Simulation3D& simulation) {
    const std::vector<Slot> slots = stable_preview_sample(
        simulation.cells(), static_cast<std::size_t>(config_.preview_max_cells),
        config_.preview_seed);
    const SimulationSnapshotView3D snapshot{
        simulation.cells(), slots, simulation.lesion_index(), simulation.clock(),
        config_.display_radius};
    const std::string name = frame_name(preview_.size());
    const std::filesystem::path path = run_directory_ / "viz" / "preview" / name;
    require_unused_path(path);
    write_vtkhdf_points_atomic(path, snapshot, config_.vtkhdf_compression_level);
    write_vessels(simulation, name);
    write_vascular_metrics(simulation.clock().time_hours, simulation.stats(),
                           simulation.snapshot_vasculature());
    preview_.push_back({"viz/preview/" + name, simulation.clock().time_hours});
    last_preview_time_ = simulation.clock().time_hours;
    update_metadata();
}

void OutputManager3D::write_preview(
    const FrozenPreviewState3D& simulation) {
    std::vector<std::size_t> indices(simulation.cells.size());
    std::iota(indices.begin(), indices.end(), std::size_t{0});
    const FrozenCellSnapshotView3D snapshot{
        simulation.cells, simulation.cell_slots,
        simulation.cell_slot_count, indices, simulation.lesion_ids,
        simulation.clock, config_.display_radius};
    const std::string name = frame_name(preview_.size());
    const std::filesystem::path path =
        run_directory_ / "viz" / "preview" / name;
    require_unused_path(path);
    write_vtkhdf_points_atomic(
        path, snapshot, config_.vtkhdf_compression_level);
    const std::filesystem::path vessel_path =
        run_directory_ / "viz" / "vessels" / name;
    require_unused_path(vessel_path);
    write_vtkhdf_vessels_atomic(
        vessel_path, simulation.vasculature.nodes,
        config_.vtkhdf_compression_level);
    vessels_.push_back(
        {"viz/vessels/" + name, simulation.clock.time_hours});
    write_vascular_metrics(simulation.clock.time_hours,
                           simulation.stats,
                           simulation.vasculature);
    preview_.push_back(
        {"viz/preview/" + name, simulation.clock.time_hours});
    last_preview_time_ = simulation.clock.time_hours;
    update_metadata();
}

void OutputManager3D::write_preview(
    const FrozenSimulationState3D& simulation) {
    const std::vector<std::size_t> indices = stable_preview_sample(
        std::span<const CellInit>(simulation.cells),
        static_cast<std::size_t>(config_.preview_max_cells),
        config_.preview_seed);
    const FrozenCellSnapshotView3D snapshot{
        simulation.cells, simulation.cell_slots, simulation.cell_slot_count,
        indices, simulation.lesion_ids, simulation.clock,
        config_.display_radius};
    const std::string name = frame_name(preview_.size());
    const std::filesystem::path path = run_directory_ / "viz" / "preview" / name;
    require_unused_path(path);
    write_vtkhdf_points_atomic(path, snapshot, config_.vtkhdf_compression_level);
    write_vessels(simulation, name);
    write_vascular_metrics(simulation.clock.time_hours, simulation.stats,
                           simulation.vasculature);
    preview_.push_back({"viz/preview/" + name, simulation.clock.time_hours});
    last_preview_time_ = simulation.clock.time_hours;
    update_metadata();
}

void OutputManager3D::write_live_preview(
    const FrozenPreviewState3D& simulation) {
    std::vector<std::size_t> indices(simulation.cells.size());
    std::iota(indices.begin(), indices.end(), std::size_t{0});
    const FrozenCellSnapshotView3D snapshot{
        simulation.cells, simulation.cell_slots,
        simulation.cell_slot_count, indices, simulation.lesion_ids,
        simulation.clock, config_.display_radius};
    const std::filesystem::path cell_path =
        run_directory_ / "viz" / "live" / "current.vtkhdf";
    const std::filesystem::path vessel_path =
        run_directory_ / "viz" / "live" / "vessels.vtkhdf";
    write_vtkhdf_points_atomic(
        cell_path, snapshot, config_.vtkhdf_compression_level);
    write_vtkhdf_vessels_atomic(
        vessel_path, simulation.vasculature.nodes,
        config_.vtkhdf_compression_level);
    write_series_atomic(
        run_directory_ / "live.vtkhdf.series",
        {{"viz/live/current.vtkhdf",
          simulation.clock.time_hours}});
    write_series_atomic(
        run_directory_ / "live-vessels.vtkhdf.series",
        {{"viz/live/vessels.vtkhdf",
          simulation.clock.time_hours}});
}

void OutputManager3D::write_vessels(const Simulation3D& simulation,
                                    const std::string& frame) {
    const std::filesystem::path path = run_directory_ / "viz" / "vessels" / frame;
    require_unused_path(path);
    write_vtkhdf_vessels_atomic(path, simulation.vessel_nodes(),
                                config_.vtkhdf_compression_level);
    vessels_.push_back({"viz/vessels/" + frame, simulation.clock().time_hours});
}

void OutputManager3D::write_vessels(
    const FrozenSimulationState3D& simulation,
    const std::string& frame) {
    const std::filesystem::path path = run_directory_ / "viz" / "vessels" / frame;
    require_unused_path(path);
    write_vtkhdf_vessels_atomic(path, simulation.vasculature.nodes,
                                config_.vtkhdf_compression_level);
    vessels_.push_back({"viz/vessels/" + frame, simulation.clock.time_hours});
}

void OutputManager3D::write_full(const Simulation3D& simulation) {
    const std::vector<Slot> slots = simulation.cells().alive_slots();
    const SimulationSnapshotView3D snapshot{
        simulation.cells(), slots, simulation.lesion_index(), simulation.clock(),
        config_.display_radius};
    const std::string name = frame_name(full_.size());
    const std::filesystem::path path = run_directory_ / "viz" / "full" / name;
    require_unused_path(path);
    write_vtkhdf_points_atomic(path, snapshot, config_.vtkhdf_compression_level);
    full_.push_back({"viz/full/" + name, simulation.clock().time_hours});
    last_full_time_ = simulation.clock().time_hours;
    write_series_atomic(run_directory_ / "full.vtkhdf.series", full_);
    update_metadata();
}

void OutputManager3D::write_full(const FrozenSimulationState3D& simulation) {
    std::vector<std::size_t> indices(simulation.cells.size());
    std::iota(indices.begin(), indices.end(), std::size_t{0});
    const FrozenCellSnapshotView3D snapshot{
        simulation.cells, simulation.cell_slots, simulation.cell_slot_count,
        indices, simulation.lesion_ids, simulation.clock,
        config_.display_radius};
    const std::string name = frame_name(full_.size());
    const std::filesystem::path path = run_directory_ / "viz" / "full" / name;
    require_unused_path(path);
    write_vtkhdf_points_atomic(path, snapshot, config_.vtkhdf_compression_level);
    full_.push_back({"viz/full/" + name, simulation.clock.time_hours});
    last_full_time_ = simulation.clock.time_hours;
    write_series_atomic(run_directory_ / "full.vtkhdf.series", full_);
    update_metadata();
}

void OutputManager3D::write_checkpoint(const Simulation3D& simulation) {
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
    FrozenSimulationState3D state = freeze(simulation);
    write_checkpoint(state);
    simulation.cells().reset_checkpoint_journal();
#else
    (void)simulation;
    throw std::runtime_error("HDF5 checkpoint support is not built");
#endif
}

void OutputManager3D::write_checkpoint(
    const FrozenSimulationState3D& simulation) {
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
    const std::filesystem::path path = run_directory_ / "checkpoints" /
                                       checkpoint_name(simulation.clock);
    require_unused_path(path);
    const CheckpointSnapshotView3D snapshot{
        simulation.cells,
        simulation.cell_slots,
        simulation.cell_slot_count,
        simulation.free_slots,
        simulation.next_uid,
        simulation.clock,
        simulation.stats,
        simulation.lineage,
        simulation.vasculature,
        simulation.state_checksum};
    write_hdf5_checkpoint(path, snapshot, config_);
    last_checkpoint_was_base_ = true;
    last_checkpoint_path_ = path;
#else
    (void)simulation;
    throw std::runtime_error("HDF5 checkpoint support is not built");
#endif
}

void OutputManager3D::write_checkpoint(
    const FrozenCheckpointJournal3D& simulation) {
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
    require_unused_path(simulation.path);
    const CheckpointJournalSnapshotView3D snapshot{
        simulation.cells.mutations,
        simulation.cells.free_list_mutations,
        simulation.cells.slot_count,
        simulation.next_uid,
        simulation.clock,
        simulation.stats,
        simulation.lineage_tail,
        simulation.lineage_prefix_count,
        simulation.vasculature,
        simulation.state_checksum};
    write_hdf5_journal_delta_checkpoint(
        simulation.path, snapshot, simulation.parent_path,
        simulation.parent_state_checksum,
        simulation.parent_time_hours, simulation.chain_length,
        config_);
    last_checkpoint_was_base_ = false;
    last_checkpoint_path_ = simulation.path;
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
    collect_if_present(artifacts, run_directory_,
                       "metrics/vascular_latest.json.tmp", true);
    collect_frame_artifacts(artifacts, run_directory_, "viz/preview", preview_prefix);
    collect_frame_artifacts(artifacts, run_directory_, "viz/full", full_prefix);
    collect_frame_artifacts(artifacts, run_directory_, "viz/vessels", vessel_prefix);
    collect_checkpoint_artifacts(artifacts, run_directory_,
                                 simulation.clock().completed_events, now);

    const bool has_non_temporary_artifact = std::any_of(
        artifacts.begin(), artifacts.end(), [](const RecoveryArtifact& artifact) {
            return !artifact.safe_before_catalog_update;
        });
    const bool state_rollback = catalogs_trimmed || lineage_trimmed ||
                                has_non_temporary_artifact;
    if (state_rollback) {
        collect_if_present(artifacts, run_directory_, "metrics/final.json", true);
        collect_if_present(artifacts, run_directory_,
                           "metrics/vascular_latest.json", true);
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
    // Constructor state may have observed catalogs newer than the selected
    // checkpoint. Recovery trims those catalogs above, so asynchronous
    // de-duplication markers must roll back with them; otherwise the first
    // post-resume preview/full boundary is incorrectly treated as already
    // scheduled while checkpoint output alone is emitted.
    last_scheduled_preview_time_ = last_preview_time_;
    last_scheduled_full_time_ = last_full_time_;

    // A resume never re-emits the restored instant. Continue each periodic
    // schedule at the first regular boundary strictly after the restored time.
    next_preview_ = advance(0.0, now, config_.preview_every_hours);
    next_full_ = advance(0.0, now, config_.full_every_hours);
    next_checkpoint_ = advance(0.0, now, config_.checkpoint_every_hours);
    // A resumed run deliberately starts a new base at its next checkpoint.
    // This keeps schema-4/5/6/7 recovery compatible without retaining or
    // materializing the restored parent solely for output bookkeeping.
    checkpoint_parent_path_.clear();
    checkpoint_parent_state_checksum_ = 0;
    checkpoint_parent_time_ = -1.0;
    checkpoint_base_time_ = -1.0;
    checkpoint_delta_chain_ = 0;
    checkpoint_lineage_scheduled_ = simulation.lineage().size();
    simulation.cells().reset_checkpoint_journal();
    resume_time_ = now;
    resume_validated_ = true;
}

void OutputManager3D::append_lineage(const Simulation3D& simulation) {
    append_lineage(simulation.lineage());
}

void OutputManager3D::append_lineage(
    const std::vector<LineageEdge>& lineage) {
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
            << "  \"lesion_count\": "
            << simulation.lesion_index().lesions().size() << ",\n"
            << "  \"vessel_nodes\": " << simulation.vessel_nodes().alive_count() << ",\n"
            << "  \"angiogenesis_seed_attempts\": "
            << simulation.stats().angiogenesis_seed_attempts << ",\n"
            << "  \"angiogenesis_roots\": "
            << simulation.stats().angiogenesis_roots << ",\n"
            << "  \"angiogenesis_seed_rejections\": "
            << simulation.stats().angiogenesis_seed_rejections << ",\n"
            << "  \"cell_store_bytes\": " << simulation.cells().allocated_bytes() << ",\n"
            << "  \"grid_bytes\": " << simulation.grid().allocated_bytes() << ",\n"
            << "  \"density_index_bytes\": "
            << simulation.density().allocated_bytes() << ",\n"
            << "  \"lesion_index_bytes\": "
            << simulation.lesion_index().allocated_bytes() << ",\n"
            << "  \"state_checksum\": " << simulation.state_checksum() << "\n"
            << "}\n";
    });
}

void OutputManager3D::write_vascular_metrics(
    double time_hours,
    const SimulationStats3D& stats,
    const VasculatureState3D& vasculature) {
    std::array<std::uint64_t, 9> status_counts{};
    std::array<std::uint64_t, 3> role_counts{};
    for (const VesselTipInit3D& tip : vasculature.tips) {
        const auto status = static_cast<std::size_t>(tip.status);
        const auto role = static_cast<std::size_t>(tip.role);
        if (status < status_counts.size()) ++status_counts[status];
        if (role < role_counts.size()) ++role_counts[role];
    }
    static constexpr std::array<std::string_view, 9> status_names{
        "dormant", "active", "blocked", "merged", "reached_target",
        "max_length", "boundary_stop", "complete", "transiting"};
    write_text_atomic(run_directory_ / "metrics" / "vascular_latest.json",
                      [&](std::ostream& out) {
        const auto& aggregate = vasculature.process;
        out << std::setprecision(17)
            << "{\n"
            << "  \"time_hours\": " << time_hours << ",\n"
            << "  \"vessel_nodes\": " << vasculature.nodes.size() << ",\n"
            << "  \"vessel_tips\": " << vasculature.tips.size() << ",\n"
            << "  \"tip_roles\": {\"root\": " << role_counts[0]
            << ", \"inward\": " << role_counts[1]
            << ", \"outward\": " << role_counts[2] << "},\n"
            << "  \"tip_status\": {";
        for (std::size_t index = 0; index < status_names.size(); ++index) {
            if (index != 0) out << ", ";
            out << '\"' << status_names[index] << "\": "
                << status_counts[index];
        }
        out << "},\n"
            << "  \"seed_attempts\": " << stats.angiogenesis_seed_attempts
            << ",\n"
            << "  \"committed_roots\": " << stats.angiogenesis_roots
            << ",\n"
            << "  \"seed_rejections\": "
            << stats.angiogenesis_seed_rejections << ",\n"
            << "  \"aggregate_density_stress\": "
            << aggregate.current_density_stress << ",\n"
            << "  \"aggregate_rate_sites_per_30_days\": "
            << aggregate.current_rate_sites_per_30_days << ",\n"
            << "  \"next_seed_time_hours\": "
            << aggregate.next_seed_time_hours << ",\n"
            << "  \"influence\": {\"profile\": \""
            << config_.angiogenesis.influence_profile
            << "\", \"max_relief_fraction\": "
            << config_.angiogenesis.influence_max_relief_fraction
            << ", \"cutoff_from_vessel_voxel_centers\": "
            << config_.angiogenesis.influence_cutoff_radius_voxels
            << ", \"maximum_centerline_extent_voxels\": "
            << config_.angiogenesis.diameter_voxels * 0.5 +
                   config_.angiogenesis.influence_cutoff_radius_voxels
            << "},\n"
            << "  \"lesions\": [";
        for (std::size_t index = 0;
             index < vasculature.lesions.processes.size(); ++index) {
            const LesionAngiogenesisState3D& entry =
                vasculature.lesions.processes[index];
            const auto& process = entry.process;
            if (index != 0) out << ',';
            out << "\n    {\"lesion_id\": " << entry.lesion_id
                << ", \"eligible\": "
                << (process.eligible ? "true" : "false")
                << ", \"density_stress\": "
                << process.current_density_stress
                << ", \"rate_sites_per_30_days\": "
                << process.current_rate_sites_per_30_days
                << ", \"next_seed_time_hours\": "
                << process.next_seed_time_hours
                << ", \"attempts\": " << process.attempted_events
                << ", \"roots\": " << process.committed_roots
                << ", \"rejections\": " << process.rejected_events << '}';
        }
        if (!vasculature.lesions.processes.empty()) out << '\n';
        out << "  ]\n}\n";
    });
}

void OutputManager3D::finalize(const Simulation3D& simulation) {
    if (!config_.output_enabled || finalized_) return;
    observe(simulation);
    const bool advanced_since_resume =
        !resume_mode_ || (simulation.clock().time_hours > resume_time_ &&
                          !same_time(simulation.clock().time_hours, resume_time_));
    if (config_.output_async_enabled) {
        const double now = simulation.clock().time_hours;
        const bool final_preview =
            !same_time(last_scheduled_preview_time_, now) &&
            advanced_since_resume &&
            (config_.preview_every_hours > 0.0 ||
             config_.full_every_hours > 0.0);
        if (final_preview) {
            OutputJob3D job;
            job.preview = final_preview;
            job.preview_state = freeze_preview(simulation);
            enqueue(std::move(job));
            last_scheduled_state_time_ = now;
            last_scheduled_preview_time_ = now;
        }
        finish_worker();
    } else if (!same_time(last_preview_time_, simulation.clock().time_hours) &&
        advanced_since_resume &&
        (config_.preview_every_hours > 0.0 || config_.full_every_hours > 0.0)) {
        write_preview(simulation);
    }
    append_lineage(simulation);
    write_vascular_metrics(simulation.clock().time_hours, simulation.stats(),
                           simulation.snapshot_vasculature());
    write_metrics(simulation);
    update_metadata();
    finalized_ = true;
}

void OutputManager3D::checkpoint_now(
    const Simulation3D& simulation) {
    if (!config_.output_enabled ||
        !config_.output_on_demand_checkpoint) {
        return;
    }
    if (checkpoint_parent_time_ >= 0.0 &&
        same_time(checkpoint_parent_time_,
                  simulation.clock().time_hours)) {
        return;
    }
    finish_worker();
    FrozenSimulationState3D state = freeze(simulation);
    const std::filesystem::path path =
        run_directory_ / "checkpoints" /
        checkpoint_name(simulation.clock());
    require_unused_path(path);
    write_checkpoint(state);
    simulation.cells().reset_checkpoint_journal();
    checkpoint_parent_path_ = path;
    checkpoint_parent_state_checksum_ = state.state_checksum;
    checkpoint_parent_time_ = state.clock.time_hours;
    checkpoint_base_time_ = state.clock.time_hours;
    checkpoint_delta_chain_ = 0;
    checkpoint_lineage_scheduled_ = state.lineage.size();
}

}  // namespace atcg3d
