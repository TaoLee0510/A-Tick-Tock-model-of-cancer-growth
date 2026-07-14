#include "io/output_manager.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

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

std::string checkpoint_name(std::uint64_t completed_events) {
    std::ostringstream value;
    value << "checkpoint_" << std::setw(16) << std::setfill('0') << completed_events << ".h5";
    return value.str();
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
    if (std::filesystem::exists(run_directory_ / "run.json")) {
        throw std::runtime_error("refusing to overwrite existing run directory: " +
                                 run_directory_.string());
    }
    std::filesystem::create_directories(run_directory_ / "metrics");
    std::filesystem::create_directories(run_directory_ / "checkpoints");
    std::filesystem::create_directories(run_directory_ / "lineage");
    std::filesystem::create_directories(run_directory_ / "viz" / "preview");
    std::filesystem::create_directories(run_directory_ / "viz" / "full");
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
    initialized_ = true;
}

void OutputManager3D::write_preview(const Simulation3D& simulation) {
    const std::vector<Slot> slots = stable_preview_sample(
        simulation.cells(), static_cast<std::size_t>(config_.preview_max_cells),
        config_.preview_seed);
    const SimulationSnapshotView3D snapshot{simulation.cells(), slots, simulation.clock()};
    const std::string name = frame_name(preview_.size());
    write_vtkhdf_points_atomic(run_directory_ / "viz" / "preview" / name, snapshot);
    preview_.push_back({"viz/preview/" + name, simulation.clock().time_hours});
    last_preview_time_ = simulation.clock().time_hours;
    write_series_atomic(run_directory_ / "preview.vtkhdf.series", preview_);
    update_metadata();
}

void OutputManager3D::write_full(const Simulation3D& simulation) {
    const std::vector<Slot> slots = simulation.cells().alive_slots();
    const SimulationSnapshotView3D snapshot{simulation.cells(), slots, simulation.clock()};
    const std::string name = frame_name(full_.size());
    write_vtkhdf_points_atomic(run_directory_ / "viz" / "full" / name, snapshot);
    full_.push_back({"viz/full/" + name, simulation.clock().time_hours});
    last_full_time_ = simulation.clock().time_hours;
    write_series_atomic(run_directory_ / "full.vtkhdf.series", full_);
    update_metadata();
}

void OutputManager3D::write_checkpoint(const Simulation3D& simulation) {
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
    write_hdf5_checkpoint(run_directory_ / "checkpoints" /
                              checkpoint_name(simulation.clock().completed_events),
                          simulation);
#else
    (void)simulation;
    throw std::runtime_error("HDF5 checkpoint support is not built");
#endif
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
    write_series_atomic(run_directory_ / "preview.vtkhdf.series", preview_);
    write_series_atomic(run_directory_ / "full.vtkhdf.series", full_);
    write_run_manifest_atomic(run_directory_, config_, preview_, full_,
                              config_.checkpoint_every_hours > 0.0);
}

void OutputManager3D::write_metrics(const Simulation3D& simulation) {
    write_text_atomic(run_directory_ / "metrics" / "final.json", [&](std::ostream& out) {
        out << std::setprecision(17)
            << "{\n"
            << "  \"time_hours\": " << simulation.clock().time_hours << ",\n"
            << "  \"completed_events\": " << simulation.clock().completed_events << ",\n"
            << "  \"alive_cells\": " << simulation.cells().alive_count() << ",\n"
            << "  \"cell_store_bytes\": " << simulation.cells().allocated_bytes() << ",\n"
            << "  \"grid_bytes\": " << simulation.grid().allocated_bytes() << ",\n"
            << "  \"state_checksum\": " << simulation.state_checksum() << "\n"
            << "}\n";
    });
}

void OutputManager3D::finalize(const Simulation3D& simulation) {
    if (!config_.output_enabled || finalized_) return;
    observe(simulation);
    if (!same_time(last_preview_time_, simulation.clock().time_hours) &&
        (config_.preview_every_hours > 0.0 || config_.full_every_hours > 0.0)) {
        write_preview(simulation);
    }
    append_lineage(simulation);
    write_metrics(simulation);
    update_metadata();
    finalized_ = true;
}

}  // namespace atcg3d
