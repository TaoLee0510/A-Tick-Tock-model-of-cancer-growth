#include "io/nutrient_output.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "engine/simulation.hpp"

namespace atcg3d::nutrient {

bool NutrientOutput3D::same_time(double lhs, double rhs) noexcept {
    return std::abs(lhs - rhs) <=
        1.0e-10 * std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

bool NutrientOutput3D::due(double now, double next, double interval) noexcept {
    return interval > 0.0 && (now > next || same_time(now, next));
}

double NutrientOutput3D::advance(double next,
                                 double now,
                                 double interval) noexcept {
    if (!(interval > 0.0)) return next;
    do {
        next += interval;
    } while (now > next || same_time(now, next));
    return next;
}

NutrientOutput3D::NutrientOutput3D(
    const NutrientModelConfig3D& config,
    NutrientEnvironment3D& environment,
    double resume_time_hours)
    : config_(config),
      environment_(environment),
      directory_(config.base.output_directory / "nutrient") {
    if (!config_.base.output_enabled) return;
    std::filesystem::create_directories(directory_ / "field");
    const std::filesystem::path metrics_path = directory_ / "metrics.csv";
    const bool append = resume_time_hours >= 0.0 &&
        std::filesystem::exists(metrics_path);
    metrics_.open(metrics_path,
                  std::ios::binary |
                      (append ? std::ios::app : std::ios::trunc));
    if (!metrics_) {
        throw std::runtime_error("unable to open nutrient metrics output");
    }
    if (!append) {
        metrics_ << "time_hours,completed_events,refresh_count,blocks,active_voxels,"
                    "source_voxels,consuming_voxels,minimum,maximum,mean,"
                    "last_max_update,field_checksum,base_checksum\n";
    }
    if (resume_time_hours >= 0.0) {
        next_metrics_ = advance(
            0.0, resume_time_hours, config_.nutrient.metrics_every_hours);
        next_field_snapshot_ = advance(
            0.0, resume_time_hours,
            config_.nutrient.field_snapshot_every_hours);
        next_checkpoint_ = advance(
            0.0, resume_time_hours, config_.base.checkpoint_every_hours);
        last_sidecar_time_ = resume_time_hours;
    }
}

std::filesystem::path NutrientOutput3D::base_checkpoint_path(
    const std::filesystem::path& output_directory,
    const SimulationClock3D& clock) {
    std::ostringstream name;
    name << "checkpoint_" << std::setw(16) << std::setfill('0')
         << clock.completed_events << "_time_" << std::hex << std::setw(16)
         << std::setfill('0') << std::bit_cast<std::uint64_t>(clock.time_hours)
         << ".h5";
    return output_directory / "checkpoints" / name.str();
}

std::filesystem::path NutrientOutput3D::sidecar_path(
    const std::filesystem::path& base_checkpoint) {
    std::filesystem::path result = base_checkpoint;
    result.replace_extension(".nutrient.bin");
    return result;
}

void NutrientOutput3D::write_metrics(const Simulation3D& simulation) {
    const NutrientFieldDiagnostics3D& field = environment_.diagnostics();
    metrics_ << std::setprecision(17)
             << simulation.clock().time_hours << ','
             << simulation.clock().completed_events << ','
             << environment_.refresh_count() << ','
             << field.block_count << ','
             << field.active_voxel_count << ','
             << field.perfused_source_voxels << ','
             << field.consuming_voxels << ','
             << field.minimum << ',' << field.maximum << ',' << field.mean << ','
             << field.last_max_update << ','
             << environment_.field_checksum() << ','
             << simulation.state_checksum() << '\n';
    metrics_.flush();
    if (!metrics_) throw std::runtime_error("unable to write nutrient metrics");
}

void NutrientOutput3D::write_field_snapshot(
    const Simulation3D& simulation) {
    std::ostringstream name;
    name << "field_" << std::setw(8) << std::setfill('0')
         << field_snapshot_index_++ << ".csv";
    const std::filesystem::path final_path = directory_ / "field" / name.str();
    const std::filesystem::path temporary = final_path.string() + ".tmp";
    std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
    if (!stream) throw std::runtime_error("unable to create nutrient field snapshot");
    stream << "time_hours,x,y,z,nutrient\n" << std::setprecision(9);
    for (const NutrientVoxelSample3D& sample : environment_.nonzero_voxels()) {
        stream << simulation.clock().time_hours << ','
               << sample.site.x << ',' << sample.site.y << ',' << sample.site.z
               << ',' << sample.value << '\n';
    }
    stream.flush();
    if (!stream) throw std::runtime_error("unable to finish nutrient field snapshot");
    stream.close();
    std::filesystem::rename(temporary, final_path);
}

void NutrientOutput3D::observe(const Simulation3D& simulation) {
    if (!config_.base.output_enabled) return;
    const double now = simulation.clock().time_hours;
    if (due(now, next_metrics_, config_.nutrient.metrics_every_hours)) {
        write_metrics(simulation);
        next_metrics_ = advance(
            next_metrics_, now, config_.nutrient.metrics_every_hours);
    }
    if (due(now, next_field_snapshot_,
            config_.nutrient.field_snapshot_every_hours)) {
        write_field_snapshot(simulation);
        next_field_snapshot_ = advance(
            next_field_snapshot_, now,
            config_.nutrient.field_snapshot_every_hours);
    }
}

void NutrientOutput3D::write_sidecar(const Simulation3D& simulation) {
    if (last_sidecar_time_ >= 0.0 &&
        same_time(last_sidecar_time_, simulation.clock().time_hours)) {
        return;
    }
    const std::filesystem::path base = base_checkpoint_path(
        config_.base.output_directory, simulation.clock());
    environment_.save_checkpoint(
        sidecar_path(base), simulation.state_checksum(),
        simulation.clock().time_hours,
        simulation.clock().completed_events);
    last_sidecar_time_ = simulation.clock().time_hours;
}

void NutrientOutput3D::checkpoint_if_due(
    const Simulation3D& simulation) {
    if (!config_.base.output_enabled ||
        !(config_.base.checkpoint_every_hours > 0.0)) return;
    const double now = simulation.clock().time_hours;
    if (!due(now, next_checkpoint_, config_.base.checkpoint_every_hours)) return;
    write_sidecar(simulation);
    next_checkpoint_ = advance(
        next_checkpoint_, now, config_.base.checkpoint_every_hours);
}

void NutrientOutput3D::checkpoint_now(const Simulation3D& simulation) {
    if (!config_.base.output_enabled ||
        !config_.base.output_on_demand_checkpoint) return;
    write_sidecar(simulation);
}

}  // namespace atcg3d::nutrient
