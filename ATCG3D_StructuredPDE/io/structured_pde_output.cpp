#include "io/structured_pde_output.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace atcg3d::structured_pde {
namespace {

bool same_time(double lhs, double rhs) noexcept {
    return std::abs(lhs - rhs) <=
        1.0e-10 * std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

void write_text_atomic(const std::filesystem::path& path,
                       const std::string& text) {
    const auto temporary = std::filesystem::path(path.string() + ".tmp");
    std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
    if (!stream) throw std::runtime_error("unable to create " + temporary.string());
    stream << text;
    stream.flush();
    if (!stream) throw std::runtime_error("unable to finish " + temporary.string());
    stream.close();
    std::filesystem::rename(temporary, path);
}

}  // namespace

bool StructuredPdeOutput3D::due(
    double now, double next, double interval) noexcept {
    return interval > 0.0 && (now > next || same_time(now, next));
}

double StructuredPdeOutput3D::advance(
    double next, double now, double interval) noexcept {
    if (!(interval > 0.0)) return next;
    do {
        next += interval;
    } while (now > next || same_time(now, next));
    return next;
}

StructuredPdeOutput3D::StructuredPdeOutput3D(
    const StructuredPdeConfig3D& config,
    double initial_time_hours)
    : config_(config), directory_(config.continuum.output.directory) {
    const auto& output = config_.continuum.output;
    if (!output.enabled) return;
    const bool resume = config_.continuum.run_mode == "resume";
    if (!resume && std::filesystem::exists(directory_) &&
        (!std::filesystem::is_directory(directory_) ||
         std::filesystem::directory_iterator(directory_) !=
             std::filesystem::directory_iterator{})) {
        throw std::runtime_error(
            "refusing to overwrite structured PDE run directory: " +
            directory_.string());
    }
    std::filesystem::create_directories(directory_ / "fields");
    std::filesystem::create_directories(directory_ / "checkpoints");
    std::filesystem::create_directories(directory_ / "config");
    const auto metrics_path = directory_ / "metrics.csv";
    const bool append = resume && std::filesystem::exists(metrics_path);
    metrics_.open(metrics_path,
                  std::ios::binary | (append ? std::ios::app : std::ios::trunc));
    if (!metrics_) throw std::runtime_error("unable to open structured metrics");
    if (!append) {
        metrics_ << "time_hours,step_count,nutrient_solves,r_normal_small,"
                    "r_normal_large,r_active_small,r_active_large,r_total,"
                    "r_active_total,active_fraction,K_small,K_large,K_total,"
                    "assembled_r_consumption_per_hour,"
                    "assembled_K_consumption_per_hour,occupied_volume,"
                    "max_occupied_fraction,mean_nutrient,"
                    "max_nutrient,r_mean_nutrient,K_mean_nutrient,r_mean_radius,"
                    "K_mean_radius,r_radius_50,r_radius_90,r_radius_99,"
                    "vessel_volume";
        if (config_.schema_version >= 6) {
            metrics_ << ",tumour_volume,tumour_front_volume,"
                        "tumour_mean_nutrient,tumour_front_mean_nutrient";
        }
        metrics_ << ",state_checksum\n";
    }
    if (resume) {
        next_metrics_ = advance(config_.continuum.start_time_hours,
            initial_time_hours, output.metrics_every_hours);
        next_field_ = advance(config_.continuum.start_time_hours,
            initial_time_hours, output.field_every_hours);
        next_checkpoint_ = advance(config_.continuum.start_time_hours,
            initial_time_hours, output.checkpoint_every_hours);
        last_checkpoint_time_ = initial_time_hours;
        for (const auto& entry :
             std::filesystem::directory_iterator(directory_ / "fields")) {
            if (entry.is_regular_file()) ++field_index_;
        }
    } else {
        next_metrics_ = initial_time_hours;
        next_field_ = initial_time_hours;
        next_checkpoint_ = initial_time_hours;
    }
    // A continuation may intentionally write into a fresh output directory.
    // Archive the exact resume request there as well, while preserving any
    // configuration already present when appending to an existing run.
    if (!std::filesystem::exists(directory_ / "config" / "requested.yaml")) {
        std::filesystem::copy_file(
            config_.source_path, directory_ / "config" / "requested.yaml");
    }
    if (!std::filesystem::exists(
            directory_ / "config" / "continuum_requested.yaml")) {
        std::filesystem::copy_file(
            config_.continuum_config_path,
            directory_ / "config" / "continuum_requested.yaml");
    }
    if (!std::filesystem::exists(directory_ / "config" / "effective.json")) {
        write_text_atomic(directory_ / "config" / "effective.json",
                          config_.to_json());
    }
}

void StructuredPdeOutput3D::write_metrics(const StructuredPdeModel3D& model) {
    const auto value = model.diagnostics();
    metrics_ << std::setprecision(17)
             << model.time_hours() << ',' << model.step_count() << ','
             << model.nutrient_solve_count() << ','
             << value.r_normal_mass[0] << ',' << value.r_normal_mass[1] << ','
             << value.r_active_mass[0] << ',' << value.r_active_mass[1] << ','
             << value.r_total << ',' << value.r_active_total << ','
             << value.active_fraction << ',' << value.K_mass[0] << ','
             << value.K_mass[1] << ',' << value.K_total << ','
             << value.assembled_r_consumption_per_hour << ','
             << value.assembled_K_consumption_per_hour << ','
             << value.occupied_volume << ',' << value.maximum_occupied_fraction
             << ',' << value.mean_nutrient << ',' << value.maximum_nutrient
             << ',' << value.r_mean_nutrient << ',' << value.K_mean_nutrient
             << ',' << value.r_mean_radius << ',' << value.K_mean_radius
             << ',' << value.r_radius_50 << ',' << value.r_radius_90
             << ',' << value.r_radius_99 << ',' << value.vessel_volume;
    if (config_.schema_version >= 6) {
        metrics_ << ',' << value.tumour_volume
                 << ',' << value.tumour_front_volume
                 << ',' << value.tumour_mean_nutrient
                 << ',' << value.tumour_front_mean_nutrient;
    }
    metrics_ << ',' << model.state_checksum() << '\n';
    metrics_.flush();
    if (!metrics_) throw std::runtime_error("unable to write structured metrics");
}

void StructuredPdeOutput3D::write_field(const StructuredPdeModel3D& model) {
    std::ostringstream name;
    name << "field_" << std::setw(8) << std::setfill('0') << field_index_++
         << ".csv";
    const auto path = directory_ / "fields" / name.str();
    if (std::filesystem::exists(path) ||
        std::filesystem::exists(path.string() + ".tmp")) {
        throw std::runtime_error("refusing to overwrite structured field");
    }
    const auto temporary = std::filesystem::path(path.string() + ".tmp");
    std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
    if (!stream) throw std::runtime_error("unable to create structured field");
    stream << "time_hours,x,y,z,r_normal_small,r_active_small,r_normal_large,"
              "r_active_large,K_small,K_large,r_total,K_total,active_fraction,"
              "activation_density_small,activation_density_large,"
              "occupied_fraction,nutrient,vessel_fraction";
    if (config_.schema_version >= 6) stream << ",tumour_mask";
    stream << '\n' << std::setprecision(10);
    const int nx = config_.continuum.grid.shape[0];
    const int ny = config_.continuum.grid.shape[1];
    const int stride = config_.continuum.output.field_stride;
    for (std::size_t location = 0; location < model.voxel_count(); ++location) {
        const int x = static_cast<int>(location % static_cast<std::size_t>(nx));
        const std::size_t yz = location / static_cast<std::size_t>(nx);
        const int y = static_cast<int>(yz % static_cast<std::size_t>(ny));
        const int z = static_cast<int>(yz / static_cast<std::size_t>(ny));
        if (x % stride != 0 || y % stride != 0 || z % stride != 0) continue;
        const double normal_small =
            model.r_normal(StructuredStage3D::small, location);
        const double active_small =
            model.r_active(StructuredStage3D::small, location);
        const double normal_large =
            model.r_normal(StructuredStage3D::large, location);
        const double active_large =
            model.r_active(StructuredStage3D::large, location);
        const double K_small = model.K(StructuredStage3D::small, location);
        const double K_large = model.K(StructuredStage3D::large, location);
        const double r_total = normal_small + active_small + normal_large + active_large;
        const double K_total = K_small + K_large;
        if (r_total + K_total <= 0.0 && model.nutrient()[location] <= 0.0 &&
            model.vessel_fraction()[location] <= 0.0) continue;
        const auto point = model.coordinate(location);
        const double active = active_small + active_large;
        stream << model.time_hours() << ',' << point[0] << ',' << point[1]
               << ',' << point[2] << ',' << normal_small << ',' << active_small
               << ',' << normal_large << ',' << active_large << ',' << K_small
               << ',' << K_large << ',' << r_total << ',' << K_total << ','
               << (r_total > 0.0 ? active / r_total : 0.0) << ','
               << model.activation_density(StructuredStage3D::small, location)
               << ','
               << model.activation_density(StructuredStage3D::large, location)
               << ',' << model.occupied_fraction(location) << ','
               << model.nutrient()[location] << ','
               << model.vessel_fraction()[location];
        if (config_.schema_version >= 6) {
            stream << ',' << static_cast<int>(model.tumour_mask()[location]);
        }
        stream << '\n';
    }
    stream.flush();
    if (!stream) throw std::runtime_error("unable to finish structured field");
    stream.close();
    std::filesystem::rename(temporary, path);
}

void StructuredPdeOutput3D::write_checkpoint(
    const StructuredPdeModel3D& model) {
    if (last_checkpoint_time_ >= 0.0 &&
        same_time(last_checkpoint_time_, model.time_hours())) return;
    std::ostringstream name;
    name << "checkpoint_step_" << std::setw(12) << std::setfill('0')
         << model.step_count() << "_time_" << std::hex << std::setw(16)
         << std::setfill('0') << std::bit_cast<std::uint64_t>(model.time_hours())
         << ".structured.bin";
    model.save_checkpoint(directory_ / "checkpoints" / name.str());
    last_checkpoint_time_ = model.time_hours();
}

void StructuredPdeOutput3D::observe(const StructuredPdeModel3D& model) {
    const auto& output = config_.continuum.output;
    if (!output.enabled) return;
    const double now = model.time_hours();
    if (due(now, next_metrics_, output.metrics_every_hours)) {
        write_metrics(model);
        next_metrics_ = advance(next_metrics_, now, output.metrics_every_hours);
    }
    if (due(now, next_field_, output.field_every_hours)) {
        write_field(model);
        next_field_ = advance(next_field_, now, output.field_every_hours);
    }
    if (due(now, next_checkpoint_, output.checkpoint_every_hours)) {
        write_checkpoint(model);
        next_checkpoint_ = advance(
            next_checkpoint_, now, output.checkpoint_every_hours);
    }
}

void StructuredPdeOutput3D::checkpoint_now(
    const StructuredPdeModel3D& model) {
    if (config_.continuum.output.enabled) write_checkpoint(model);
}

void StructuredPdeOutput3D::finalize(const StructuredPdeModel3D& model) {
    if (!config_.continuum.output.enabled || finalized_) return;
    observe(model);
    const auto value = model.diagnostics();
    std::ostringstream json;
    json << std::setprecision(17)
         << "{\n  \"schema\": \"atcg3d.structured_pde.final\",\n"
         << "  \"time_hours\": " << model.time_hours() << ",\n"
         << "  \"step_count\": " << model.step_count() << ",\n"
         << "  \"r_total\": " << value.r_total << ",\n"
         << "  \"r_active_total\": " << value.r_active_total << ",\n"
         << "  \"active_fraction\": " << value.active_fraction << ",\n"
         << "  \"K_total\": " << value.K_total << ",\n"
         << "  \"assembled_r_consumption_per_hour\": "
         << value.assembled_r_consumption_per_hour << ",\n"
         << "  \"assembled_K_consumption_per_hour\": "
         << value.assembled_K_consumption_per_hour << ",\n"
         << "  \"r_radius_90\": " << value.r_radius_90 << ",\n"
         << "  \"r_radius_99\": " << value.r_radius_99 << ",\n"
         << "  \"maximum_occupied_fraction\": "
         << value.maximum_occupied_fraction << ",\n"
         << "  \"state_checksum\": " << model.state_checksum() << "\n}\n";
    write_text_atomic(directory_ / "final.json", json.str());
    finalized_ = true;
}

}  // namespace atcg3d::structured_pde
