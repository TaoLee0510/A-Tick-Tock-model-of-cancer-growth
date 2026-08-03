#include "io/continuum_output.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace atcg3d::continuum {
namespace {

void require_unused(const std::filesystem::path& path) {
    if (std::filesystem::exists(path) ||
        std::filesystem::exists(path.string() + ".tmp")) {
        throw std::runtime_error("refusing to overwrite continuum output: " +
                                 path.string());
    }
}

void write_text_atomic(const std::filesystem::path& path,
                       const std::string& text) {
    const std::filesystem::path temporary = path.string() + ".tmp";
    std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
    if (!stream) throw std::runtime_error("unable to create " + temporary.string());
    stream << text;
    stream.flush();
    if (!stream) throw std::runtime_error("unable to finish " + temporary.string());
    stream.close();
    std::filesystem::rename(temporary, path);
}

}  // namespace

bool ContinuumOutput3D::same_time(double lhs, double rhs) noexcept {
    return std::abs(lhs - rhs) <=
        1.0e-10 * std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

bool ContinuumOutput3D::due(double now, double next, double interval) noexcept {
    return interval > 0.0 && (now > next || same_time(now, next));
}

double ContinuumOutput3D::advance(double next,
                                  double now,
                                  double interval) noexcept {
    if (!(interval > 0.0)) return next;
    do {
        next += interval;
    } while (now > next || same_time(now, next));
    return next;
}

ContinuumOutput3D::ContinuumOutput3D(
    const ContinuumModelConfig3D& config,
    double initial_time_hours)
    : config_(config), directory_(config.output.directory) {
    if (!config_.output.enabled) return;
    const bool resume = config_.run_mode == "resume";
    if (!resume && std::filesystem::exists(directory_) &&
        (!std::filesystem::is_directory(directory_) ||
         std::filesystem::directory_iterator(directory_) !=
             std::filesystem::directory_iterator{})) {
        throw std::runtime_error(
            "refusing to overwrite existing continuum run directory: " +
            directory_.string());
    }
    std::filesystem::create_directories(directory_ / "fields");
    std::filesystem::create_directories(directory_ / "profiles");
    std::filesystem::create_directories(directory_ / "checkpoints");
    std::filesystem::create_directories(directory_ / "config");

    const std::filesystem::path metrics_path = directory_ / "metrics.csv";
    const bool append = resume && std::filesystem::exists(metrics_path);
    metrics_.open(metrics_path,
                  std::ios::binary | (append ? std::ios::app : std::ios::trunc));
    if (!metrics_) throw std::runtime_error("unable to open continuum metrics");
    if (!append) {
        metrics_ << "time_hours,step_count,nutrient_solves,r_small,r_large,K_small,"
                    "K_large,r_total,K_total,occupied_volume,max_occupied_fraction,"
                    "mean_nutrient,max_nutrient,r_mean_nutrient,K_mean_nutrient,"
                    "r_mean_radius,K_mean_radius,vessel_volume,state_checksum\n";
    }
    if (resume) {
        next_metrics_ = advance(config_.start_time_hours, initial_time_hours,
                                config_.output.metrics_every_hours);
        next_field_ = advance(config_.start_time_hours, initial_time_hours,
                              config_.output.field_every_hours);
        next_profile_ = advance(config_.start_time_hours, initial_time_hours,
                                config_.output.radial_profile_every_hours);
        next_checkpoint_ = advance(config_.start_time_hours, initial_time_hours,
                                   config_.output.checkpoint_every_hours);
        last_checkpoint_time_ = initial_time_hours;
        for (const auto& entry : std::filesystem::directory_iterator(directory_ / "fields")) {
            if (entry.is_regular_file()) ++field_index_;
        }
        for (const auto& entry : std::filesystem::directory_iterator(directory_ / "profiles")) {
            if (entry.is_regular_file()) ++profile_index_;
        }
    } else {
        next_metrics_ = initial_time_hours;
        next_field_ = initial_time_hours;
        next_profile_ = initial_time_hours;
        next_checkpoint_ = initial_time_hours;
        if (!std::filesystem::exists(directory_ / "config" / "requested.yaml")) {
            std::filesystem::copy_file(
                config_.source_path, directory_ / "config" / "requested.yaml");
        }
        write_text_atomic(directory_ / "config" / "effective.json",
                          config_.to_json());
    }
}

std::filesystem::path ContinuumOutput3D::checkpoint_path(
    const std::filesystem::path& directory,
    const ContinuumModel3D& model) {
    std::ostringstream name;
    name << "checkpoint_step_" << std::setw(12) << std::setfill('0')
         << model.step_count() << "_time_" << std::hex << std::setw(16)
         << std::setfill('0') << std::bit_cast<std::uint64_t>(model.time_hours())
         << ".continuum.bin";
    return directory / "checkpoints" / name.str();
}

void ContinuumOutput3D::write_metrics(const ContinuumModel3D& model) {
    const ContinuumDiagnostics3D value = model.diagnostics();
    metrics_ << std::setprecision(17)
             << model.time_hours() << ',' << model.step_count() << ','
             << model.nutrient_solve_count();
    for (const double mass : value.population_mass) metrics_ << ',' << mass;
    metrics_ << ',' << value.type_mass[0] << ',' << value.type_mass[1]
             << ',' << value.occupied_volume
             << ',' << value.maximum_occupied_fraction
             << ',' << value.mean_nutrient << ',' << value.maximum_nutrient
             << ',' << value.mean_nutrient_by_type[0]
             << ',' << value.mean_nutrient_by_type[1]
             << ',' << value.mean_radius_by_type[0]
             << ',' << value.mean_radius_by_type[1]
             << ',' << value.vessel_volume
             << ',' << model.state_checksum() << '\n';
    metrics_.flush();
    if (!metrics_) throw std::runtime_error("unable to write continuum metrics");
}

void ContinuumOutput3D::write_field(const ContinuumModel3D& model) {
    std::ostringstream name;
    name << "field_" << std::setw(8) << std::setfill('0') << field_index_++
         << ".csv";
    const std::filesystem::path path = directory_ / "fields" / name.str();
    require_unused(path);
    const std::filesystem::path temporary = path.string() + ".tmp";
    std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
    if (!stream) throw std::runtime_error("unable to create continuum field output");
    stream << "time_hours,x,y,z,r_small,r_large,K_small,K_large,r_total,K_total,"
              "occupied_fraction,nutrient,vessel_fraction\n"
           << std::setprecision(10);
    for (std::size_t location = 0; location < model.voxel_count(); ++location) {
        const double total = model.population(PopulationField3D::r_small)[location] +
            model.population(PopulationField3D::r_large)[location] +
            model.population(PopulationField3D::K_small)[location] +
            model.population(PopulationField3D::K_large)[location];
        if (total <= 0.0 && model.nutrient()[location] <= 0.0 &&
            model.vessel_fraction()[location] <= 0.0) continue;
        const auto point = model.coordinate(location);
        const double r_total =
            model.population(PopulationField3D::r_small)[location] +
            model.population(PopulationField3D::r_large)[location];
        const double K_total =
            model.population(PopulationField3D::K_small)[location] +
            model.population(PopulationField3D::K_large)[location];
        stream << model.time_hours() << ',' << point[0] << ',' << point[1]
               << ',' << point[2];
        for (std::size_t field = 0; field < kPopulationFieldCount3D; ++field) {
            stream << ',' << model.population(
                static_cast<PopulationField3D>(field))[location];
        }
        stream << ',' << r_total << ',' << K_total
               << ',' << model.occupied_fraction(location)
               << ',' << model.nutrient()[location]
               << ',' << model.vessel_fraction()[location] << '\n';
    }
    stream.flush();
    if (!stream) throw std::runtime_error("unable to finish continuum field output");
    stream.close();
    std::filesystem::rename(temporary, path);
}

void ContinuumOutput3D::write_radial_profile(const ContinuumModel3D& model) {
    struct Bin {
        std::size_t voxels{};
        std::array<long double, 4> population{};
        long double nutrient{};
        long double vessel{};
    };
    double maximum_radius = 0.0;
    for (std::size_t location = 0; location < model.voxel_count(); ++location) {
        const auto point = model.coordinate(location);
        maximum_radius = std::max(maximum_radius,
            std::sqrt(point[0] * point[0] + point[1] * point[1] +
                (config_.base.thin_layer ? 0.0 : point[2] * point[2])));
    }
    const double width = config_.grid.spacing_voxels;
    std::vector<Bin> bins(static_cast<std::size_t>(std::floor(maximum_radius / width)) + 1);
    for (std::size_t location = 0; location < model.voxel_count(); ++location) {
        const auto point = model.coordinate(location);
        const double radius = std::sqrt(point[0] * point[0] + point[1] * point[1] +
            (config_.base.thin_layer ? 0.0 : point[2] * point[2]));
        Bin& bin = bins[static_cast<std::size_t>(std::floor(radius / width))];
        ++bin.voxels;
        for (std::size_t field = 0; field < 4; ++field) {
            bin.population[field] += model.population(
                static_cast<PopulationField3D>(field))[location];
        }
        bin.nutrient += model.nutrient()[location];
        bin.vessel += model.vessel_fraction()[location];
    }
    std::ostringstream name;
    name << "profile_" << std::setw(8) << std::setfill('0') << profile_index_++
         << ".csv";
    const std::filesystem::path path = directory_ / "profiles" / name.str();
    require_unused(path);
    const std::filesystem::path temporary = path.string() + ".tmp";
    std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
    if (!stream) throw std::runtime_error("unable to create continuum profile output");
    stream << "time_hours,radius_inner,radius_outer,voxels,r_small_mean,r_large_mean,"
              "K_small_mean,K_large_mean,nutrient_mean,vessel_fraction_mean\n"
           << std::setprecision(10);
    for (std::size_t index = 0; index < bins.size(); ++index) {
        const Bin& bin = bins[index];
        if (bin.voxels == 0) continue;
        stream << model.time_hours() << ',' << index * width << ','
               << (index + 1) * width << ',' << bin.voxels;
        for (const long double population : bin.population) {
            stream << ',' << static_cast<double>(population / bin.voxels);
        }
        stream << ',' << static_cast<double>(bin.nutrient / bin.voxels)
               << ',' << static_cast<double>(bin.vessel / bin.voxels) << '\n';
    }
    stream.flush();
    if (!stream) throw std::runtime_error("unable to finish continuum profile output");
    stream.close();
    std::filesystem::rename(temporary, path);
}

void ContinuumOutput3D::write_checkpoint(const ContinuumModel3D& model) {
    if (last_checkpoint_time_ >= 0.0 &&
        same_time(last_checkpoint_time_, model.time_hours())) return;
    model.save_checkpoint(checkpoint_path(directory_, model));
    last_checkpoint_time_ = model.time_hours();
}

void ContinuumOutput3D::observe(const ContinuumModel3D& model) {
    if (!config_.output.enabled) return;
    const double now = model.time_hours();
    if (due(now, next_metrics_, config_.output.metrics_every_hours)) {
        write_metrics(model);
        next_metrics_ = advance(next_metrics_, now,
                                config_.output.metrics_every_hours);
    }
    if (due(now, next_field_, config_.output.field_every_hours)) {
        write_field(model);
        next_field_ = advance(next_field_, now, config_.output.field_every_hours);
    }
    if (due(now, next_profile_, config_.output.radial_profile_every_hours)) {
        write_radial_profile(model);
        next_profile_ = advance(next_profile_, now,
                                config_.output.radial_profile_every_hours);
    }
    if (due(now, next_checkpoint_, config_.output.checkpoint_every_hours)) {
        write_checkpoint(model);
        next_checkpoint_ = advance(next_checkpoint_, now,
                                   config_.output.checkpoint_every_hours);
    }
}

void ContinuumOutput3D::checkpoint_now(const ContinuumModel3D& model) {
    if (config_.output.enabled) write_checkpoint(model);
}

void ContinuumOutput3D::finalize(const ContinuumModel3D& model) {
    if (!config_.output.enabled || finalized_) return;
    observe(model);
    const ContinuumDiagnostics3D value = model.diagnostics();
    std::ostringstream json;
    json << std::setprecision(17)
         << "{\n"
         << "  \"schema\": \"atcg3d.continuum.final\",\n"
         << "  \"time_hours\": " << model.time_hours() << ",\n"
         << "  \"step_count\": " << model.step_count() << ",\n"
         << "  \"r_total\": " << value.type_mass[0] << ",\n"
         << "  \"K_total\": " << value.type_mass[1] << ",\n"
         << "  \"maximum_occupied_fraction\": "
         << value.maximum_occupied_fraction << ",\n"
         << "  \"mean_nutrient\": " << value.mean_nutrient << ",\n"
         << "  \"state_checksum\": " << model.state_checksum() << "\n"
         << "}\n";
    write_text_atomic(directory_ / "final.json", json.str());
    finalized_ = true;
}

}  // namespace atcg3d::continuum
