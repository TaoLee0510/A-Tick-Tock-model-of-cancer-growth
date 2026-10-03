#include "config/output_paths.hpp"
#include "model/hybrid_model.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <numeric>
#include <cmath>
#include <sstream>
#include <limits>
#include "geometry/footprint.hpp"

namespace {
using atcg3d::hybrid::HybridConfig3D;
using atcg3d::hybrid::HybridModel3D;

double quantile(const std::vector<double>& radial, double total, double fraction) {
    if (!(total > 0.0))
        return 0.0;
    double sum = 0.0;
    for (std::size_t i = 0; i < radial.size(); ++i) {
        sum += radial[i];
        if (sum >= total * fraction)
            return double(i + 1);
    }
    return 0.0;
}

double boundary_mass(const HybridModel3D& simulation, const HybridConfig3D& config) {
    const auto& grid = config.rules.continuum.grid;
    const bool thin = config.rules.continuum.base.thin_layer;
    double mass = 0.0;
    if (config.mode != "all_pde")
        for (auto slot : simulation.abm().cells().alive_slots()) {
            const auto on_boundary = [&](atcg3d::Vec3i point) {
                const std::array<int, 3> coordinates{point.x, point.y, point.z};
                for (int axis = 0; axis < (thin ? 2 : 3); ++axis) {
                    const int index = int(std::floor(coordinates[axis] - grid.origin[axis]));
                    if (index <= 0 || index >= grid.shape[axis] - 1)
                        return true;
                }
                return false;
            };
            const auto anchor = simulation.abm().cells().anchor(slot);
            if (simulation.abm().cells().stage(slot) == atcg3d::CellStage::large) {
                for (auto point : atcg3d::large_footprint(anchor))
                    if ((!thin || point.z == anchor.z) && on_boundary(point))
                        mass += thin ? 0.25 : 0.125;
            } else if (on_boundary(anchor)) {
                mass += 1.0;
            }
        }
    if (config.mode != "all_abm")
        for (std::size_t i = 0; i < simulation.pde().voxel_count(); ++i) {
            const int x = int(i % grid.shape[0]);
            const int y = int((i / grid.shape[0]) % grid.shape[1]);
            const int z = int(i / (std::size_t(grid.shape[0]) * grid.shape[1]));
            if (x != 0 && x != grid.shape[0] - 1 && y != 0 && y != grid.shape[1] - 1 &&
                (thin || (z != 0 && z != grid.shape[2] - 1)))
                continue;
            for (auto stage : {atcg3d::structured_pde::StructuredStage3D::small,
                               atcg3d::structured_pde::StructuredStage3D::large})
                mass += simulation.pde().voxel_measure() *
                    (simulation.pde().r_normal(stage, i) + simulation.pde().r_active(stage, i) +
                     simulation.pde().K(stage, i));
        }
    return mass;
}

std::string summary(const HybridModel3D& simulation, const HybridConfig3D& config,
                    bool validation, bool sampling) {
    const auto d = simulation.diagnostics();
    std::ostringstream out;
    out << std::setprecision(17)
        << "{\"time_hours\":" << simulation.time_hours()
        << ",\"total_mass\":" << d.total_mass
        << ",\"abm_mass\":" << d.abm_mass << ",\"pde_mass\":" << d.pde_mass
        << ",\"active_mass\":" << d.active_mass
        << ",\"to_pde\":" << d.to_pde << ",\"to_abm\":" << d.to_abm
        << ",\"exchanges\":" << d.exchanges
        << ",\"state_checksum\":" << simulation.state_checksum();
    if (validation || config.schema_version >= 2 || sampling) {
        const auto radial = simulation.radial_mass();
        const double total = std::accumulate(radial.begin(), radial.end(), 0.0);
        const double r99 = quantile(radial, total, 0.99);
        out << ",\"r_mass\":" << d.r_mass << ",\"K_mass\":" << d.K_mass
            << ",\"r_K_ratio\":" << (d.K_mass > 0.0 ? d.r_mass / d.K_mass : 0.0)
            << ",\"active_fraction\":" << (d.r_mass > 0.0 ? d.active_mass / d.r_mass : 0.0)
            << ",\"r50\":" << quantile(radial, total, 0.5)
            << ",\"r90\":" << quantile(radial, total, 0.9)
            << ",\"r99\":" << r99 << ",\"radial_mass\":[";
        for (std::size_t i = 0; i < radial.size(); ++i) {
            if (i)
                out << ',';
            out << radial[i];
        }
        out << ']';
        if (sampling) {
            const auto& grid = config.rules.continuum.grid;
            const int dimensions = config.rules.continuum.base.thin_layer ? 2 : 3;
            double half_width = std::numeric_limits<double>::infinity();
            for (int axis = 0; axis < dimensions; ++axis)
                half_width = std::min(half_width, std::min(-grid.origin[axis],
                    grid.origin[axis] + grid.shape[axis] * grid.spacing_voxels));
            out << ",\"spatial_dimensions\":" << dimensions << ",\"grid_shape\":["
                << grid.shape[0] << ',' << grid.shape[1] << ',' << grid.shape[2] << ']'
                << ",\"domain_half_width\":" << half_width
                << ",\"r99_to_half_width\":" << r99 / half_width
                << ",\"boundary_mass\":" << boundary_mass(simulation, config);
        }
    }
    if (config.schema_version >= 4) {
        out << ",\"model\":\"" << config.model << "\",\"mode\":\"" << config.mode
            << "\",\"seed\":" << config.rules.continuum.base.seed
            << ",\"abm_fraction\":" << (d.total_mass > 0.0 ? d.abm_mass / d.total_mass : 0.0)
            << ",\"pde_fraction\":" << (d.total_mass > 0.0 ? d.pde_mass / d.total_mass : 0.0)
            << ",\"minimum_abm_fraction\":" << d.minimum_abm_fraction
            << ",\"minimum_active_fraction\":" << d.minimum_active_fraction
            << ",\"maximum_pde_fraction\":" << d.maximum_pde_fraction
            << ",\"front_pde_mass\":" << d.front_pde_mass
            << ",\"agent_migration_attempts\":" << simulation.abm().stats().migration_attempts
            << ",\"agent_migration_commits\":" << simulation.abm().stats().migration_commits
            << ",\"agent_swap_commits\":" << simulation.abm().stats().migration_swap_commits;
    }
    out << '}';
    return out.str();
}
}

int main(int argc, char **argv) {
    try {
        std::filesystem::path yaml, checkpoint, resume, report, output_root;
        std::string mode;
        int threads = 0;
        bool dry = false, no_output = false, validation_report = false;
        std::optional<std::uint64_t> seed;
        std::optional<double> step_hours, exchange_hours;
        double sample_interval = 0.0;
        std::vector<std::string> trace;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--help") {
                std::cout << "Usage: atcg3d_hybrid --config YAML [--mode "
                             "adaptive|all_abm|all_pde] [--threads N] "
                             "[--checkpoint FILE] "
                             "[--resume-checkpoint FILE] [--report JSON] "
                             "[--output-root PATH] [--dry-run] [--seed N] "
                             "[--step-hours H] [--exchange-hours H] "
                             "[--no-output] [--validation-report] [--sample-every-hours H]\n";
                return 0;
            }
            if (arg == "--dry-run") {
                dry = true;
                continue;
            }
            if (arg == "--no-output") { no_output = true; continue; }
            if (arg == "--validation-report") { validation_report = true; continue; }
            if (i + 1 == argc)
                throw std::invalid_argument("missing argument value");
            const std::string value = argv[++i];
            if (arg == "--config")
                yaml = value;
            else if (arg == "--mode")
                mode = value;
            else if (arg == "--threads")
                threads = std::stoi(value);
            else if (arg == "--checkpoint")
                checkpoint = value;
            else if (arg == "--resume-checkpoint")
                resume = value;
            else if (arg == "--output-root")
                output_root = value;
            else if (arg == "--report")
                report = value;
            else if (arg == "--seed") seed = std::stoull(value);
            else if (arg == "--step-hours") step_hours = std::stod(value);
            else if (arg == "--exchange-hours") exchange_hours = std::stod(value);
            else if (arg == "--sample-every-hours") sample_interval = std::stod(value);
            else
                throw std::invalid_argument("unknown argument: " + arg);
        }
        auto config = atcg3d::hybrid::HybridConfig3D::load(yaml);
        if (!mode.empty())
            config.mode = mode;
        if (threads)
            config.rules.continuum.base.threads = threads;
        if (seed) config.rules.continuum.base.seed = *seed;
        if (step_hours) {
            config.rules.continuum.time_step_hours = *step_hours;
            config.rules.continuum.nutrient.refresh_every_hours = *step_hours;
        }
        if (exchange_hours) config.exchange_every_hours = *exchange_hours;
        if (sample_interval != 0.0 &&
            (!std::isfinite(sample_interval) || sample_interval < 4.0 || sample_interval > 24.0))
            throw std::invalid_argument("sample interval must be between 4 and 24 hours");
        config.sample_every_hours = sample_interval;
        if (sample_interval > 0.0 &&
            std::abs(sample_interval / config.rules.continuum.time_step_hours -
                     std::round(sample_interval / config.rules.continuum.time_step_hours)) > 1e-10)
            throw std::invalid_argument("sampling interval must be a multiple of the hybrid time step");
        config.validate();
        if (dry) {
            std::cout << "{\"model\":\"" << config.model
                      << "\",\"fingerprint\":" << config.fingerprint();
            if (config.schema_version >= 4)
                std::cout << ",\"rules_config\":" << config.rules.to_json();
            std::cout << "}\n";
            return 0;
        }
        const auto directory = atcg3d::resolve_output_directory(
            config.output_directory, output_root);
        if (!no_output && resume.empty() && std::filesystem::exists(directory) &&
            !std::filesystem::is_empty(directory))
            throw std::runtime_error(
                "refusing to overwrite hybrid output directory");
        std::ofstream metrics;
        if (!no_output) {
            std::filesystem::create_directories(directory);
            metrics.open(directory / (resume.empty() ? "metrics.csv" : "metrics_resumed.csv"));
        }
        if (!no_output && !metrics)
            throw std::runtime_error("cannot write hybrid metrics");
        metrics << "time_hours,total_mass,abm_mass,pde_mass,active_mass,mean_"
                   "nutrient,to_pde,to_abm,state_checksum";
        if (config.schema_version >= 4)
            metrics << ",abm_fraction,pde_fraction,minimum_abm_fraction,minimum_active_fraction,maximum_pde_fraction";
        metrics << '\n' << std::setprecision(17);
        atcg3d::hybrid::HybridModel3D simulation(config);
        if (!resume.empty())
            simulation.load_checkpoint(resume);
        double next_sample = resume.empty() ? 0.0 : simulation.time_hours();
        const auto observe = [&] {
            if (sample_interval > 0.0 && simulation.time_hours() + 1e-9 >= next_sample) {
                trace.push_back(summary(simulation, config, validation_report, true));
                next_sample = simulation.time_hours() + sample_interval;
            }
            if (no_output) return;
            const auto d = simulation.diagnostics();
            metrics << simulation.time_hours() << ',' << d.total_mass << ','
                    << d.abm_mass << ',' << d.pde_mass << ',' << d.active_mass
                    << ',' << d.mean_nutrient << ',' << d.to_pde << ','
                    << d.to_abm << ',' << simulation.state_checksum();
            if (config.schema_version >= 4)
                metrics << ',' << (d.total_mass > 0.0 ? d.abm_mass / d.total_mass : 0.0)
                        << ',' << (d.total_mass > 0.0 ? d.pde_mass / d.total_mass : 0.0)
                        << ',' << d.minimum_abm_fraction << ',' << d.minimum_active_fraction
                        << ',' << d.maximum_pde_fraction;
            metrics << '\n';
        };
        if (resume.empty())
            simulation.initialize();
        observe();
        if (sample_interval > 0.0 && config.mode == "all_abm") {
            while (simulation.time_hours() < config.rules.continuum.end_time_hours) {
                simulation.run_abm_until(std::min(config.rules.continuum.end_time_hours, next_sample));
                observe();
            }
        } else {
            while (simulation.step())
                observe();
        }
        const auto endpoint = summary(simulation, config, validation_report, sample_interval > 0.0);
        if (sample_interval > 0.0 && (trace.empty() || trace.back() != endpoint))
            trace.push_back(endpoint);
        if (!checkpoint.empty())
            simulation.save_checkpoint(checkpoint);
        std::ofstream file;
        if (!report.empty()) {
            if (!report.parent_path().empty())
                std::filesystem::create_directories(report.parent_path());
            file.open(report);
            if (!file)
                throw std::runtime_error("cannot write hybrid report");
        }
        auto &out = report.empty() ? std::cout : file;
        out << endpoint.substr(0, endpoint.size() - 1);
        if (sample_interval > 0.0) {
            out << ",\"time_series\":[";
            for (std::size_t i = 0; i < trace.size(); ++i) {
                if (i)
                    out << ',';
                out << trace[i];
            }
            out << ']';
        }
        out << "}\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "atcg3d_hybrid: " << e.what() << '\n';
        return 1;
    }
}
