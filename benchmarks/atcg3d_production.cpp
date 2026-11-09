#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <sys/resource.h>
#include <vector>

#include "engine/simulation.hpp"
#include "io/checkpoint_hdf5.hpp"
#include "model/hybrid_model.hpp"
#include "model/shared_angiogenesis.hpp"
#include "model/shared_resource_environment.hpp"
#include "model/structured_pde_model.hpp"

namespace {
struct Options {
    std::filesystem::path config, checkpoint, resume;
    std::string model;
    double stop_hours{};
    int threads{1};
};

struct Result {
    std::uint64_t checksum{}, field_checksum{};
    double time{}, mass{}, active_mass{}, radius_99{};
    double abm_fraction{}, vascular_roots{}, vascular_length{};
};

bool same_time(double lhs, double rhs) {
    return std::abs(lhs - rhs) <= 1.0e-10 * std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

std::uint64_t peak_memory() {
    struct rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
        throw std::runtime_error("unable to read peak resident memory");
    }
#ifdef __APPLE__
    return static_cast<std::uint64_t>(usage.ru_maxrss);
#else
    return static_cast<std::uint64_t>(usage.ru_maxrss) * 1024ULL;
#endif
}

std::uint64_t checkpoint_size(const std::filesystem::path& path) {
    if (path.empty()) {
        return 0;
    }
    std::uint64_t bytes = std::filesystem::file_size(path);
    for (const auto* suffix : {".resource.bin", ".pde.bin"}) {
        const std::filesystem::path sidecar = path.string() + suffix;
        if (std::filesystem::exists(sidecar)) {
            bytes += std::filesystem::file_size(sidecar);
        }
    }
    return bytes;
}

void configure_observation(atcg3d::Model3DConfig& base, double interval) {
    base.output_enabled = true;
    base.preview_every_hours = interval;
    base.full_every_hours = interval;
    base.checkpoint_every_hours = interval;
    base.preview_keyframe_every_hours = interval;
    base.full_keyframe_every_hours = std::max(base.full_keyframe_every_hours, interval);
    base.checkpoint_base_every_hours = std::max(base.checkpoint_base_every_hours, interval);
}

double radius_99(const std::vector<double>& radial_mass, double total) {
    if (total <= 0.0) {
        return 0.0;
    }
    double cumulative = 0.0;
    for (std::size_t bin = 0; bin < radial_mass.size(); ++bin) {
        cumulative += radial_mass[bin];
        if (cumulative >= 0.99 * total) {
            return static_cast<double>(bin);
        }
    }
    throw std::logic_error("radial mass does not contain the requested quantile");
}

void add_radial_mass(std::vector<double>& radial_mass, double x, double y,
                     double z, bool thin, double amount) {
    const double radius = std::sqrt(x * x + y * y + (thin ? 0.0 : z * z));
    const auto bin = static_cast<std::size_t>(std::floor(radius));
    if (bin >= radial_mass.size()) {
        radial_mass.resize(bin + 1, 0.0);
    }
    radial_mass[bin] += amount;
}

double agent_radius_99(const atcg3d::Simulation3D& simulation) {
    std::vector<double> radial_mass;
    for (const auto slot : simulation.cells().alive_slots()) {
        const auto point = simulation.cells().anchor(slot);
        add_radial_mass(radial_mass, point.x + 0.5, point.y + 0.5, point.z + 0.5,
            simulation.config().thin_layer, 1.0);
    }
    return radius_99(radial_mass, static_cast<double>(simulation.cells().alive_count()));
}

double pde_radius_99(const atcg3d::structured_pde::StructuredPdeModel3D& simulation,
                    double total) {
    std::vector<double> radial_mass;
    for (std::size_t i = 0; i < simulation.voxel_count(); ++i) {
        double amount = 0.0;
        for (const auto stage : {atcg3d::structured_pde::StructuredStage3D::small,
                                atcg3d::structured_pde::StructuredStage3D::large}) {
            amount += simulation.r_normal(stage, i) + simulation.r_active(stage, i) +
                simulation.K(stage, i);
        }
        const auto center = simulation.coordinate(i);
        add_radial_mass(radial_mass, center[0], center[1], center[2],
            simulation.config().continuum.base.thin_layer, amount * simulation.voxel_measure());
    }
    return radius_99(radial_mass, total);
}

Result run_abm(const Options& options,
               const atcg3d::structured_pde::StructuredPdeConfig3D& rules) {
#ifndef ATCG3D_HAS_HDF5_CHECKPOINT
    if (!options.checkpoint.empty() || !options.resume.empty()) {
        throw std::invalid_argument("production ABM restart requires HDF5 checkpoints");
    }
#endif
    auto base = atcg3d::shared_rules::abm_config(rules);
    configure_observation(base, rules.continuum.time_step_hours);
    auto environment = std::make_unique<atcg3d::shared_rules::SharedResourceEnvironment3D>(rules);
    auto* resource = environment.get();
    std::unique_ptr<atcg3d::Simulation3D> simulation;
    if (options.resume.empty()) {
        simulation = std::make_unique<atcg3d::Simulation3D>(base, std::move(environment));
        simulation->initialize();
    } else {
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
        base.run_mode = "resume";
        base.resume_checkpoint = options.resume;
        const auto saved = atcg3d::read_hdf5_checkpoint(options.resume, base);
        resource->load_checkpoint(options.resume.string() + ".resource.bin", saved.state_checksum);
        simulation = std::make_unique<atcg3d::Simulation3D>(base, std::move(environment));
        simulation->restore(saved.cells, saved.next_uid, saved.clock, saved.stats, saved.lineage,
            saved.vasculature, saved.cell_slot_count, saved.cell_slots, saved.cell_free_slots);
#endif
    }
    // Stop only after the resource event at this macro boundary. Keeping the
    // original end time preserves migration proposal and biological horizons.
    simulation->run([&](const atcg3d::Simulation3D& current) {
        if (same_time(current.clock().time_hours, options.stop_hours) &&
            resource->next_refresh_time_hours() > options.stop_hours) {
            simulation->request_stop();
        }
    });
    if (!same_time(simulation->clock().time_hours, options.stop_hours)) {
        throw std::runtime_error("ABM did not reach the requested resource barrier");
    }
    if (!options.checkpoint.empty()) {
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
        atcg3d::write_hdf5_checkpoint(options.checkpoint, *simulation);
        resource->save_checkpoint(options.checkpoint.string() + ".resource.bin", simulation->state_checksum());
#endif
    }
    Result result;
    result.time = simulation->clock().time_hours;
    result.checksum = simulation->state_checksum();
    result.field_checksum = resource->field_checksum();
    result.mass = simulation->cells().alive_count();
    result.abm_fraction = result.mass > 0.0 ? 1.0 : 0.0;
    result.radius_99 = agent_radius_99(*simulation);
    for (const auto slot : simulation->cells().alive_slots()) {
        if (simulation->cells().flags(slot) & atcg3d::kMigrationActive) {
            result.active_mass += 1.0;
        }
    }
    const auto* vascular = resource->angiogenesis()->shared_diagnostics();
    result.vascular_roots = vascular->seeded_tips;
    result.vascular_length = vascular->centerline_growth;
    return result;
}

Result run_pde(const Options& options,
               const atcg3d::structured_pde::StructuredPdeConfig3D& rules) {
    atcg3d::structured_pde::StructuredPdeModel3D simulation(rules);
    if (options.resume.empty()) {
        auto resource = std::make_unique<atcg3d::shared_rules::SharedResourceEnvironment3D>(rules);
        atcg3d::Simulation3D seed(atcg3d::shared_rules::abm_config(rules), std::move(resource));
        seed.initialize();
        double maximum_work = 0.0, maximum_active_hours = 0.0;
        for (const auto slot : seed.cells().alive_slots()) {
            maximum_work = std::max(maximum_work, double(seed.cells().division_work_remaining(slot)));
            maximum_active_hours = std::max(maximum_active_hours,
                (seed.cells().flags(slot) & atcg3d::kMigrationActive)
                    ? seed.cells().migration_activation_end_time(slot) : 0.0);
        }
        std::cerr << "Initial maximum division work=" << maximum_work
                  << " active duration=" << maximum_active_hours << '\n';
        simulation.initialize_from_abm(seed);
    } else {
        simulation.load_checkpoint(options.resume);
    }
    while (simulation.time_hours() < options.stop_hours && simulation.step()) {
    }
    if (!same_time(simulation.time_hours(), options.stop_hours)) {
        throw std::runtime_error("PDE did not reach the requested macro boundary");
    }
    if (!options.checkpoint.empty()) {
        simulation.save_checkpoint(options.checkpoint);
    }
    const auto diagnostic = simulation.diagnostics();
    Result result;
    result.time = simulation.time_hours();
    result.checksum = simulation.state_checksum();
    result.field_checksum = simulation.angiogenesis()->checksum();
    result.mass = diagnostic.r_total + diagnostic.K_total;
    result.active_mass = diagnostic.r_active_total;
    result.radius_99 = pde_radius_99(simulation, result.mass);
    const auto* vascular = simulation.angiogenesis()->shared_diagnostics();
    result.vascular_roots = vascular->seeded_tips;
    result.vascular_length = vascular->centerline_growth;
    return result;
}

Result run_hybrid(const Options& options, atcg3d::hybrid::HybridConfig3D config) {
    config.rules.continuum.base.threads = options.threads;
    config.rules.continuum.output.enabled = false;
    atcg3d::hybrid::HybridModel3D simulation(config);
    if (options.resume.empty()) {
        simulation.initialize();
    } else {
        simulation.load_checkpoint(options.resume);
    }
    while (simulation.time_hours() < options.stop_hours && simulation.step()) {
    }
    if (!same_time(simulation.time_hours(), options.stop_hours)) {
        throw std::runtime_error("hybrid did not reach the requested macro boundary");
    }
    if (!options.checkpoint.empty()) {
        simulation.save_checkpoint(options.checkpoint);
    }
    const auto diagnostic = simulation.diagnostics();
    Result result;
    result.time = simulation.time_hours();
    result.checksum = simulation.state_checksum();
    result.field_checksum = simulation.pde().angiogenesis()->checksum();
    result.mass = diagnostic.total_mass;
    result.active_mass = diagnostic.active_mass;
    result.abm_fraction = result.mass > 0.0 ? diagnostic.abm_mass / result.mass : 0.0;
    result.radius_99 = radius_99(simulation.radial_mass(), result.mass);
    const auto* vascular = simulation.pde().angiogenesis()->shared_diagnostics();
    result.vascular_roots = vascular->seeded_tips;
    result.vascular_length = vascular->centerline_growth;
    return result;
}

Options parse(int argc, char** argv) {
    Options result;
    for (int i = 1; i < argc; ++i) {
        if (i + 1 == argc) {
            throw std::invalid_argument("production benchmark option needs a value");
        }
        const std::string key = argv[i], value = argv[++i];
        if (key == "--config") {
            result.config = value;
        } else if (key == "--model") {
            result.model = value;
        } else if (key == "--stop-hours") {
            result.stop_hours = std::stod(value);
        } else if (key == "--threads") {
            result.threads = std::stoi(value);
        } else if (key == "--checkpoint") {
            result.checkpoint = value;
        } else if (key == "--resume") {
            result.resume = value;
        } else {
            throw std::invalid_argument("unknown production benchmark option");
        }
    }
    if (result.config.empty() || result.threads < 1 ||
        !std::isfinite(result.stop_hours) || result.stop_hours <= 0.0 ||
        (result.model != "abm" && result.model != "pde" && result.model != "hybrid")) {
        throw std::invalid_argument("invalid production benchmark model, time or threads");
    }
    return result;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        const auto options = parse(argc, argv);
        auto hybrid = options.model == "hybrid"
            ? atcg3d::hybrid::HybridConfig3D::load(options.config)
            : atcg3d::hybrid::HybridConfig3D{};
        auto rules = options.model == "hybrid" ? hybrid.rules
            : atcg3d::structured_pde::StructuredPdeConfig3D::load(options.config);
        rules.continuum.base.threads = options.threads;
        rules.continuum.output.enabled = false;
        if (rules.continuum.angiogenesis.model != "shared_vegf_lattice_v2" ||
            rules.schema_version < 16 || options.stop_hours > rules.continuum.end_time_hours ||
            !same_time(options.stop_hours / rules.continuum.time_step_hours,
                std::round(options.stop_hours / rules.continuum.time_step_hours))) {
            throw std::invalid_argument("benchmark requires shared angiogenesis and a valid macro barrier");
        }
        const auto start = std::chrono::steady_clock::now();
        const auto result = options.model == "hybrid" ? run_hybrid(options, hybrid)
            : options.model == "abm" ? run_abm(options, rules) : run_pde(options, rules);
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        const auto shape = rules.continuum.grid.shape;
        const double half_width = 0.5 * std::min(shape[0], shape[1]) * rules.continuum.grid.spacing_voxels;
        std::cout << std::setprecision(17)
                  << "{\"model\":\"" << options.model << "\",\"threads\":" << options.threads
                  << ",\"grid_shape\":[" << shape[0] << ',' << shape[1] << ',' << shape[2] << ']'
                  << ",\"configured_hours\":" << rules.continuum.end_time_hours
                  << ",\"time_hours\":" << result.time
                  << ",\"wall_seconds\":" << elapsed << ",\"peak_resident_bytes\":" << peak_memory()
                  << ",\"checkpoint_bytes\":" << checkpoint_size(options.checkpoint)
                  << ",\"mass\":" << result.mass << ",\"active_mass\":" << result.active_mass
                  << ",\"abm_fraction\":" << result.abm_fraction << ",\"radius_99\":" << result.radius_99
                  << ",\"radius_99_half_width_ratio\":" << result.radius_99 / half_width
                  << ",\"vascular_roots\":" << result.vascular_roots << ",\"vascular_length\":" << result.vascular_length
                  << ",\"state_checksum\":" << result.checksum
                  << ",\"field_checksum\":" << result.field_checksum << "}\n";
    } catch (const std::exception& error) {
        std::cerr << "production benchmark: " << error.what() << '\n';
        return 1;
    }
}
