#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <sys/resource.h>

#include "engine/simulation.hpp"
#include "model/shared_angiogenesis.hpp"
#include "model/shared_resource_environment.hpp"
#include "model/structured_pde_model.hpp"

namespace {
std::uint64_t peak_memory() {
    struct rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0)
        throw std::runtime_error("unable to read peak resident memory");
#ifdef __APPLE__
    return static_cast<std::uint64_t>(usage.ru_maxrss);
#else
    return static_cast<std::uint64_t>(usage.ru_maxrss) * 1024ULL;
#endif
}

struct Result {
    std::uint64_t checksum{}, field_checksum{};
    double mass{}, active_mass{}, roots{}, centerline{};
};

Result run_abm(const atcg3d::structured_pde::StructuredPdeConfig3D& config) {
    auto environment = std::make_unique<atcg3d::shared_rules::SharedResourceEnvironment3D>(config);
    auto* resource = environment.get();
    atcg3d::Simulation3D simulation(atcg3d::shared_rules::abm_config(config), std::move(environment));
    simulation.initialize();
    simulation.run();
    Result result;
    result.checksum = simulation.state_checksum();
    result.field_checksum = resource->field_checksum();
    result.mass = simulation.cells().alive_count();
    for (const auto slot : simulation.cells().alive_slots())
        if (simulation.cells().flags(slot) & atcg3d::kMigrationActive)
            result.active_mass += 1.0;
    result.roots = resource->angiogenesis()->shared_diagnostics()->seeded_tips;
    result.centerline = resource->angiogenesis()->shared_diagnostics()->centerline_growth;
    return result;
}

Result run_pde(const atcg3d::structured_pde::StructuredPdeConfig3D& config) {
    atcg3d::structured_pde::StructuredPdeModel3D simulation(config);
    {
        // Use the same resource-limited initial work and activation clocks as
        // the standalone shared ABM, then release its derived resource state.
        auto environment = std::make_unique<atcg3d::shared_rules::SharedResourceEnvironment3D>(config);
        atcg3d::Simulation3D seed(atcg3d::shared_rules::abm_config(config), std::move(environment));
        seed.initialize();
        simulation.initialize_from_abm(seed);
    }
    const auto start = std::chrono::steady_clock::now();
    while (simulation.step()) {
        std::cerr << "PDE time=" << simulation.time_hours() << " wall_seconds="
                  << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()
                  << '\n';
    }
    const auto diagnostic = simulation.diagnostics();
    Result result;
    result.checksum = simulation.state_checksum();
    result.field_checksum = simulation.angiogenesis()->checksum();
    result.mass = diagnostic.r_total + diagnostic.K_total;
    result.active_mass = diagnostic.r_active_total;
    result.roots = simulation.angiogenesis()->shared_diagnostics()->seeded_tips;
    result.centerline = simulation.angiogenesis()->shared_diagnostics()->centerline_growth;
    return result;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        std::filesystem::path path;
        std::string model = "abm";
        int threads = 1;
        for (int i = 1; i < argc; ++i) {
            if (i + 1 == argc) throw std::invalid_argument("benchmark option needs a value");
            const std::string option = argv[i], value = argv[++i];
            if (option == "--config") path = value;
            else if (option == "--model") model = value;
            else if (option == "--threads") threads = std::stoi(value);
            else throw std::invalid_argument("unknown 3D benchmark option");
        }
        if (threads < 1 || (model != "abm" && model != "pde"))
            throw std::invalid_argument("invalid 3D benchmark model or threads");
        auto config = atcg3d::structured_pde::StructuredPdeConfig3D::load(path);
        config.continuum.base.threads = threads;
        config.continuum.output.enabled = false;
        if (config.continuum.base.thin_layer || config.schema_version < 16 ||
            config.continuum.angiogenesis.model != "shared_vegf_lattice_v2" ||
            config.sector_mean_model != "prepared_prefix_fft_v2")
            throw std::invalid_argument("benchmark requires 3D prefix sectors and shared angiogenesis");
        config.validate();
        const auto start = std::chrono::steady_clock::now();
        const auto result = model == "abm" ? run_abm(config) : run_pde(config);
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        const auto shape = config.continuum.grid.shape;
        std::cout << std::setprecision(17)
                  << "{\"model\":\"" << model << "\",\"threads\":" << threads
                  << ",\"grid_shape\":[" << shape[0] << ',' << shape[1] << ',' << shape[2] << ']'
                  << ",\"hours\":" << config.continuum.end_time_hours
                  << ",\"wall_seconds\":" << elapsed << ",\"peak_resident_bytes\":" << peak_memory()
                  << ",\"mass\":" << result.mass << ",\"active_mass\":" << result.active_mass
                  << ",\"vascular_roots\":" << result.roots << ",\"vascular_length\":" << result.centerline
                  << ",\"state_checksum\":" << result.checksum
                  << ",\"field_checksum\":" << result.field_checksum << "}\n";
    } catch (const std::exception& error) {
        std::cerr << "3D benchmark: " << error.what() << '\n';
        return 1;
    }
}
