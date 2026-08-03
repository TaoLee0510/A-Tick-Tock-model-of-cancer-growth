#include <exception>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include "config/continuum_config.hpp"
#include "engine/simulation.hpp"
#include "io/continuum_output.hpp"
#include "model/continuum_model.hpp"
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
#include "io/checkpoint_hdf5.hpp"
#endif

namespace {

void help(const char* executable) {
    std::cout
        << "Usage: " << executable << " --config PATH [--dry-run]\n\n"
        << "Runs the four-population ATCG3D continuum reduction coupled to "
           "an effective-nutrient PDE.\n";
}

std::unique_ptr<atcg3d::Simulation3D> source_abm(
    const atcg3d::continuum::ContinuumModelConfig3D& config) {
    atcg3d::Model3DConfig base = config.base;
    base.output_enabled = false;
    base.control_enabled = false;
    auto simulation = std::make_unique<atcg3d::Simulation3D>(base);
    if (config.initialization_mode == "base_model") {
        simulation->initialize();
        return simulation;
    }
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
    const atcg3d::CheckpointData3D state =
        atcg3d::read_hdf5_checkpoint(config.abm_checkpoint, base);
    simulation->restore(
        state.cells, state.next_uid, state.clock, state.stats, state.lineage,
        state.vasculature, state.cell_slot_count, state.cell_slots,
        state.cell_free_slots);
    if (simulation->state_checksum() != state.state_checksum) {
        throw std::runtime_error("imported ABM checkpoint checksum mismatch");
    }
    return simulation;
#else
    throw std::runtime_error(
        "ABM checkpoint initialization requires ATCG3D_ENABLE_HDF5_CHECKPOINT=ON");
#endif
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::filesystem::path config_path;
        bool dry_run = false;
        for (int index = 1; index < argc; ++index) {
            const std::string argument = argv[index];
            if (argument == "--help" || argument == "-h") {
                help(argv[0]);
                return 0;
            }
            if (argument == "--config" && index + 1 < argc) {
                config_path = argv[++index];
            } else if (argument == "--dry-run") {
                dry_run = true;
            } else {
                throw std::invalid_argument(
                    "unknown or incomplete argument: " + argument);
            }
        }
        if (config_path.empty()) {
            throw std::invalid_argument("--config PATH is required");
        }

        atcg3d::continuum::ContinuumModelConfig3D config =
            atcg3d::continuum::ContinuumModelConfig3D::load(config_path);
        if (dry_run) {
            std::cout << config.to_json();
            return 0;
        }

        atcg3d::continuum::ContinuumModel3D model(config);
        if (config.run_mode == "resume") {
            model.load_checkpoint(config.resume_checkpoint);
        } else {
            const auto simulation = source_abm(config);
            model.initialize_from_abm(*simulation);
        }

        std::cout << std::unitbuf
                  << "ATCG3D Continuum started\n"
                  << "profile=" << config.profile << '\n'
                  << "initialization=" << config.initialization_mode << '\n'
                  << "grid=" << config.grid.shape[0] << 'x'
                  << config.grid.shape[1] << 'x' << config.grid.shape[2] << '\n'
                  << "spacing_voxels=" << config.grid.spacing_voxels << '\n';

        atcg3d::continuum::ContinuumOutput3D output(config, model.time_hours());
        output.observe(model);
        while (model.step()) output.observe(model);
        output.checkpoint_now(model);
        output.finalize(model);

        const auto value = model.diagnostics();
        std::cout << "ATCG3D Continuum completed\n"
                  << "time_hours=" << model.time_hours() << '\n'
                  << "steps=" << model.step_count() << '\n'
                  << "nutrient_solves=" << model.nutrient_solve_count() << '\n'
                  << "r_total=" << value.type_mass[0] << '\n'
                  << "K_total=" << value.type_mass[1] << '\n'
                  << "occupied_volume=" << value.occupied_volume << '\n'
                  << "max_occupied_fraction="
                  << value.maximum_occupied_fraction << '\n'
                  << "mean_nutrient=" << value.mean_nutrient << '\n'
                  << "r_mean_nutrient=" << value.mean_nutrient_by_type[0] << '\n'
                  << "K_mean_nutrient=" << value.mean_nutrient_by_type[1] << '\n'
                  << "r_mean_radius=" << value.mean_radius_by_type[0] << '\n'
                  << "K_mean_radius=" << value.mean_radius_by_type[1] << '\n'
                  << "state_checksum=" << model.state_checksum() << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "atcg3d_continuum: " << error.what() << '\n';
        return 1;
    }
}
