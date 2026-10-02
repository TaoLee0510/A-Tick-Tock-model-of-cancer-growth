#include <exception>
#include <filesystem>
#include <iostream>
#include "io/pde_vtkhdf_writer.hpp"
#include <memory>
#include <stdexcept>
#include <string>

#include "config/output_paths.hpp"

#include "config/structured_config.hpp"
#include "engine/simulation.hpp"
#include "io/structured_pde_output.hpp"
#include "model/structured_pde_model.hpp"
#include "model/shared_resource_environment.hpp"
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
#include "io/checkpoint_hdf5.hpp"
#endif

namespace {

void help(const char* executable) {
    std::cout
        << "Usage: " << executable
        << " --config PATH [--dry-run] [--output-root PATH] "
           "[--vtkhdf-fields] [--field-vtkhdf FILE]\n\n"
        << "Runs the ABM-aligned activation-clock structured PDE.\n";
}

std::unique_ptr<atcg3d::Simulation3D> source_abm(
    const atcg3d::structured_pde::StructuredPdeConfig3D& config) {
    atcg3d::Model3DConfig base = config.continuum.base;
    if (base.direction_guidance_model == "nutrient_gradient_shared_resource_v4") {
        if (config.continuum.initialization_mode != "base_model")
            throw std::invalid_argument("shared-resource ABM checkpoint import requires the shared-resource executable");
        auto environment = std::make_unique<atcg3d::shared_rules::SharedResourceEnvironment3D>(config);
        auto simulation = std::make_unique<atcg3d::Simulation3D>(
            atcg3d::shared_rules::abm_config(config), std::move(environment));
        simulation->initialize();
        return simulation;
    }
    base.output_enabled = false;
    base.control_enabled = false;
    if(base.angiogenesis.seed_process_model=="hypoxia_modulated_poisson_v2") {
        if(config.continuum.initialization_mode != "base_model") throw std::invalid_argument("hypoxic ABM checkpoint import requires the shared-resource executable");
        base.angiogenesis.enabled=false;
    }
    auto simulation = std::make_unique<atcg3d::Simulation3D>(base);
    if (config.continuum.initialization_mode == "base_model") {
        simulation->initialize();
        return simulation;
    }
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
    const atcg3d::CheckpointData3D state = atcg3d::read_hdf5_checkpoint(
        config.continuum.abm_checkpoint, base);
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
        std::filesystem::path output_root;
        bool dry_run = false;
        bool vtkhdf_fields=false;std::filesystem::path field_path;
        for (int index = 1; index < argc; ++index) {
            const std::string argument = argv[index];
            if (argument == "--help" || argument == "-h") {
                help(argv[0]);
                return 0;
            }
            if (argument == "--config" && index + 1 < argc) {
                config_path = argv[++index];
            } else if (argument == "--output-root" && index + 1 < argc) {
                output_root = argv[++index];
                if (output_root.empty()) {
                    throw std::invalid_argument("--output-root must not be empty");
                }
            } else if(argument=="--vtkhdf-fields") {vtkhdf_fields=true;
            } else if(argument=="--field-vtkhdf"&&index+1<argc) {field_path=argv[++index];
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

        auto config =
            atcg3d::structured_pde::StructuredPdeConfig3D::load(config_path);
        config.continuum.output.directory = atcg3d::resolve_output_directory(
            config.continuum.output.directory, output_root);
        if((vtkhdf_fields||!field_path.empty())&&!atcg3d::continuum::pde_vtkhdf_available())throw std::runtime_error("PDE VTK-HDF requires an HDF5-enabled build");
        config.continuum.output.vtkhdf_fields=vtkhdf_fields;
        if (dry_run) {
            std::cout << config.to_json();
            return 0;
        }
        atcg3d::structured_pde::StructuredPdeModel3D model(config);
        if (config.continuum.run_mode == "resume") {
            model.load_checkpoint(config.continuum.resume_checkpoint);
        } else {
            const auto simulation = source_abm(config);
            model.initialize_from_abm(*simulation);
        }

        std::cout << std::unitbuf
                  << "ATCG3D Structured PDE started\n"
                  << "profile=" << config.profile << '\n'
                  << "activation_rate_rule=normal_multiplier\n"
                  << "activation_multiplier="
                  << config.continuum.base.activated_r_normal_multiplier << '\n'
                  << "direction_guidance="
                  << (config.schema_version >= 5
                          ? "nutrient_gradient_only_v4"
                          : config.continuum.base.direction_guidance_model)
                  << '\n'
                  << (config.schema_version >= 5
                          ? "direction_nutrient_window_edge="
                          : "direction_density_window_edge=")
                  << (config.schema_version >= 5
                          ? config.migration.direction_nutrient_window_edge
                          : config.migration.direction_density_window_edge)
                  << '\n'
                  << "activation_stop="
                  << config.migration.activation_stop << '\n'
                  << "crowding_exchange="
                  << config.migration.crowding_exchange << '\n'
                  << "vessel_exclusion="
                  << (config.migration.vessel_exclusion ? "true" : "false")
                  << '\n'
                  << "resource_consumption="
                  << config.continuum.nutrient.consumption_model << '\n'
                  << "nutrient_solver="
                  << config.continuum.nutrient.solver << '\n'
                  << "nutrient_boundary_mode="
                  << config.continuum.nutrient.boundary_mode << '\n'
                  << "r_to_K_consumption_ratio="
                  << (config.continuum.nutrient.K_consumption_rate_per_hour > 0.0
                          ? config.continuum.nutrient.r_consumption_rate_per_hour /
                                config.continuum.nutrient.K_consumption_rate_per_hour
                          : 0.0)
                  << '\n'
                  << "grid=" << config.continuum.grid.shape[0] << 'x'
                  << config.continuum.grid.shape[1] << 'x'
                  << config.continuum.grid.shape[2] << '\n'
                  << "moving_directions=" << model.moving_direction_count() << '\n';

        atcg3d::structured_pde::StructuredPdeOutput3D output(
            config, model.time_hours());
        output.observe(model);
        while (model.step()) output.observe(model);
        output.checkpoint_now(model);
        output.finalize(model);

        if(!field_path.empty())atcg3d::continuum::write_pde_vtkhdf(field_path,model);
        const auto value = model.diagnostics();
        std::cout << "ATCG3D Structured PDE completed\n"
                  << "time_hours=" << model.time_hours() << '\n'
                  << "steps=" << model.step_count() << '\n'
                  << "nutrient_solves=" << model.nutrient_solve_count() << '\n'
                  << "r_total=" << value.r_total << '\n'
                  << "r_active_total=" << value.r_active_total << '\n'
                  << "active_fraction=" << value.active_fraction << '\n'
                  << "K_total=" << value.K_total << '\n'
                  << "assembled_r_consumption_per_hour="
                  << value.assembled_r_consumption_per_hour << '\n'
                  << "assembled_K_consumption_per_hour="
                  << value.assembled_K_consumption_per_hour << '\n'
                  << "occupied_volume=" << value.occupied_volume << '\n'
                  << "max_occupied_fraction="
                  << value.maximum_occupied_fraction << '\n'
                  << "r_mean_radius=" << value.r_mean_radius << '\n'
                  << "K_mean_radius=" << value.K_mean_radius << '\n'
                  << "r_radius_90=" << value.r_radius_90 << '\n'
                  << "r_radius_99=" << value.r_radius_99 << '\n'
                  << "tumour_volume=" << value.tumour_volume << '\n'
                  << "tumour_front_volume=" << value.tumour_front_volume << '\n'
                  << "tumour_mean_nutrient="
                  << value.tumour_mean_nutrient << '\n'
                  << "tumour_front_mean_nutrient="
                  << value.tumour_front_mean_nutrient << '\n'
                  << "state_checksum=" << model.state_checksum() << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "atcg3d_structured_pde: " << error.what() << '\n';
        return 1;
    }
}
