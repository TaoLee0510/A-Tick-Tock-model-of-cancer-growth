#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

#include "app/run_controller.hpp"
#include "config/nutrient_config.hpp"
#include "engine/simulation.hpp"
#include "field/nutrient_field.hpp"
#include "io/nutrient_output.hpp"
#include "io/output_manager.hpp"
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
#include "io/checkpoint_hdf5.hpp"
#endif

namespace {

void print_help(const char* executable) {
    std::cout
        << "Usage: " << executable << " --config PATH [--dry-run]\n\n"
        << "Runs the individual-cell ATCG3D PDMP coupled to the effective-"
           "nutrient reaction-diffusion field.\n";
}

void write_text_atomic(const std::filesystem::path& path,
                       const std::string& text) {
    std::filesystem::create_directories(path.parent_path());
    const std::filesystem::path temporary = path.string() + ".tmp";
    {
        std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
        if (!stream) throw std::runtime_error("unable to create " + temporary.string());
        stream << text;
        stream.flush();
        if (!stream) throw std::runtime_error("unable to finish " + temporary.string());
    }
    std::filesystem::rename(temporary, path);
}

void persist_config(const atcg3d::nutrient::NutrientModelConfig3D& config) {
    const std::filesystem::path directory =
        config.base.output_directory / "config";
    std::filesystem::create_directories(directory);
    const std::filesystem::path requested =
        directory / "nutrient_requested.yaml";
    if (!std::filesystem::exists(requested)) {
        std::filesystem::copy_file(config.source_path, requested);
    }
    const std::filesystem::path effective =
        directory / "nutrient_effective.json";
    if (!std::filesystem::exists(effective)) {
        write_text_atomic(effective, config.to_json());
    }
}

std::uint64_t combined_checksum(std::uint64_t base,
                                std::uint64_t nutrient) noexcept {
    nutrient += 0x9e3779b97f4a7c15ULL;
    nutrient = (nutrient ^ (nutrient >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    nutrient = (nutrient ^ (nutrient >> 27U)) * 0x94d049bb133111ebULL;
    nutrient ^= nutrient >> 31U;
    return base ^ nutrient;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::filesystem::path config_path;
        bool dry_run = false;
        for (int index = 1; index < argc; ++index) {
            const std::string argument = argv[index];
            if (argument == "--help" || argument == "-h") {
                print_help(argv[0]);
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

        atcg3d::nutrient::NutrientModelConfig3D config =
            atcg3d::nutrient::NutrientModelConfig3D::load(config_path);
        if (dry_run) {
            std::cout << config.to_json();
            return 0;
        }

        auto environment =
            std::make_unique<atcg3d::nutrient::NutrientEnvironment3D>(
                config.nutrient, config.base.thin_layer);
        auto* nutrient = environment.get();

#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
        std::optional<atcg3d::CheckpointData3D> restored;
        if (config.base.run_mode == "resume") {
            restored = atcg3d::read_hdf5_checkpoint(
                config.base.resume_checkpoint, config.base);
            nutrient->load_checkpoint(
                atcg3d::nutrient::NutrientOutput3D::sidecar_path(
                    config.base.resume_checkpoint),
                restored->state_checksum, restored->clock.time_hours,
                restored->clock.completed_events);
        }
#else
        if (config.base.run_mode == "resume") {
            throw std::runtime_error(
                "run.mode=resume requires ATCG3D_ENABLE_HDF5_CHECKPOINT=ON");
        }
#endif

        atcg3d::Simulation3D simulation(
            config.base, std::move(environment));
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
        if (restored) {
            simulation.restore(
                restored->cells, restored->next_uid, restored->clock,
                restored->stats, restored->lineage, restored->vasculature,
                restored->cell_slot_count, restored->cell_slots,
                restored->cell_free_slots);
            if (simulation.state_checksum() != restored->state_checksum) {
                throw std::runtime_error(
                    "restored base checkpoint checksum mismatch");
            }
        }
#endif

        std::cout << std::unitbuf
                  << "ATCG3D Nutrient started\n"
                  << "profile=" << config.profile << '\n'
                  << "base_profile=" << config.base.profile << '\n'
                  << "nutrient_solver=" << config.nutrient.solver << '\n'
                  << "nutrient_refresh_hours="
                  << config.nutrient.refresh_every_hours << '\n';

        atcg3d::OutputManager3D output(config.base);
        const double resume_time = config.base.run_mode == "resume"
            ? simulation.clock().time_hours : -1.0;
        atcg3d::nutrient::NutrientOutput3D nutrient_output(
            config, *nutrient, resume_time);
        persist_config(config);
        atcg3d::RunController3D controller(config.base);

        bool checkpoint_stop_requested = false;
        bool stop_requested = false;
        try {
            simulation.run([&](const atcg3d::Simulation3D& current) {
                output.observe(current);
                nutrient_output.observe(current);
                nutrient_output.checkpoint_if_due(current);
                const atcg3d::RunControlAction3D action =
                    controller.observe(current);
                if (action != atcg3d::RunControlAction3D::none) {
                    checkpoint_stop_requested =
                        action ==
                        atcg3d::RunControlAction3D::checkpoint_and_stop;
                    stop_requested = true;
                    simulation.request_stop();
                }
            });
            output.finalize(simulation);
            if (checkpoint_stop_requested) {
                output.checkpoint_now(simulation);
                nutrient_output.checkpoint_now(simulation);
            }
            if (stop_requested) controller.mark_stopped(simulation);
            else controller.mark_completed(simulation);
        } catch (const std::exception& error) {
            controller.mark_failed(error.what());
            throw;
        }

        const auto& field = nutrient->diagnostics();
        const std::uint64_t base_checksum = simulation.state_checksum();
        const std::uint64_t nutrient_checksum = nutrient->field_checksum();
        std::cout << "ATCG3D Nutrient completed\n"
                  << "time_hours=" << simulation.clock().time_hours << '\n'
                  << "completed_events="
                  << simulation.clock().completed_events << '\n'
                  << "alive_cells=" << simulation.cells().alive_count() << '\n'
                  << "nutrient_refreshes=" << nutrient->refresh_count() << '\n'
                  << "nutrient_blocks=" << field.block_count << '\n'
                  << "nutrient_active_voxels=" << field.active_voxel_count << '\n'
                  << "nutrient_min=" << field.minimum << '\n'
                  << "nutrient_max=" << field.maximum << '\n'
                  << "nutrient_mean=" << field.mean << '\n'
                  << "nutrient_bytes=" << nutrient->allocated_bytes() << '\n'
                  << "base_checksum=" << base_checksum << '\n'
                  << "nutrient_checksum=" << nutrient_checksum << '\n'
                  << "combined_checksum="
                  << combined_checksum(base_checksum, nutrient_checksum) << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "atcg3d_nutrient: " << error.what() << '\n';
        return 1;
    }
}
