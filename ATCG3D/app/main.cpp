#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "config/model_config.hpp"
#include "app/run_controller.hpp"
#include "engine/simulation.hpp"
#include "io/output_manager.hpp"
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
#include "io/checkpoint_hdf5.hpp"
#endif

namespace {

void print_help(const char* executable) {
    std::cout << "Usage: " << executable
              << " --config PATH [--dry-run]\n"
              << "       " << executable
              << " --print-config-schema\n"
              << "\n"
              << "Runs the sparse, event-driven ATCG 3D model. The 3D target never emits PNG.\n"
              << "All model/run parameters, including resume input, come from the YAML file.\n";
}

void print_config_schema() {
    const atcg3d::Model3DConfig defaults;
    std::cout
        << "{\n"
        << "  \"schema\": \"atcg3d.config-ui-schema\",\n"
        << "  \"schema_version\": 1,\n"
        << "  \"model_schema_name\": \"atcg3d.model_config\",\n"
        << "  \"model_schema_version\": 3,\n"
        << "  \"parameter_transport\": \"strict_yaml_only\",\n"
        << "  \"ui_groups\": [\"run\", \"space\", \"direction\", "
           "\"migration\", \"density\", \"biology\", \"stage\", "
           "\"division\", \"initial\", \"simulation\", \"parallel\", "
           "\"scheduler\", \"angiogenesis\", \"output\", \"control\"],\n"
        << "  \"defaults\": " << defaults.to_json()
        << "}\n";
}

void persist_run_config(
    const std::filesystem::path& source,
    const atcg3d::Model3DConfig& config) {
    const std::filesystem::path directory =
        config.output_directory / "config";
    std::filesystem::create_directories(directory);
    const std::filesystem::path requested =
        directory / "requested.yaml";
    if (!std::filesystem::exists(requested)) {
        std::filesystem::copy_file(source, requested);
    }
    const std::filesystem::path effective =
        directory / "effective.json";
    const std::filesystem::path temporary =
        effective.string() + ".tmp";
    {
        std::ofstream stream(
            temporary, std::ios::binary | std::ios::trunc);
        if (!stream) {
            throw std::runtime_error(
                "unable to create effective run config");
        }
        stream << config.to_json();
        stream.flush();
        if (!stream) {
            throw std::runtime_error(
                "unable to finish effective run config");
        }
    }
    std::filesystem::rename(temporary, effective);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::filesystem::path config_path;
        bool dry_run = false;
        bool print_schema = false;
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
            } else if (argument == "--print-config-schema") {
                print_schema = true;
            } else {
                throw std::invalid_argument("unknown or incomplete argument: " + argument);
            }
        }

        if (print_schema) {
            if (!config_path.empty() || dry_run) {
                throw std::invalid_argument(
                    "--print-config-schema cannot be combined with a run");
            }
            print_config_schema();
            return 0;
        }
        if (config_path.empty()) {
            throw std::invalid_argument("--config PATH is required");
        }

        atcg3d::Model3DConfig config = atcg3d::Model3DConfig::load(config_path);
        if (dry_run) {
            std::cout << config.to_json();
            return 0;
        }

        atcg3d::Simulation3D simulation(config);
        if (config.run_mode == "resume") {
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
            const atcg3d::CheckpointData3D checkpoint =
                atcg3d::read_hdf5_checkpoint(config.resume_checkpoint, config);
            simulation.restore(checkpoint.cells, checkpoint.next_uid, checkpoint.clock,
                               checkpoint.stats, checkpoint.lineage,
                               checkpoint.vasculature, checkpoint.cell_slot_count,
                               checkpoint.cell_slots, checkpoint.cell_free_slots);
            if (simulation.state_checksum() != checkpoint.state_checksum) {
                throw std::runtime_error(
                    "restored checkpoint state checksum does not match stored checksum");
            }
#else
            throw std::runtime_error(
                "run.mode=resume requires a build configured with "
                "ATCG3D_ENABLE_HDF5_CHECKPOINT=ON");
#endif
        }
        std::cout << std::unitbuf
                  << "ATCG3D started\n"
                  << "profile=" << config.profile << '\n'
                  << "threads=" << config.threads << '\n'
                  << "scheduler_backend=" << config.scheduler_backend << '\n'
                  << "proposal_window_hours=" << config.proposal_window_hours << '\n'
                  << "async_output=" << (config.output_async_enabled ? 1 : 0) << '\n';
        atcg3d::OutputManager3D output(config);
        persist_run_config(config_path, config);
        atcg3d::RunController3D controller(config);
        try {
        bool checkpoint_stop_requested = false;
        bool stop_requested = false;
        simulation.run([&](const atcg3d::Simulation3D& current) {
            output.observe(current);
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
        }
        if (stop_requested) {
            controller.mark_stopped(simulation);
        } else {
            controller.mark_completed(simulation);
        }
        const auto& stats = simulation.stats();
        const auto& proposal_diagnostics =
            simulation.proposal_window_diagnostics();
        const double cell_store_bytes_per_slot =
            simulation.cells().slot_count() == 0
                ? 0.0
                : static_cast<double>(simulation.cells().allocated_bytes()) /
                      static_cast<double>(simulation.cells().slot_count());
        std::cout << "ATCG3D completed\n"
                  << "time_hours=" << simulation.clock().time_hours << '\n'
                  << "completed_events=" << simulation.clock().completed_events << '\n'
                  << "alive_cells=" << simulation.cells().alive_count() << '\n'
                  << "chunks=" << simulation.grid().chunk_count() << '\n'
                  << "cell_store_bytes=" << simulation.cells().allocated_bytes() << '\n'
                  << "cell_store_bytes_per_slot=" << cell_store_bytes_per_slot << '\n'
                  << "cell_store_logical_bytes_per_slot="
                  << atcg3d::CellStore3D::logical_bytes_per_slot() << '\n'
                  << "grid_bytes=" << simulation.grid().allocated_bytes() << '\n'
                  << "density_index_bytes="
                  << simulation.density().allocated_bytes() << '\n'
                  << "migration_attempts=" << stats.migration_attempts << '\n'
                  << "migration_commits=" << stats.migration_commits << '\n'
                  << "migration_swap_waits="
                  << stats.migration_swap_waits << '\n'
                  << "migration_swap_attempts="
                  << stats.migration_swap_attempts << '\n'
                  << "migration_swap_commits="
                  << stats.migration_swap_commits << '\n'
                  << "migration_swap_rejections="
                  << stats.migration_swap_rejections << '\n'
                  << "divisions=" << stats.divisions << '\n'
                  << "deaths=" << stats.deaths << '\n'
                  << "angiogenesis_seed_attempts="
                  << stats.angiogenesis_seed_attempts << '\n'
                  << "angiogenesis_roots=" << stats.angiogenesis_roots << '\n'
                  << "vessel_growth_commits="
                  << stats.vessel_growth_commits << '\n'
                  << "vascular_displacements="
                  << stats.vascular_displacements << '\n'
                  << "vessel_nodes="
                  << simulation.vessel_nodes().alive_count() << '\n'
                  << "vessel_occupied_voxels="
                  << simulation.vessel_grid().occupied_voxel_count() << '\n'
                  << "proposal_windows=" << proposal_diagnostics.windows << '\n'
                  << "prefetched_migration_proposals="
                  << proposal_diagnostics.migration_proposals << '\n'
                  << "proposal_cache_hits="
                  << proposal_diagnostics.migration_cache_hits << '\n'
                  << "proposal_cache_invalidations="
                  << proposal_diagnostics.migration_cache_invalidations << '\n'
                  << "proposal_cache_misses="
                  << proposal_diagnostics.migration_cache_misses << '\n'
                  << "proposal_maximum_workers="
                  << proposal_diagnostics.maximum_workers << '\n'
                  << "checksum=" << simulation.state_checksum() << '\n';
        return 0;
        } catch (const std::exception& error) {
            controller.mark_failed(error.what());
            throw;
        }
    } catch (const std::exception& error) {
        std::cerr << "atcg3d: " << error.what() << '\n';
        return 1;
    }
}
