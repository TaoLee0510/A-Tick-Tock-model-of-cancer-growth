#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

#include "config/model_config.hpp"
#include "engine/simulation.hpp"
#include "io/output_manager.hpp"
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
#include "io/checkpoint_hdf5.hpp"
#endif

namespace {

void print_help(const char* executable) {
    std::cout << "Usage: " << executable
              << " --config PATH [--dry-run]\n"
              << "\n"
              << "Runs the sparse, event-driven ATCG 3D model. The 3D target never emits PNG.\n"
              << "All model/run parameters, including resume input, come from the YAML file.\n";
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
                throw std::invalid_argument("unknown or incomplete argument: " + argument);
            }
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
        atcg3d::OutputManager3D output(config);
        simulation.run([&output](const atcg3d::Simulation3D& current) { output.observe(current); });
        output.finalize(simulation);
        const auto& stats = simulation.stats();
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
                  << "migration_attempts=" << stats.migration_attempts << '\n'
                  << "migration_commits=" << stats.migration_commits << '\n'
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
                  << "checksum=" << simulation.state_checksum() << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "atcg3d: " << error.what() << '\n';
        return 1;
    }
}
