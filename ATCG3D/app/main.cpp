#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "config/model_config.hpp"
#include "engine/simulation.hpp"
#include "io/output_manager.hpp"
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
#include "io/checkpoint_hdf5.hpp"
#endif

namespace {

void print_help(const char* executable) {
    std::cout << "Usage: " << executable
              << " [--config PATH] [--set key=value] [--resume CHECKPOINT] [--dry-run]\n"
              << "\n"
              << "Runs the sparse, event-driven ATCG 3D model. The 3D target never emits PNG.\n"
              << "--set may be repeated and is applied after the configuration file.\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::filesystem::path config_path = "configs/atcg3d_legacy_like_v1.cfg";
        std::vector<std::string> overrides;
        std::filesystem::path resume_path;
        bool dry_run = false;
        for (int index = 1; index < argc; ++index) {
            const std::string argument = argv[index];
            if (argument == "--help" || argument == "-h") {
                print_help(argv[0]);
                return 0;
            }
            if (argument == "--config" && index + 1 < argc) {
                config_path = argv[++index];
            } else if (argument == "--set" && index + 1 < argc) {
                overrides.emplace_back(argv[++index]);
            } else if (argument == "--resume" && index + 1 < argc) {
                resume_path = argv[++index];
            } else if (argument == "--dry-run") {
                dry_run = true;
            } else {
                throw std::invalid_argument("unknown or incomplete argument: " + argument);
            }
        }

        atcg3d::Model3DConfig config = atcg3d::Model3DConfig::load(config_path);
        for (const std::string& override_value : overrides) {
            config.apply_override(override_value);
        }
        if (dry_run) {
            std::cout << config.to_json();
            return 0;
        }

        atcg3d::Simulation3D simulation(config);
        if (!resume_path.empty()) {
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
            const atcg3d::CheckpointData3D checkpoint =
                atcg3d::read_hdf5_checkpoint(resume_path, config);
            simulation.restore(checkpoint.cells, checkpoint.next_uid, checkpoint.clock,
                               checkpoint.stats, checkpoint.lineage);
#else
            throw std::runtime_error(
                "--resume requires a build configured with ATCG3D_ENABLE_HDF5_CHECKPOINT=ON");
#endif
        }
        atcg3d::OutputManager3D output(config);
        simulation.run([&output](const atcg3d::Simulation3D& current) { output.observe(current); });
        output.finalize(simulation);
        const auto& stats = simulation.stats();
        std::cout << "ATCG3D completed\n"
                  << "time_hours=" << simulation.clock().time_hours << '\n'
                  << "completed_events=" << simulation.clock().completed_events << '\n'
                  << "alive_cells=" << simulation.cells().alive_count() << '\n'
                  << "chunks=" << simulation.grid().chunk_count() << '\n'
                  << "cell_store_bytes=" << simulation.cells().allocated_bytes() << '\n'
                  << "grid_bytes=" << simulation.grid().allocated_bytes() << '\n'
                  << "migration_attempts=" << stats.migration_attempts << '\n'
                  << "migration_commits=" << stats.migration_commits << '\n'
                  << "divisions=" << stats.divisions << '\n'
                  << "deaths=" << stats.deaths << '\n'
                  << "checksum=" << simulation.state_checksum() << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "atcg3d: " << error.what() << '\n';
        return 1;
    }
}
