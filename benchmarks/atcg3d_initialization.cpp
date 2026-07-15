#include <chrono>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>

#include "config/model_config.hpp"
#include "engine/simulation.hpp"

namespace {

std::uint64_t peak_rss_bytes() {
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) return 0;
#ifdef __APPLE__
    return static_cast<std::uint64_t>(usage.ru_maxrss);
#else
    return static_cast<std::uint64_t>(usage.ru_maxrss) * 1024ULL;
#endif
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3 || std::string(argv[1]) != "--config") {
            throw std::invalid_argument("usage: atcg3d_initialization_benchmark --config PATH");
        }
        atcg3d::Model3DConfig config = atcg3d::Model3DConfig::load(
            std::filesystem::path(argv[2]));
        config.output_enabled = false;
        config.run_mode = "new";
        config.resume_checkpoint.clear();
        const auto start = std::chrono::steady_clock::now();
        atcg3d::Simulation3D simulation(config);
        simulation.initialize();
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();

        std::uint64_t r_cells = 0;
        std::uint64_t K_cells = 0;
        for (const atcg3d::Slot slot : simulation.cells().alive_slots()) {
            if (simulation.cells().type(slot) == atcg3d::CellType::r) ++r_cells;
            else ++K_cells;
        }
        std::cout << std::setprecision(17)
                  << "{\n"
                  << "  \"schema\": \"atcg3d.initialization.v1\",\n"
                  << "  \"profile\": \"" << config.profile << "\",\n"
                  << "  \"initialization_mode\": \"" << config.initialization_mode << "\",\n"
                  << "  \"cells\": " << simulation.cells().alive_count() << ",\n"
                  << "  \"r_cells\": " << r_cells << ",\n"
                  << "  \"K_cells\": " << K_cells << ",\n"
                  << "  \"stage0_large\": "
                  << simulation.cells().stage_count(atcg3d::CellStage::large) << ",\n"
                  << "  \"stage1_small\": "
                  << simulation.cells().stage_count(atcg3d::CellStage::small) << ",\n"
                  << "  \"stage2_ultrasmall\": "
                  << simulation.cells().stage_count(atcg3d::CellStage::ultrasmall) << ",\n"
                  << "  \"biological_volume_voxels3\": "
                  << simulation.biological_tumor_volume() << ",\n"
                  << "  \"cell_chunks\": " << simulation.grid().chunk_count() << ",\n"
                  << "  \"surface_faces\": " << simulation.tumor_surface().size() << ",\n"
                  << "  \"cell_store_bytes\": " << simulation.cells().allocated_bytes() << ",\n"
                  << "  \"grid_bytes\": " << simulation.grid().allocated_bytes() << ",\n"
                  << "  \"density_bytes\": " << simulation.density().allocated_bytes() << ",\n"
                  << "  \"peak_rss_bytes\": " << peak_rss_bytes() << ",\n"
                  << "  \"wall_seconds\": " << seconds << ",\n"
                  << "  \"checksum\": " << simulation.state_checksum() << "\n"
                  << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "atcg3d_initialization_benchmark: " << error.what() << '\n';
        return 1;
    }
}
