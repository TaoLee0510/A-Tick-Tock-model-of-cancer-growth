#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>

#include "config/model_config.hpp"
#include "engine/simulation.hpp"
#include "io/preview_sampler.hpp"
#include "io/snapshot.hpp"
#include "io/vtkhdf_writer.hpp"
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
#include "io/checkpoint_hdf5.hpp"
#endif

namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

std::uint64_t peak_rss_bytes() {
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) return 0;
#ifdef __APPLE__
    return static_cast<std::uint64_t>(usage.ru_maxrss);
#else
    return static_cast<std::uint64_t>(usage.ru_maxrss) * 1024ULL;
#endif
}

struct Options {
    std::uint64_t cells{100000};
    std::uint64_t preview_max{1000000};
    std::filesystem::path directory{"atcg3d_scale_result"};
    std::filesystem::path result{"atcg3d_scale_result/result.json"};
    bool write_vtk{true};
    bool write_checkpoint{true};
    bool require_vtk{};
    bool require_checkpoint{};
    int compression_level{1};
};

std::uint64_t parse_u64(const std::string& value, const char* name) {
    std::size_t consumed = 0;
    const std::uint64_t result = std::stoull(value, &consumed);
    if (consumed != value.size() || result == 0) {
        throw std::invalid_argument(std::string(name) + " must be a positive integer");
    }
    return result;
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--cells" && index + 1 < argc) {
            options.cells = parse_u64(argv[++index], "--cells");
        } else if (argument == "--preview-max" && index + 1 < argc) {
            options.preview_max = parse_u64(argv[++index], "--preview-max");
        } else if (argument == "--directory" && index + 1 < argc) {
            options.directory = argv[++index];
        } else if (argument == "--result" && index + 1 < argc) {
            options.result = argv[++index];
        } else if (argument == "--skip-vtk") {
            options.write_vtk = false;
        } else if (argument == "--skip-checkpoint") {
            options.write_checkpoint = false;
        } else if (argument == "--require-vtk") {
            options.require_vtk = true;
        } else if (argument == "--require-checkpoint") {
            options.require_checkpoint = true;
        } else if (argument == "--compression-level" && index + 1 < argc) {
            const std::string value = argv[++index];
            std::size_t consumed = 0;
            options.compression_level = std::stoi(value, &consumed);
            if (consumed != value.size() || options.compression_level < 0 ||
                options.compression_level > 9) {
                throw std::invalid_argument(
                    "--compression-level must be an integer in [0,9]");
            }
        } else {
            throw std::invalid_argument("unknown or incomplete argument: " + argument);
        }
    }
    return options;
}

template <class Writer>
void write_atomic(const std::filesystem::path& path, Writer writer) {
    std::filesystem::create_directories(path.parent_path());
    if (std::filesystem::exists(path)) {
        throw std::runtime_error("refusing to overwrite benchmark result: " + path.string());
    }
    const auto temporary = path.string() + ".tmp";
    std::filesystem::remove(temporary);
    {
        std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
        if (!stream) throw std::runtime_error("unable to create " + temporary);
        writer(stream);
        stream.flush();
        if (!stream) throw std::runtime_error("unable to finish " + temporary);
    }
    std::filesystem::rename(temporary, path);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        using namespace atcg3d;
        const Options options = parse_options(argc, argv);
        if (std::filesystem::exists(options.directory)) {
            throw std::runtime_error("refusing to overwrite benchmark directory: " +
                                     options.directory.string());
        }
        std::filesystem::create_directories(options.directory);

        Model3DConfig config;
        config.output_enabled = false;
        config.initial_r_cells = 0;
        config.initial_K_cells = 0;
        config.max_events = 1;
        config.end_time_hours = 0.0;
        config.chunk_edge = 32;
        config.density_block_edge = 4;
        config.vtkhdf_compression_level = options.compression_level;
        config.hdf5_compression_level = options.compression_level;

        const auto build_start = Clock::now();
        std::vector<CellInit> records;
        records.reserve(static_cast<std::size_t>(options.cells));
        const std::uint64_t side = static_cast<std::uint64_t>(
            std::ceil(std::cbrt(static_cast<long double>(options.cells))));
        for (std::uint64_t index = 0; index < options.cells; ++index) {
            CellInit cell;
            cell.uid = index + 1;
            cell.clone_id = static_cast<std::uint32_t>((index % 1000000ULL) + 1ULL);
            cell.type = index % 2 == 0 ? CellType::r : CellType::K;
            cell.stage = CellStage::small;
            cell.anchor = {static_cast<std::int32_t>(index % side),
                           static_cast<std::int32_t>((index / side) % side),
                           static_cast<std::int32_t>(index / (side * side))};
            records.push_back(cell);
        }
        Simulation3D simulation(config);
        simulation.restore(records, options.cells + 1, {}, {}, {});
        records.clear();
        records.shrink_to_fit();
        const double build_seconds = seconds_since(build_start);

        const auto preview_start = Clock::now();
        const std::vector<Slot> preview = stable_preview_sample(
            simulation.cells(), static_cast<std::size_t>(
                std::min(options.preview_max, options.cells)), 0x3dULL);
        const double preview_sample_seconds = seconds_since(preview_start);

        bool vtk_written = false;
        double vtk_write_seconds = 0.0;
        double vtk_read_seconds = 0.0;
        std::uint64_t full_bytes = 0;
        std::uint64_t preview_bytes = 0;
        if (options.write_vtk && vtkhdf_output_available()) {
            const std::vector<Slot> all = simulation.cells().alive_slots();
            const SimulationSnapshotView3D full_view{
                simulation.cells(), all, simulation.lesion_index(), {}};
            const auto full_start = Clock::now();
            const auto full_path = options.directory / "full.vtkhdf";
            write_vtkhdf_points_atomic(full_path, full_view,
                                        options.compression_level);
            vtk_write_seconds = seconds_since(full_start);
            full_bytes = std::filesystem::file_size(full_path);
            const SimulationSnapshotView3D preview_view{
                simulation.cells(), preview, simulation.lesion_index(), {}};
            const auto preview_path = options.directory / "preview.vtkhdf";
            write_vtkhdf_points_atomic(preview_path, preview_view,
                                        options.compression_level);
            preview_bytes = std::filesystem::file_size(preview_path);
            vtk_written = true;
            const auto read_start = Clock::now();
            if (read_vtkhdf_point_count(full_path) != options.cells) {
                throw std::runtime_error("VTK-HDF full reader round-trip count mismatch");
            }
            vtk_read_seconds = seconds_since(read_start);
        } else if (options.require_vtk) {
            throw std::runtime_error("--require-vtk specified but VTK-HDF is unavailable");
        }

        bool checkpoint_written = false;
        double checkpoint_write_seconds = 0.0;
        double checkpoint_read_restore_seconds = 0.0;
        std::uint64_t checkpoint_bytes = 0;
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
        if (options.write_checkpoint) {
            const auto checkpoint_path = options.directory / "checkpoint.h5";
            const auto write_start = Clock::now();
            write_hdf5_checkpoint(checkpoint_path, simulation);
            checkpoint_write_seconds = seconds_since(write_start);
            checkpoint_bytes = std::filesystem::file_size(checkpoint_path);
            const auto restore_start = Clock::now();
            const CheckpointData3D data = read_hdf5_checkpoint(checkpoint_path, config);
            Simulation3D restored(config);
            restored.restore(data.cells, data.next_uid, data.clock, data.stats,
                             data.lineage, data.vasculature, data.cell_slot_count,
                             data.cell_slots, data.cell_free_slots);
            if (restored.cells().alive_count() != options.cells ||
                restored.state_checksum() != simulation.state_checksum()) {
                throw std::runtime_error("checkpoint restore checksum mismatch");
            }
            checkpoint_read_restore_seconds = seconds_since(restore_start);
            checkpoint_written = true;
        }
#else
        if (options.require_checkpoint) {
            throw std::runtime_error(
                "--require-checkpoint specified but HDF5 checkpoint support is unavailable");
        }
#endif

        const std::uint64_t cell_store_bytes = simulation.cells().allocated_bytes();
        const std::uint64_t grid_bytes = simulation.grid().allocated_bytes();
        const std::uint64_t density_bytes = simulation.density().allocated_bytes();
        const double bytes_per_cell = static_cast<double>(cell_store_bytes + grid_bytes + density_bytes) /
                                      static_cast<double>(options.cells);
        const std::uint64_t peak = peak_rss_bytes();

        std::ostringstream json;
        json << std::setprecision(17)
             << "{\n"
             << "  \"schema\": \"atcg3d.synthetic_scale.v1\",\n"
             << "  \"cells\": " << options.cells << ",\n"
             << "  \"compression_level\": "
             << options.compression_level << ",\n"
             << "  \"cell_store_logical_bytes_per_slot\": "
             << CellStore3D::logical_bytes_per_slot() << ",\n"
             << "  \"cell_store_bytes\": " << cell_store_bytes << ",\n"
             << "  \"grid_bytes\": " << grid_bytes << ",\n"
             << "  \"density_index_bytes\": " << density_bytes << ",\n"
             << "  \"core_bytes_per_cell\": " << bytes_per_cell << ",\n"
             << "  \"peak_rss_bytes\": " << peak << ",\n"
             << "  \"build_seconds\": " << build_seconds << ",\n"
             << "  \"preview_sample_seconds\": " << preview_sample_seconds << ",\n"
             << "  \"preview_cells\": " << preview.size() << ",\n"
             << "  \"vtk_written\": " << (vtk_written ? "true" : "false") << ",\n"
             << "  \"vtk_full_write_seconds\": " << vtk_write_seconds << ",\n"
             << "  \"vtk_full_read_seconds\": " << vtk_read_seconds << ",\n"
             << "  \"vtk_full_bytes\": " << full_bytes << ",\n"
             << "  \"vtk_preview_bytes\": " << preview_bytes << ",\n"
             << "  \"checkpoint_written\": " << (checkpoint_written ? "true" : "false") << ",\n"
             << "  \"checkpoint_write_seconds\": " << checkpoint_write_seconds << ",\n"
             << "  \"checkpoint_read_restore_seconds\": " << checkpoint_read_restore_seconds << ",\n"
             << "  \"checkpoint_bytes\": " << checkpoint_bytes << "\n"
             << "}\n";
        write_atomic(options.result, [&](std::ostream& out) { out << json.str(); });
        std::cout << json.str();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "atcg3d_scale_benchmark: " << error.what() << '\n';
        return 1;
    }
}
