#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>

#include "config/model_config.hpp"
#include "core/stateless_rng.hpp"
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
        std::uint64_t cells = 1000;
        std::uint64_t max_events = 100000;
        int threads = 1;
        for (int index = 1; index < argc; ++index) {
            const std::string argument = argv[index];
            if (argument == "--cells" && index + 1 < argc) cells = std::stoull(argv[++index]);
            else if (argument == "--events" && index + 1 < argc) max_events = std::stoull(argv[++index]);
            else if (argument == "--threads" && index + 1 < argc) threads = std::stoi(argv[++index]);
            else throw std::invalid_argument("unknown or incomplete argument: " + argument);
        }
        if (cells < 2 || max_events == 0 || threads <= 0) {
            throw std::invalid_argument("cells/events/threads are outside valid ranges");
        }
        atcg3d::Model3DConfig config;
        config.output_enabled = false;
        config.initial_r_cells = 0;
        config.initial_K_cells = 0;
        config.migration_activation_enabled = false;
        config.angiogenesis.enabled = false;
        config.end_time_hours = 1000000.0;
        config.max_events = max_events;
        config.threads = threads;

        // Build real typed cells and a real sparse occupancy/density index.
        // A three-voxel lattice spacing leaves migration room, and deterministic
        // sub-hour event times avoid a synthetic all-at-once batch. Division
        // work is intentionally far in the future: this scale benchmark
        // advances migration biology rather than merely changing metadata.
        const auto build_start = std::chrono::steady_clock::now();
        std::vector<atcg3d::CellInit> records;
        records.reserve(static_cast<std::size_t>(cells));
        const std::uint64_t side = static_cast<std::uint64_t>(
            std::ceil(std::cbrt(static_cast<long double>(cells))));
        for (std::uint64_t index = 0; index < cells; ++index) {
            atcg3d::CellInit cell;
            cell.uid = index + 1;
            cell.clone_id = static_cast<std::uint32_t>((index % 1000000ULL) + 1ULL);
            cell.type = index % 2 == 0 ? atcg3d::CellType::r : atcg3d::CellType::K;
            cell.stage = atcg3d::CellStage::small;
            cell.anchor = {
                static_cast<std::int32_t>(3ULL * (index % side)),
                static_cast<std::int32_t>(3ULL * ((index / side) % side)),
                static_cast<std::int32_t>(3ULL * (index / (side * side)))};
            cell.flags = static_cast<std::uint8_t>(atcg3d::kDirtyDensity);
            cell.inherent_growth_rate = 1.0F;
            cell.density_growth_rate = 1.0F;
            cell.migration_rate = 1.0F;
            cell.normal_migration_rate = 1.0F;
            cell.division_work_remaining = 1000000.0F;
            cell.next_migration_time = 0.01 + 0.99 * atcg3d::rng_unit(
                config.seed, cell.uid, 0x45564e5442454e43ULL, 0, 0);
            cell.next_division_time = 1000000.0;
            cell.migration_schedule_generation = 1;
            cell.division_schedule_generation = 1;
            cell.death_schedule_generation = 1;
            records.push_back(cell);
        }
        atcg3d::Simulation3D simulation(config);
        simulation.restore(records, cells + 1, {}, {}, {});
        records.clear();
        records.shrink_to_fit();
        const double build_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - build_start).count();

        const auto run_start = std::chrono::steady_clock::now();
        simulation.run();
        const double run_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - run_start).count();
        if (simulation.clock().completed_events == 0) {
            throw std::runtime_error("event benchmark advanced no biological events");
        }
        const double rate = run_seconds > 0.0
            ? simulation.clock().completed_events / run_seconds : 0.0;
        std::cout << std::setprecision(17)
                  << "{\n"
                  << "  \"schema\": \"atcg3d.events.v2\",\n"
                  << "  \"initial_cells\": " << cells << ",\n"
                  << "  \"final_cells\": " << simulation.cells().alive_count() << ",\n"
                  << "  \"completed_events\": " << simulation.clock().completed_events << ",\n"
                  << "  \"simulated_hours\": " << simulation.clock().time_hours << ",\n"
                  << "  \"threads\": " << threads << ",\n"
                  << "  \"build_seconds\": " << build_seconds << ",\n"
                  << "  \"run_seconds\": " << run_seconds << ",\n"
                  << "  \"events_per_second\": " << rate << ",\n"
                  << "  \"migration_attempts\": "
                  << simulation.stats().migration_attempts << ",\n"
                  << "  \"migration_commits\": "
                  << simulation.stats().migration_commits << ",\n"
                  << "  \"pending_events\": "
                  << simulation.pending_event_count() << ",\n"
                  << "  \"event_queue_rebuilds\": "
                  << simulation.event_queue_rebuild_count() << ",\n"
                  << "  \"peak_rss_bytes\": " << peak_rss_bytes() << ",\n"
                  << "  \"checksum\": " << simulation.state_checksum() << "\n"
                  << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "atcg3d_event_benchmark: " << error.what() << '\n';
        return 1;
    }
}
