#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>

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

std::uint64_t parse_unsigned(const std::string& text, const char* name) {
    if (text.empty() || text.front() == '-') {
        throw std::invalid_argument(std::string(name) + " must be unsigned");
    }
    std::size_t consumed = 0;
    const std::uint64_t value = std::stoull(text, &consumed);
    if (consumed != text.size()) {
        throw std::invalid_argument(std::string(name) + " is not an integer");
    }
    return value;
}

int parse_threads(const std::string& text) {
    const std::uint64_t value = parse_unsigned(text, "threads");
    if (value == 0 || value > static_cast<std::uint64_t>(
            std::numeric_limits<int>::max())) {
        throw std::invalid_argument("threads are outside the supported range");
    }
    return static_cast<int>(value);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::uint64_t cell_count = 1000;
        std::uint64_t max_events = 256;
        int threads = 1;
        std::string profile = "synthetic";
        for (int index = 1; index < argc; ++index) {
            const std::string argument = argv[index];
            if (argument == "--cells" && index + 1 < argc) {
                cell_count = parse_unsigned(argv[++index], "cells");
            } else if (argument == "--events" && index + 1 < argc) {
                max_events = parse_unsigned(argv[++index], "events");
            } else if (argument == "--threads" && index + 1 < argc) {
                threads = parse_threads(argv[++index]);
            } else if (argument == "--profile" && index + 1 < argc) {
                profile = argv[++index];
            } else {
                throw std::invalid_argument(
                    "unknown or incomplete argument: " + argument);
            }
        }
        if (cell_count == 0 ||
            cell_count >= static_cast<std::uint64_t>(atcg3d::kEmptySlot) ||
            max_events == 0 || threads <= 0 ||
            (profile != "synthetic" && profile != "production")) {
            throw std::invalid_argument(
                "cells/events/threads/profile are outside valid ranges");
        }
        const bool production_profile = profile == "production";

        atcg3d::Model3DConfig config;
        config.output_enabled = false;
        config.initial_r_cells = 0;
        config.initial_K_cells = 0;
        config.migration_activation_enabled = false;
        config.end_time_hours = production_profile ? 720.0 : 24.0;
        config.max_events = max_events;
            config.threads = threads;
            config.parallel_min_events_per_thread = 1;
            config.parallel_thread_thresholds = {{0, 1.0}};

        config.angiogenesis.enabled = true;
        if (production_profile) {
            // Keep the production legacy-mapped vascular load: true Poisson
            // intensity, root/tip limits, diameter, speeds, stopping lengths,
            // collision policy, relief radius, density mapping, and death
            // parameters. Only the synthetic cells' initial schedules below
            // are neutralized so this remains an interpretable vascular scale
            // benchmark rather than a 30-day whole-model calibration. Vessel
            // displacement still removes real cells from the real sparse
            // occupancy/density indexes.
            config.angiogenesis.seed_rate_sites_per_30_days = 10.0;
            config.angiogenesis.seed_rate_sites_per_hour = 10.0 / 720.0;
        } else {
            config.angiogenesis.lesion_block_edge = 8;
            config.angiogenesis.lesion_core_activation_occupied_fraction =
                1.0 / 512.0;
            config.angiogenesis.lesion_core_deactivation_occupied_fraction = 0.0;
            config.angiogenesis.lesion_minimum_cells_per_core_block = 1;
            config.angiogenesis.trigger_minimum_core_blocks = 1;
            config.angiogenesis.trigger_activation_volume_voxels3 = 1.0;
            config.angiogenesis.trigger_deactivation_volume_voxels3 = 0.0;
            config.angiogenesis.trigger_delay_hours = 0.0;
            config.angiogenesis.seed_rate_sites_per_30_days = 720000.0;
            config.angiogenesis.seed_rate_sites_per_hour = 1000.0;
            config.angiogenesis.roots_per_event = 1;
            config.angiogenesis.surface_min_separation_voxels = 0;
            config.angiogenesis.surface_max_sampling_attempts = 4096;
            config.angiogenesis.max_total_roots = 1;
            config.angiogenesis.max_active_tips = 2;
            config.angiogenesis.diameter_voxels = 1.0;
            config.angiogenesis.inward_speed_voxels_per_hour = 4.0;
            config.angiogenesis.outward_speed_voxels_per_hour = 4.0;
            config.angiogenesis.inward_max_length_voxels = 12;
            config.angiogenesis.outward_max_length_voxels = 12;
            config.angiogenesis.inward_target_tolerance_voxels = 0.0;
            config.angiogenesis.outward_external_connection_distance_voxels = 4.0;
            config.angiogenesis.influence_activation = "immediate";
        }
        config.validate();

        const auto build_start = std::chrono::steady_clock::now();
        std::vector<atcg3d::CellInit> records;
        records.reserve(static_cast<std::size_t>(cell_count));
        const std::uint64_t side = static_cast<std::uint64_t>(
            std::ceil(std::cbrt(static_cast<long double>(cell_count))));
        for (std::uint64_t index = 0; index < cell_count; ++index) {
            atcg3d::CellInit cell;
            cell.uid = index + 1;
            cell.clone_id = static_cast<std::uint32_t>(
                (index % 1000000ULL) + 1ULL);
            cell.type = index % 2 == 0
                ? atcg3d::CellType::r : atcg3d::CellType::K;
            cell.stage = atcg3d::CellStage::small;
            cell.anchor = {
                static_cast<std::int32_t>(index % side),
                static_cast<std::int32_t>((index / side) % side),
                static_cast<std::int32_t>(index / (side * side))};
            cell.flags = static_cast<std::uint8_t>(atcg3d::kDirtyDensity);
            cell.inherent_growth_rate = 1.0F;
            cell.density_growth_rate = production_profile ? 0.0F : 1.0F;
            cell.migration_rate = 0.0F;
            cell.normal_migration_rate = 0.0F;
            cell.division_work_remaining = 1000000000.0F;
            cell.next_migration_time = 0.0;
            cell.next_division_time = production_profile ? 0.0 : 1000000000.0;
            cell.death_deadline = production_profile ? 1000000000.0 : 0.0;
            cell.migration_schedule_generation = 1;
            cell.division_schedule_generation = 1;
            cell.death_schedule_generation = 1;
            records.push_back(cell);
        }

        atcg3d::Simulation3D simulation(config);
        simulation.restore(records, cell_count + 1, {}, {}, {});
        const auto process_state = simulation.angiogenesis_state();
        if (!process_state.eligible ||
            !(process_state.next_seed_time_hours > 0.0) ||
            process_state.next_seed_time_hours >= config.end_time_hours) {
            throw std::runtime_error(
                "per-lesion Poisson seed falls outside the benchmark horizon");
        }
        const double scheduled_seed_time_hours =
            process_state.next_seed_time_hours;
        const std::size_t initial_surface_faces = simulation.tumor_surface().size();
        const std::size_t initial_chunks = simulation.grid().chunk_count();
        std::vector<atcg3d::CellInit>().swap(records);
        const double build_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - build_start).count();

        const auto run_start = std::chrono::steady_clock::now();
        simulation.run();
        const double run_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - run_start).count();
        const atcg3d::SimulationStats3D& stats = simulation.stats();
        std::size_t root_nodes = 0;
        std::size_t inward_nodes = 0;
        std::size_t outward_nodes = 0;
        for (const atcg3d::VesselNodeSlot slot :
             simulation.vessel_nodes().alive_slots()) {
            switch (simulation.vessel_nodes().role(slot)) {
                case atcg3d::VesselBranchRole::root:
                    ++root_nodes;
                    break;
                case atcg3d::VesselBranchRole::inward:
                    ++inward_nodes;
                    break;
                case atcg3d::VesselBranchRole::outward:
                    ++outward_nodes;
                    break;
            }
        }
        const bool valid_seed_progress = production_profile
            ? stats.angiogenesis_seed_attempts > 0 &&
                  stats.angiogenesis_roots > 0 &&
                  stats.angiogenesis_seed_attempts ==
                      stats.angiogenesis_roots +
                          stats.angiogenesis_seed_rejections &&
                  std::abs(simulation.clock().time_hours -
                           config.end_time_hours) <= 1.0e-9
            : stats.angiogenesis_seed_attempts == 1 &&
                  stats.angiogenesis_roots == 1 &&
                  stats.angiogenesis_seed_rejections == 0;
        const std::size_t alive_tips = simulation.vessel_tips().alive_count();
        const std::size_t active_tips = simulation.active_vessel_tip_count();
        const std::size_t occupied_vessel_voxels =
            simulation.vessel_grid().occupied_voxel_count();
        const std::size_t influenced_voxels =
            simulation.vascular_influence().influenced_voxel_count();
        if (!valid_seed_progress ||
            stats.vessel_growth_attempts == 0 ||
            stats.vessel_growth_commits == 0 ||
            stats.vascular_displacements == 0 ||
            root_nodes != stats.angiogenesis_roots ||
            inward_nodes == 0 || outward_nodes == 0 ||
            (!production_profile && alive_tips != root_nodes * 2) ||
            occupied_vessel_voxels == 0 || influenced_voxels == 0) {
            std::ostringstream detail;
            detail << "benchmark did not advance the required seed process, "
                      "bidirectional growth, displacement, and vascular influence"
                   << " [completed_events=" << simulation.clock().completed_events
                   << ", time_hours=" << simulation.clock().time_hours
                   << ", seed_attempts=" << stats.angiogenesis_seed_attempts
                   << ", roots=" << stats.angiogenesis_roots
                   << ", seed_rejections=" << stats.angiogenesis_seed_rejections
                   << ", growth_attempts=" << stats.vessel_growth_attempts
                   << ", growth_commits=" << stats.vessel_growth_commits
                   << ", displacements=" << stats.vascular_displacements
                   << ", root_nodes=" << root_nodes
                   << ", inward_nodes=" << inward_nodes
                   << ", outward_nodes=" << outward_nodes
                   << ", alive_tips=" << alive_tips
                   << ", active_tips=" << active_tips
                   << ", occupied_voxels=" << occupied_vessel_voxels
                   << ", influenced_voxels=" << influenced_voxels;
            for (const atcg3d::VesselTipSlot slot :
                 simulation.vessel_tips().alive_slots()) {
                const atcg3d::Vec3i position =
                    simulation.vessel_tips().position(slot);
                detail << ", tip{uid=" << simulation.vessel_tips().uid(slot)
                       << ",role=" << static_cast<int>(
                              simulation.vessel_tips().role(slot))
                       << ",status=" << static_cast<int>(
                              simulation.vessel_tips().status(slot))
                       << ",position=" << position.x << '/' << position.y
                       << '/' << position.z
                       << ",grown="
                       << simulation.vessel_tips().grown_length_voxels(slot)
                       << '}';
            }
            detail << ']';
            throw std::runtime_error(detail.str());
        }

        const double event_rate = run_seconds > 0.0
            ? static_cast<double>(simulation.clock().completed_events) /
                  run_seconds
            : 0.0;
        const std::uint64_t cell_store_bytes = simulation.cells().allocated_bytes();
        const std::uint64_t cell_grid_bytes = simulation.grid().allocated_bytes();
        const std::uint64_t density_index_bytes = simulation.density().allocated_bytes();
        const std::uint64_t vessel_grid_bytes = simulation.vessel_grid().allocated_bytes();
        const std::uint64_t vessel_node_store_bytes =
            simulation.vessel_nodes().allocated_bytes();
        const std::uint64_t vessel_tip_store_bytes =
            simulation.vessel_tips().allocated_bytes();
        const std::uint64_t vascular_influence_bytes =
            simulation.vascular_influence().allocated_bytes();
        const std::uint64_t lesion_index_bytes =
            simulation.lesion_index().allocated_bytes();
        // This subtotal deliberately covers only components with explicit
        // allocation accounting. Surface hash tables, scheduler storage, and
        // allocator overhead are represented by peak_rss_bytes instead.
        const std::uint64_t tracked_component_bytes =
            cell_store_bytes + cell_grid_bytes + density_index_bytes +
            vessel_grid_bytes + vessel_node_store_bytes +
            vessel_tip_store_bytes + vascular_influence_bytes +
            lesion_index_bytes;
        std::cout << std::setprecision(17)
                  << "{\n"
                  << "  \"schema\": \"atcg3d.angiogenesis_scale.v1\",\n"
                  << "  \"profile\": \"" << profile << "\",\n"
                  << "  \"initial_cells\": " << cell_count << ",\n"
                  << "  \"final_cells\": "
                  << simulation.cells().alive_count() << ",\n"
                  << "  \"threads\": " << threads << ",\n"
                  << "  \"initial_surface_faces\": "
                  << initial_surface_faces << ",\n"
                  << "  \"initial_chunks\": " << initial_chunks << ",\n"
                  << "  \"lesions\": "
                  << simulation.lesion_index().lesions().size() << ",\n"
                  << "  \"completed_events\": "
                  << simulation.clock().completed_events << ",\n"
                  << "  \"simulated_hours\": "
                  << simulation.clock().time_hours << ",\n"
                  << "  \"background_cell_schedules_neutralized\": "
                  << (production_profile ? "true" : "false") << ",\n"
                  << "  \"benchmark_initial_death_deadline_hours\": "
                  << (production_profile ? 1000000000.0 : 0.0) << ",\n"
                  << "  \"seed_rate_sites_per_30_days\": "
                  << config.angiogenesis.seed_rate_sites_per_30_days << ",\n"
                  << "  \"trigger_activation_volume_voxels3\": "
                  << config.angiogenesis.trigger_activation_volume_voxels3 << ",\n"
                  << "  \"trigger_deactivation_volume_voxels3\": "
                  << config.angiogenesis.trigger_deactivation_volume_voxels3 << ",\n"
                  << "  \"max_total_roots\": "
                  << config.angiogenesis.max_total_roots << ",\n"
                  << "  \"max_active_tips\": "
                  << config.angiogenesis.max_active_tips << ",\n"
                  << "  \"diameter_voxels\": "
                  << config.angiogenesis.diameter_voxels << ",\n"
                  << "  \"inward_speed_voxels_per_hour\": "
                  << config.angiogenesis.inward_speed_voxels_per_hour << ",\n"
                  << "  \"outward_speed_voxels_per_hour\": "
                  << config.angiogenesis.outward_speed_voxels_per_hour << ",\n"
                  << "  \"inward_max_length_voxels\": "
                  << config.angiogenesis.inward_max_length_voxels << ",\n"
                  << "  \"outward_max_length_voxels\": "
                  << config.angiogenesis.outward_max_length_voxels << ",\n"
                  << "  \"surface_min_separation_voxels\": "
                  << config.angiogenesis.surface_min_separation_voxels << ",\n"
                  << "  \"outward_external_connection_distance_voxels\": "
                  << config.angiogenesis.outward_external_connection_distance_voxels << ",\n"
                  << "  \"influence_cutoff_radius_voxels\": "
                  << config.angiogenesis.influence_cutoff_radius_voxels << ",\n"
                  << "  \"influence_max_relief_fraction\": "
                  << config.angiogenesis.influence_max_relief_fraction << ",\n"
                  << "  \"scheduled_seed_time_hours\": "
                  << scheduled_seed_time_hours << ",\n"
                  << "  \"seed_attempts\": "
                  << stats.angiogenesis_seed_attempts << ",\n"
                  << "  \"committed_roots\": "
                  << stats.angiogenesis_roots << ",\n"
                  << "  \"seed_rejections\": "
                  << stats.angiogenesis_seed_rejections << ",\n"
                  << "  \"vessel_growth_attempts\": "
                  << stats.vessel_growth_attempts << ",\n"
                  << "  \"vessel_growth_commits\": "
                  << stats.vessel_growth_commits << ",\n"
                  << "  \"cell_divisions\": " << stats.divisions << ",\n"
                  << "  \"cell_deaths\": " << stats.deaths << ",\n"
                  << "  \"vascular_displacements\": "
                  << stats.vascular_displacements << ",\n"
                  << "  \"vessel_nodes\": "
                  << simulation.vessel_nodes().alive_count() << ",\n"
                  << "  \"root_nodes\": " << root_nodes << ",\n"
                  << "  \"inward_nodes\": " << inward_nodes << ",\n"
                  << "  \"outward_nodes\": " << outward_nodes << ",\n"
                  << "  \"vessel_tips\": "
                  << alive_tips << ",\n"
                  << "  \"active_vessel_tips\": "
                  << active_tips << ",\n"
                  << "  \"vessel_occupied_voxels\": "
                  << occupied_vessel_voxels << ",\n"
                  << "  \"vascular_influenced_voxels\": "
                  << influenced_voxels << ",\n"
                  << "  \"build_seconds\": " << build_seconds << ",\n"
                  << "  \"run_seconds\": " << run_seconds << ",\n"
                  << "  \"events_per_second\": " << event_rate << ",\n"
                  << "  \"cell_store_bytes\": "
                  << cell_store_bytes << ",\n"
                  << "  \"cell_grid_bytes\": "
                  << cell_grid_bytes << ",\n"
                  << "  \"density_index_bytes\": "
                  << density_index_bytes << ",\n"
                  << "  \"vessel_grid_bytes\": "
                  << vessel_grid_bytes << ",\n"
                  << "  \"vessel_node_store_bytes\": "
                  << vessel_node_store_bytes << ",\n"
                  << "  \"vessel_tip_store_bytes\": "
                  << vessel_tip_store_bytes << ",\n"
                  << "  \"vascular_influence_bytes\": "
                  << vascular_influence_bytes << ",\n"
                  << "  \"lesion_index_bytes\": "
                  << lesion_index_bytes << ",\n"
                  << "  \"tracked_component_bytes\": "
                  << tracked_component_bytes << ",\n"
                  << "  \"tracked_bytes_per_initial_cell\": "
                  << static_cast<double>(tracked_component_bytes) /
                         static_cast<double>(cell_count) << ",\n"
                  << "  \"peak_rss_bytes\": " << peak_rss_bytes() << ",\n"
                  << "  \"checksum\": " << simulation.state_checksum() << "\n"
                  << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "atcg3d_angiogenesis_scale_benchmark: "
                  << error.what() << '\n';
        return 1;
    }
}
