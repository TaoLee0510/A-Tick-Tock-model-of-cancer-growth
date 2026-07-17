#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>

#include "config/model_config.hpp"
#include "engine/simulation.hpp"

namespace {

atcg3d::Model3DConfig test_config(int threads) {
    atcg3d::Model3DConfig config;
    config.output_enabled = false;
    config.initial_r_cells = 20;
    config.initial_K_cells = 20;
    config.initial_radius = 8;
    config.initial_shell_inner_radius = 6;
    config.initial_shell_thickness = 2;
    config.initial_inner_small_radius = 6;
    config.initial_r_migration_rate = 0.0;
    config.initial_K_migration_rate = 0.0;
    config.end_time_hours = 2.0;
    config.max_events = 100000;
    config.threads = threads;
    config.parallel_min_events_per_thread = 1;
    config.parallel_thread_thresholds = {{0, 1.0}};

    auto& vessels = config.angiogenesis;
    vessels.enabled = true;
    vessels.trigger_activation_volume_voxels3 = 1.0;
    vessels.trigger_deactivation_volume_voxels3 = 0.0;
    vessels.seed_rate_sites_per_30_days = 720000.0;
    vessels.seed_rate_sites_per_hour = vessels.seed_rate_sites_per_30_days / 720.0;
    vessels.max_total_roots = 1;
    vessels.max_active_tips = 2;
    vessels.surface_min_separation_voxels = 0;
    vessels.diameter_voxels = 1.0;
    vessels.inward_speed_voxels_per_hour = 10.0;
    vessels.outward_speed_voxels_per_hour = 10.0;
    vessels.inward_max_length_voxels = 12;
    vessels.outward_max_length_voxels = 12;
    vessels.outward_external_connection_distance_voxels = 1.0;
    vessels.influence_activation = "immediate";
    vessels.influence_cutoff_radius_voxels = 3.0;
    vessels.influence_max_relief_fraction = 0.5;
    config.validate();
    return config;
}

atcg3d::Model3DConfig seed_arrival_config(std::uint64_t max_total_roots,
                                          std::uint64_t max_active_tips,
                                          int minimum_separation,
                                          std::uint64_t max_events,
                                          int threads = 1) {
    atcg3d::Model3DConfig config = test_config(threads);
    auto& vessels = config.angiogenesis;
    // Keep successive site arrivals far ahead of the first tip-growth event,
    // making these tests specifically about the global seed process.
    vessels.seed_rate_sites_per_30_days = 720000000.0;
    vessels.seed_rate_sites_per_hour = vessels.seed_rate_sites_per_30_days / 720.0;
    vessels.roots_per_event = 1;
    vessels.max_total_roots = max_total_roots;
    vessels.max_active_tips = max_active_tips;
    vessels.surface_min_separation_voxels = minimum_separation;
    config.max_events = max_events;
    config.validate();
    return config;
}

std::vector<atcg3d::Vec3i> root_positions(const atcg3d::Simulation3D& simulation) {
    std::vector<atcg3d::Vec3i> roots;
    for (const atcg3d::VesselNodeSlot slot : simulation.vessel_nodes().alive_slots()) {
        if (simulation.vessel_nodes().role(slot) == atcg3d::VesselBranchRole::root) {
            roots.push_back(simulation.vessel_nodes().position(slot));
        }
    }
    std::sort(roots.begin(), roots.end());
    return roots;
}

atcg3d::CellInit scheduled_cell(atcg3d::Vec3i anchor,
                                atcg3d::CellUid uid) {
    atcg3d::CellInit cell;
    cell.anchor = anchor;
    cell.uid = uid;
    cell.clone_id = static_cast<std::uint32_t>(uid);
    cell.type = atcg3d::CellType::r;
    cell.stage = atcg3d::CellStage::small;
    cell.flags = atcg3d::kNoFlags;
    cell.inherent_growth_rate = 1.0F;
    cell.density_growth_rate = 1.0F;
    cell.migration_rate = 1.0F;
    cell.normal_migration_rate = 0.25F;
    cell.division_work_remaining = 100.0F;
    cell.next_migration_time = 100.0;
    cell.next_division_time = 100.0;
    cell.event_sequence = 11;
    return cell;
}

}  // namespace

int main() {
    using namespace atcg3d;

    Model3DConfig exponential_config = test_config(1);
    exponential_config.angiogenesis.influence_profile = "exponential";
    exponential_config.angiogenesis.influence_decay_length_voxels = 2.5;
    Simulation3D exponential_selection(exponential_config);
    assert(exponential_selection.vascular_influence().profile() ==
           VascularInfluenceProfile3D::exponential);
    assert(exponential_selection.vascular_influence().decay_length_voxels() == 2.5F);

    Simulation3D first(test_config(1));
    first.run();
    assert(first.stats().angiogenesis_seed_attempts == 1);
    assert(first.stats().angiogenesis_roots == 1);
    assert(first.vessel_nodes().alive_count() > 1);
    assert(first.vessel_tips().alive_count() == 2);
    assert(first.vessel_grid().occupied_voxel_count() > 0);
    assert(first.stats().vessel_growth_commits > 0);
    assert(first.stats().vascular_displacements > 0);
    assert(first.vascular_influence().influenced_voxel_count() > 0);
    assert(!first.angiogenesis_state().eligible);  // max_total_roots reached
    for (const Vec3i site : first.vessel_grid().occupied_sites()) {
        assert(!first.grid().available(site));
        assert(first.vessel_grid().vessel_id(site) != 0);
        assert(first.vessel_grid().perfused(site));
    }
    for (const VesselNodeSlot slot : first.vessel_nodes().alive_slots()) {
        assert(first.vessel_nodes().perfused(slot));
    }
    for (const VesselTipSlot slot : first.vessel_tips().alive_slots()) {
        assert(first.vessel_tips().perfused(slot));
    }

    Simulation3D second(test_config(4));
    second.run();
    assert(second.state_checksum() == first.state_checksum());
    assert(second.stats().vascular_displacements ==
           first.stats().vascular_displacements);

    const VasculatureState3D snapshot = first.snapshot_vasculature();
    assert(snapshot.process.committed_roots == 1);
    assert(snapshot.next_vessel_id == 2);
    assert(!snapshot.nodes.empty() && snapshot.tips.size() == 2);
    assert(snapshot.perfused_vessels == std::vector<VesselId>{1});

    // The configured unit is sites/30 days: three roots require three Poisson
    // arrivals. One arrival is never multiplied by roots_per_event.
    Simulation3D multiple(seed_arrival_config(3, 6, 0, 3));
    multiple.run();
    assert(multiple.stats().angiogenesis_seed_attempts == 3);
    assert(multiple.stats().angiogenesis_roots == 3);
    assert(multiple.stats().angiogenesis_seed_rejections == 0);
    assert(multiple.angiogenesis_state().attempted_events == 3);
    assert(multiple.angiogenesis_state().committed_roots == 3);
    assert(multiple.angiogenesis_state().rejected_events == 0);
    const std::vector<Vec3i> multiple_roots = root_positions(multiple);
    assert(multiple_roots.size() == 3);
    assert(std::adjacent_find(multiple_roots.begin(), multiple_roots.end()) ==
           multiple_roots.end());
    assert(!multiple.angiogenesis_state().eligible);  // max_total_roots reached

    Model3DConfig multiple_threads_config = seed_arrival_config(3, 6, 0, 3, 4);
    Simulation3D multiple_threads(multiple_threads_config);
    multiple_threads.run();
    assert(multiple_threads.state_checksum() == multiple.state_checksum());

    // A root creates two active tips. While both remain active, the following
    // arrival is consumed but rejected by max_active_tips.
    Simulation3D tip_limited(seed_arrival_config(3, 2, 0, 2));
    tip_limited.run();
    assert(tip_limited.stats().angiogenesis_seed_attempts == 2);
    assert(tip_limited.stats().angiogenesis_roots == 1);
    assert(tip_limited.stats().angiogenesis_seed_rejections == 1);
    assert(tip_limited.active_vessel_tip_count() == 2);

    // The first site succeeds, but every subsequent surface site lies within
    // this deliberately oversized separation radius and is rejected.
    Simulation3D separated(seed_arrival_config(3, 6, 1000, 2));
    separated.run();
    assert(separated.stats().angiogenesis_seed_attempts == 2);
    assert(separated.stats().angiogenesis_roots == 1);
    assert(separated.stats().angiogenesis_seed_rejections == 1);
    assert(root_positions(separated).size() == 1);

    // A thick new segment can hit the side of another vessel capsule without
    // landing on an existing centreline node.  This must commit as an
    // anastomosis and stop the arriving tip; overlap with the tip's own root
    // capsule remains allowed.
    Model3DConfig side_collision_config = test_config(1);
    side_collision_config.angiogenesis.enabled = true;
    side_collision_config.end_time_hours = 2.0;
    side_collision_config.max_events = 1;
    VasculatureState3D side_collision_state;
    VesselNodeInit3D moving_root;
    moving_root.position = {0, 0, 0};
    moving_root.uid = 1;
    moving_root.vessel_id = 1;
    moving_root.role = VesselBranchRole::root;
    moving_root.diameter_voxels = 2.0F;
    side_collision_state.nodes.push_back(moving_root);

    VesselNodeInit3D side_root;
    side_root.position = {1, 1, 0};
    side_root.uid = 2;
    side_root.vessel_id = 2;
    side_root.role = VesselBranchRole::root;
    side_root.diameter_voxels = 1.0F;
    side_collision_state.nodes.push_back(side_root);

    VesselTipInit3D moving_tip;
    moving_tip.position = moving_root.position;
    moving_tip.bias_axis = {1, 0, 0};
    moving_tip.target = {10, 0, 0};
    moving_tip.uid = 1;
    moving_tip.vessel_id = moving_root.vessel_id;
    moving_tip.current_node_uid = moving_root.uid;
    moving_tip.current_node_slot = 0;
    moving_tip.role = VesselBranchRole::outward;
    moving_tip.status = VesselTipStatus::active;
    moving_tip.pending_direction = 6;  // (+1, 0, 0)
    moving_tip.diameter_voxels = 2.0F;
    moving_tip.speed_voxels_per_hour = 1.0F;
    moving_tip.max_length_voxels = 10.0F;
    moving_tip.next_growth_time = 1.0;
    moving_tip.schedule_generation = 1;
    side_collision_state.tips.push_back(moving_tip);
    side_collision_state.next_vessel_id = 3;
    side_collision_state.next_node_uid = 3;
    side_collision_state.next_tip_uid = 2;

    Simulation3D side_collision(side_collision_config);
    side_collision.restore({}, 1, {}, {}, {}, side_collision_state);
    assert(side_collision.step());
    assert(side_collision.stats().vessel_growth_commits == 1);
    assert(side_collision.stats().vessel_anastomoses == 1);
    assert(side_collision.vessel_tips().status(0) == VesselTipStatus::merged);
    assert(side_collision.vessel_nodes().alive_count() == 3);
    assert(side_collision.vessel_grid().vessel_id({1, 1, 0}) == 2);

    // Two thick segments proposed at the same time have distinct endpoints,
    // but their capsules both claim the initially empty voxel (1,1,0).  The
    // first deterministic winner commits atomically and the second proposal
    // must be rejected even though that voxel is occupied by commit time.
    Model3DConfig batch_collision_config = test_config(1);
    batch_collision_config.angiogenesis.enabled = true;
    batch_collision_config.end_time_hours = 2.0;
    batch_collision_config.max_events = 2;
    VasculatureState3D batch_collision_state;

    VesselNodeInit3D first_root;
    first_root.position = {0, 0, 0};
    first_root.uid = 1;
    first_root.vessel_id = 1;
    first_root.role = VesselBranchRole::root;
    first_root.diameter_voxels = 2.0F;
    batch_collision_state.nodes.push_back(first_root);

    VesselNodeInit3D second_root;
    second_root.position = {2, 2, 0};
    second_root.uid = 2;
    second_root.vessel_id = 2;
    second_root.role = VesselBranchRole::root;
    second_root.diameter_voxels = 2.0F;
    batch_collision_state.nodes.push_back(second_root);

    VesselTipInit3D first_tip;
    first_tip.position = first_root.position;
    first_tip.bias_axis = {1, 0, 0};
    first_tip.target = {10, 0, 0};
    first_tip.uid = 1;
    first_tip.vessel_id = first_root.vessel_id;
    first_tip.current_node_uid = first_root.uid;
    first_tip.current_node_slot = 0;
    first_tip.role = VesselBranchRole::outward;
    first_tip.status = VesselTipStatus::active;
    first_tip.pending_direction = 6;  // (+1, 0, 0)
    first_tip.diameter_voxels = 2.0F;
    first_tip.speed_voxels_per_hour = 1.0F;
    first_tip.max_length_voxels = 10.0F;
    first_tip.next_growth_time = 1.0;
    first_tip.schedule_generation = 1;
    batch_collision_state.tips.push_back(first_tip);

    VesselTipInit3D second_tip;
    second_tip.position = second_root.position;
    second_tip.bias_axis = {0, -1, 0};
    second_tip.target = {2, -8, 0};
    second_tip.uid = 2;
    second_tip.vessel_id = second_root.vessel_id;
    second_tip.current_node_uid = second_root.uid;
    second_tip.current_node_slot = 1;
    second_tip.role = VesselBranchRole::outward;
    second_tip.status = VesselTipStatus::active;
    second_tip.pending_direction = 8;  // (0, -1, 0)
    second_tip.diameter_voxels = 2.0F;
    second_tip.speed_voxels_per_hour = 1.0F;
    second_tip.max_length_voxels = 10.0F;
    second_tip.next_growth_time = 1.0;
    second_tip.schedule_generation = 1;
    batch_collision_state.tips.push_back(second_tip);
    batch_collision_state.next_vessel_id = 3;
    batch_collision_state.next_node_uid = 3;
    batch_collision_state.next_tip_uid = 3;

    Simulation3D batch_collision(batch_collision_config);
    batch_collision.restore({}, 1, {}, {}, {}, batch_collision_state);
    assert(batch_collision.step());
    assert(batch_collision.stats().vessel_growth_attempts == 2);
    assert(batch_collision.stats().vessel_growth_commits == 1);
    assert(batch_collision.stats().conflict_rejections == 1);
    assert(batch_collision.stats().vessel_anastomoses == 0);
    assert(batch_collision.vessel_nodes().alive_count() == 3);

    // Perfusing this L-shaped vessel used to refresh the entire expanded AABB.
    // The remote cell at (10,10,0) is inside that AABB but more than the
    // configured cutoff from either arm; only the near cell may be visited and
    // refreshed, otherwise migration activation would be triggered nonlocally.
    Model3DConfig local_refresh_config = test_config(1);
    local_refresh_config.end_time_hours = 2.0;
    local_refresh_config.max_events = 1;
    local_refresh_config.angiogenesis.trigger_activation_volume_voxels3 = 1.0e9;
    local_refresh_config.angiogenesis.outward_external_connection_distance_voxels = 1.0;
    local_refresh_config.angiogenesis.influence_cutoff_radius_voxels = 3.0;
    local_refresh_config.migration_activation_threshold = 0.0;
    local_refresh_config.validate();

    VasculatureState3D bent_state;
    VesselNodeInit3D bend_root;
    bend_root.position = {0, 0, 0};
    bend_root.uid = 1;
    bend_root.vessel_id = 1;
    bend_root.role = VesselBranchRole::root;
    bent_state.nodes.push_back(bend_root);

    VesselNodeInit3D bend_vertical;
    bend_vertical.position = {0, 20, 0};
    bend_vertical.uid = 2;
    bend_vertical.parent_uid = 1;
    bend_vertical.parent_node_slot = 0;
    bend_vertical.vessel_id = 1;
    bend_vertical.role = VesselBranchRole::outward;
    bent_state.nodes.push_back(bend_vertical);

    VesselNodeInit3D bend_horizontal;
    bend_horizontal.position = {20, 20, 0};
    bend_horizontal.uid = 3;
    bend_horizontal.parent_uid = 2;
    bend_horizontal.parent_node_slot = 1;
    bend_horizontal.vessel_id = 1;
    bend_horizontal.role = VesselBranchRole::outward;
    bent_state.nodes.push_back(bend_horizontal);

    VesselTipInit3D connection_tip;
    connection_tip.position = bend_horizontal.position;
    connection_tip.bias_axis = {1, 0, 0};
    connection_tip.target = {40, 20, 0};
    connection_tip.uid = 1;
    connection_tip.vessel_id = 1;
    connection_tip.current_node_uid = bend_horizontal.uid;
    connection_tip.current_node_slot = 2;
    connection_tip.role = VesselBranchRole::outward;
    connection_tip.status = VesselTipStatus::active;
    connection_tip.pending_direction = 6;
    connection_tip.diameter_voxels = 1.0F;
    connection_tip.speed_voxels_per_hour = 1.0F;
    connection_tip.max_length_voxels = 100.0F;
    connection_tip.next_growth_time = 1.0;
    connection_tip.schedule_generation = 1;
    bent_state.tips.push_back(connection_tip);
    bent_state.next_vessel_id = 2;
    bent_state.next_node_uid = 4;
    bent_state.next_tip_uid = 2;

    const CellInit near_cell = scheduled_cell({2, 2, 0}, 101);
    const CellInit remote_aabb_cell = scheduled_cell({10, 10, 0}, 102);
    Simulation3D local_refresh(local_refresh_config);
    local_refresh.restore(
        {near_cell, remote_aabb_cell}, 103, {}, {}, {}, bent_state);
    assert(local_refresh.step());
    assert(local_refresh.vascular_influence().relief(near_cell.anchor) > 0.0F);
    assert(local_refresh.vascular_influence().relief(remote_aabb_cell.anchor) == 0.0F);
    const VascularRefreshDiagnostics3D diagnostics =
        local_refresh.vascular_refresh_diagnostics();
    assert(diagnostics.queried_density_blocks > 0);
    assert(diagnostics.visited_density_slots == 1);
    assert(diagnostics.refreshed_cell_slots == 1);
    assert((local_refresh.cells().flags(0) &
            static_cast<std::uint8_t>(kMigrationActive)) != 0);
    assert(local_refresh.cells().last_update_time(0) == 1.0);
    assert(local_refresh.cells().event_sequence(0) > near_cell.event_sequence);
    assert((local_refresh.cells().flags(1) &
            static_cast<std::uint8_t>(kMigrationActive)) == 0);
    assert(local_refresh.cells().last_update_time(1) == 0.0);
    assert(local_refresh.cells().event_sequence(1) == remote_aabb_cell.event_sequence);

    return 0;
}
