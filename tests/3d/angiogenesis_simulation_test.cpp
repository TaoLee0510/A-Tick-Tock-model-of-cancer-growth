#include <algorithm>
#include <cassert>
#include <cstdint>
#include <set>
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
    config.initial_K_migration_rate = 0.0;
    config.end_time_hours = 2.0;
    config.max_events = 100000;
    config.threads = threads;
    config.parallel_min_events_per_thread = 1;
    config.parallel_thread_thresholds = {{0, 1.0}};

    auto& vessels = config.angiogenesis;
    vessels.enabled = true;
    // Most fixtures isolate geometry/scheduling behavior. Dedicated tests
    // below exercise the density-modulated hazard independently.
    vessels.seed_process_model = "homogeneous_poisson";
    // Small deterministic fixtures need a permissive lesion detector.  The
    // production profile intentionally requires a substantially denser core.
    vessels.lesion_block_edge = 8;
    vessels.lesion_core_activation_occupied_fraction = 1.0 / 512.0;
    vessels.lesion_core_deactivation_occupied_fraction = 0.0;
    vessels.lesion_minimum_cells_per_core_block = 1;
    vessels.lesion_minimum_biological_volume_per_core_block = 0.0;
    vessels.trigger_minimum_core_blocks = 1;
    vessels.trigger_activation_volume_voxels3 = 1.0;
    vessels.trigger_deactivation_volume_voxels3 = 0.0;
    vessels.seed_rate_sites_per_30_days = 720000.0;
    vessels.seed_rate_sites_per_hour = vessels.seed_rate_sites_per_30_days / 720.0;
    vessels.max_total_roots = 1;
    vessels.max_active_tips = 2;
    vessels.surface_min_separation_voxels = 0;
    vessels.diameter_voxels = 1.0;
    vessels.inward_speed_voxels_per_hour = 10.0;
    vessels.outward_speed_voxels_per_hour = 20.0;
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

std::vector<atcg3d::CellInit> two_by_two_by_two_cluster(
    atcg3d::Vec3i origin, atcg3d::CellUid& next_uid) {
    std::vector<atcg3d::CellInit> cells;
    cells.reserve(8);
    for (int z = 0; z < 2; ++z) {
        for (int y = 0; y < 2; ++y) {
            for (int x = 0; x < 2; ++x) {
                cells.push_back(scheduled_cell(
                    origin + atcg3d::Vec3i{x, y, z}, next_uid++));
            }
        }
    }
    return cells;
}

atcg3d::Model3DConfig explicit_lesion_config(
    std::uint64_t roots, std::uint64_t events, int threads = 1) {
    atcg3d::Model3DConfig config = test_config(threads);
    config.initial_r_cells = 0;
    config.initial_K_cells = 0;
    config.max_events = events;
    config.end_time_hours = 2.0;
    auto& vessels = config.angiogenesis;
    vessels.lesion_block_edge = 1;
    vessels.lesion_core_activation_occupied_fraction = 1.0;
    vessels.lesion_core_deactivation_occupied_fraction = 1.0;
    vessels.lesion_minimum_cells_per_core_block = 1;
    vessels.lesion_halo_blocks = 0;
    vessels.trigger_activation_volume_voxels3 = 4.0;
    vessels.trigger_deactivation_volume_voxels3 = 0.0;
    vessels.trigger_minimum_core_blocks = 1;
    vessels.seed_rate_sites_per_30_days = 720000000.0;
    vessels.seed_rate_sites_per_hour =
        vessels.seed_rate_sites_per_30_days / 720.0;
    vessels.max_total_roots = roots;
    vessels.max_roots_per_lesion = 1;
    vessels.max_active_tips = 2 * roots;
    vessels.max_active_tips_per_lesion = 2;
    vessels.surface_min_separation_voxels = 0;
    config.validate();
    return config;
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
    assert(snapshot.lesions.processes.size() == 1);
    assert(snapshot.next_vessel_id == 2);
    assert(!snapshot.nodes.empty() && snapshot.tips.size() == 2);
    assert(snapshot.perfused_vessels == std::vector<VesselId>{1});
    const LesionId first_source = snapshot.nodes.front().source_lesion_id;
    assert(first_source != kNoLesionId);
    for (const VesselNodeInit3D& node : snapshot.nodes) {
        assert(node.source_lesion_id == first_source);
    }
    for (const VesselTipInit3D& tip : snapshot.tips) {
        assert(tip.source_lesion_id == first_source);
    }

    // A disconnected outlier or small metastatic focus is indexed but does
    // not share the primary lesion's angiogenesis clock.  Until it crosses its
    // own configured biological-volume threshold, only the mature lesion can
    // create a root, and that root replaces a cell at the source lesion's
    // inside surface voxel.
    CellUid next_lesion_uid = 1000;
    std::vector<CellInit> primary_cells =
        two_by_two_by_two_cluster({0, 0, 0}, next_lesion_uid);
    primary_cells.push_back(scheduled_cell({30, 0, 0}, next_lesion_uid++));
    primary_cells.push_back(scheduled_cell({31, 0, 0}, next_lesion_uid++));
    Simulation3D primary_only(explicit_lesion_config(1, 1));
    primary_only.restore(primary_cells, next_lesion_uid, {}, {}, {});
    assert(primary_only.lesion_index().lesions().size() == 2);
    const LesionId primary_id =
        *primary_only.lesion_index().lesion_for_site({0, 0, 0});
    const LesionId small_metastasis_id =
        *primary_only.lesion_index().lesion_for_site({30, 0, 0});
    assert(primary_id != small_metastasis_id);
    primary_only.run();
    assert(primary_only.stats().angiogenesis_roots == 1);
    const std::vector<Vec3i> primary_roots = root_positions(primary_only);
    assert(primary_roots.size() == 1);
    assert(primary_roots.front().x >= 0 && primary_roots.front().x <= 1);
    assert(primary_roots.front().y >= 0 && primary_roots.front().y <= 1);
    assert(primary_roots.front().z >= 0 && primary_roots.front().z <= 1);
    const VasculatureState3D primary_state =
        primary_only.snapshot_vasculature();
    assert(primary_state.lesions.next_refresh_time_hours >
           primary_only.clock().time_hours);
    assert(primary_state.lesions.next_refresh_time_hours <=
           primary_only.clock().time_hours +
               primary_only.config().angiogenesis.lesion_refresh_interval_hours +
               1e-12);
    assert(primary_state.lesions.refresh_schedule_generation > 0);
    assert(primary_state.nodes.front().source_lesion_id == primary_id);
    const auto small_process = std::find_if(
        primary_state.lesions.processes.begin(),
        primary_state.lesions.processes.end(),
        [small_metastasis_id](const LesionAngiogenesisState3D& state) {
            return state.lesion_id == small_metastasis_id;
        });
    // Ineligible lesions no longer allocate empty Poisson-process state. The
    // process is created deterministically when this focus reaches its own
    // activation threshold.
    assert(small_process == primary_state.lesions.processes.end());
    assert(primary_state.lesions.processes.size() == 1);

    // Once two spatially disconnected lesions independently cross the same
    // threshold, each owns a distinct Poisson process and creates at most one
    // root from its own surface.  No centre-primary special case is used.
    CellUid next_two_lesion_uid = 2000;
    std::vector<CellInit> two_lesion_cells =
        two_by_two_by_two_cluster({0, 0, 0}, next_two_lesion_uid);
    std::vector<CellInit> distant_cells =
        two_by_two_by_two_cluster({30, 0, 0}, next_two_lesion_uid);
    two_lesion_cells.insert(two_lesion_cells.end(),
                            distant_cells.begin(), distant_cells.end());
    Simulation3D two_lesions(explicit_lesion_config(2, 2));
    two_lesions.restore(two_lesion_cells, next_two_lesion_uid, {}, {}, {});
    assert(two_lesions.lesion_index().lesions().size() == 2);
    two_lesions.run();
    assert(two_lesions.stats().angiogenesis_seed_attempts == 2);
    assert(two_lesions.stats().angiogenesis_roots == 2);
    std::set<LesionId> source_lesions;
    for (const VesselNodeSlot slot : two_lesions.vessel_nodes().alive_slots()) {
        if (two_lesions.vessel_nodes().role(slot) == VesselBranchRole::root) {
            source_lesions.insert(
                two_lesions.vessel_nodes().source_lesion_id(slot));
        }
    }
    assert(source_lesions.size() == 2);
    Simulation3D two_lesions_parallel(explicit_lesion_config(2, 2, 4));
    two_lesions_parallel.restore(
        two_lesion_cells, next_two_lesion_uid, {}, {}, {});
    two_lesions_parallel.run();
    assert(two_lesions_parallel.state_checksum() ==
           two_lesions.state_checksum());

    // One migration can split lesion A while its detached fragment merges
    // into lesion B in the same coarse-topology refresh.  Both A and B remain
    // current IDs, so neither independent Poisson process may be erased or
    // merged into the other merely because A appears as B's predecessor.
    Model3DConfig split_merge_config = explicit_lesion_config(10, 10);
    split_merge_config.thin_layer = true;
    split_merge_config.bounded_domain = true;
    split_merge_config.domain_policy = "bounded";
    split_merge_config.domain_min = {0, 0, 0};
    split_merge_config.domain_max = {4, 0, 0};
    split_merge_config.angiogenesis.lesion_refresh_interval_hours = 0.15;
    split_merge_config.angiogenesis.trigger_activation_volume_voxels3 = 1.0;
    split_merge_config.angiogenesis.seed_rate_sites_per_30_days = 1.0e-9;
    split_merge_config.angiogenesis.seed_rate_sites_per_hour =
        split_merge_config.angiogenesis.seed_rate_sites_per_30_days / 720.0;
    split_merge_config.end_time_hours = 1.0;
    split_merge_config.max_events = 2;
    split_merge_config.validate();

    std::vector<CellInit> split_merge_cells;
    for (int x = 0; x <= 2; ++x) {
        split_merge_cells.push_back(scheduled_cell({x, 0, 0}, 4000 + x));
    }
    split_merge_cells.push_back(scheduled_cell({4, 0, 0}, 4003));
    CellInit& bridge_cell = split_merge_cells[2];
    bridge_cell.type = CellType::K;
    bridge_cell.migration_rate = 10.0F;
    bridge_cell.normal_migration_rate = 10.0F;
    bridge_cell.next_migration_time = 0.1;

    Simulation3D split_merge(split_merge_config);
    split_merge.restore(split_merge_cells, 4004, {}, {}, {});
    const LesionId split_source =
        *split_merge.lesion_index().lesion_for_site({0, 0, 0});
    const LesionId merge_target =
        *split_merge.lesion_index().lesion_for_site({4, 0, 0});
    assert(split_source != merge_target);
    const VasculatureState3D before_split_merge =
        split_merge.snapshot_vasculature();
    const auto process_for = [](const VasculatureState3D& state,
                                LesionId id) -> const AngiogenesisProcessState3D& {
        const auto found = std::find_if(
            state.lesions.processes.begin(), state.lesions.processes.end(),
            [id](const LesionAngiogenesisState3D& entry) {
                return entry.lesion_id == id;
            });
        assert(found != state.lesions.processes.end());
        return found->process;
    };
    const double split_seed_time =
        process_for(before_split_merge, split_source).next_seed_time_hours;
    const double target_seed_time =
        process_for(before_split_merge, merge_target).next_seed_time_hours;

    assert(split_merge.step());  // x=2 migrates to the only feasible site x=3
    assert(split_merge.cells().anchor(2) == (Vec3i{3, 0, 0}));
    assert(split_merge.step());  // the scheduled lesion refresh at t=0.15
    assert(split_merge.clock().time_hours == 0.15);
    assert(*split_merge.lesion_index().lesion_for_site({0, 0, 0}) ==
           split_source);
    assert(*split_merge.lesion_index().lesion_for_site({3, 0, 0}) ==
           merge_target);
    const VasculatureState3D after_split_merge =
        split_merge.snapshot_vasculature();
    assert(process_for(after_split_merge, split_source).next_seed_time_hours ==
           split_seed_time);
    assert(process_for(after_split_merge, merge_target).next_seed_time_hours ==
           target_seed_time);
    assert(after_split_merge.lesions.source_ownership.empty());

    // Checkpoint-like in-memory restore between a root displacement and the
    // next coarse-lesion refresh must preserve both the pre-refresh block
    // observations and the pending refresh event.  Rebuilding directly from
    // current cells would otherwise apply the displacement too early and
    // diverge after resume.
    CellUid next_resume_uid = 3000;
    std::vector<CellInit> resume_cells =
        two_by_two_by_two_cluster({0, 0, 0}, next_resume_uid);
    Model3DConfig resume_config = explicit_lesion_config(1, 8);
    Simulation3D uninterrupted(resume_config);
    uninterrupted.restore(resume_cells, next_resume_uid, {}, {}, {});
    assert(uninterrupted.step());  // the root arrival
    const VasculatureState3D resume_vessels =
        uninterrupted.snapshot_vasculature();
    assert(!resume_vessels.lesions.dirty_blocks.empty());
    assert(resume_vessels.lesions.next_refresh_time_hours >
           uninterrupted.clock().time_hours);
    Simulation3D resumed(resume_config);
    resumed.restore(
        uninterrupted.snapshot_cells(), uninterrupted.next_uid(),
        uninterrupted.clock(), uninterrupted.stats(), uninterrupted.lineage(),
        resume_vessels, uninterrupted.cells().slot_count(),
        uninterrupted.snapshot_cell_slots(), uninterrupted.cells().free_slots());
    assert(resumed.state_checksum() == uninterrupted.state_checksum());
    uninterrupted.run();
    resumed.run();
    assert(resumed.state_checksum() == uninterrupted.state_checksum());

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

    // The production density model remains a repeated Poisson process when
    // its intensity is scaled by approximate lesion surface area. With a
    // fixed seed, three arrivals still create three distinct roots rather
    // than collapsing into a one-shot pair of tips.
    Model3DConfig surface_scaled_config =
        seed_arrival_config(3, 6, 0, 3);
    surface_scaled_config.angiogenesis.seed_process_model =
        "density_modulated_poisson_v1";
    surface_scaled_config.angiogenesis.seed_density_stress_on_fraction = 0.0;
    surface_scaled_config.angiogenesis.seed_density_stress_full_fraction =
        0.001;
    surface_scaled_config.angiogenesis.seed_volume_reference_voxels3 = 1.0;
    surface_scaled_config.angiogenesis.seed_volume_exponent = 2.0 / 3.0;
    surface_scaled_config.angiogenesis.seed_minimum_rate_multiplier = 0.25;
    surface_scaled_config.angiogenesis.seed_maximum_rate_multiplier = 64.0;
    surface_scaled_config.validate();
    Simulation3D surface_scaled(surface_scaled_config);
    surface_scaled.initialize();
    assert(surface_scaled.angiogenesis_state()
               .current_rate_sites_per_30_days >
           surface_scaled_config.angiogenesis.seed_rate_sites_per_30_days);
    surface_scaled.run();
    assert(surface_scaled.stats().angiogenesis_seed_attempts == 3);
    assert(surface_scaled.stats().angiogenesis_roots == 3);
    assert(root_positions(surface_scaled).size() == 3);

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

    // The legacy minimum path length is only a floor. A source lesion wider
    // than that value receives a deterministic scale-aware budget sufficient
    // to reach its centroid, cross its far half, and leave a configured
    // exterior margin.
    Model3DConfig scaled_length_config = test_config(1);
    scaled_length_config.initial_r_cells = 0;
    scaled_length_config.initial_K_cells = 0;
    scaled_length_config.max_events = 1;
    scaled_length_config.angiogenesis.inward_max_length_voxels = 32;
    scaled_length_config.angiogenesis.inward_length_tortuosity_factor = 1.5;
    scaled_length_config.angiogenesis.inward_exit_margin_voxels = 16.0;
    scaled_length_config.angiogenesis.inward_hard_max_length_voxels = 1024;
    scaled_length_config.validate();
    std::vector<CellInit> scaled_length_cells;
    for (int x = 0; x <= 300; ++x) {
        scaled_length_cells.push_back(
            scheduled_cell({x, 0, 0}, 8000 + x));
    }
    Simulation3D scaled_length(scaled_length_config);
    scaled_length.restore(scaled_length_cells, 8301, {}, {}, {});
    assert(scaled_length.step());
    bool found_scaled_inward = false;
    for (const VesselTipSlot slot :
         scaled_length.vessel_tips().alive_slots()) {
        if (scaled_length.vessel_tips().role(slot) !=
            VesselBranchRole::inward) {
            continue;
        }
        found_scaled_inward = true;
        assert(scaled_length.vessel_tips().max_length_voxels(slot) > 128.0F);
        assert(scaled_length.vessel_tips().max_length_voxels(slot) <=
               1024.0F);
    }
    assert(found_scaled_inward);

    // Reaching the initial centroid is a phase transition, not a terminal
    // state. An empty pocket beyond the centre must not be mistaken for the
    // far tumour surface: the inward branch crosses it, displaces cells on
    // the far side, and continues outside until its path budget is exhausted.
    Model3DConfig through_config = test_config(1);
    through_config.initial_r_cells = 0;
    through_config.initial_K_cells = 0;
    through_config.thin_layer = true;
    through_config.bounded_domain = true;
    through_config.domain_policy = "bounded";
    through_config.domain_min = {-1, 0, 0};
    through_config.domain_max = {10, 0, 0};
    through_config.end_time_hours = 30.0;
    through_config.max_events = 20;
    through_config.angiogenesis.trigger_activation_volume_voxels3 = 1.0e9;
    through_config.angiogenesis.direction_persistence_probability = 1.0;
    through_config.angiogenesis.inward_path_policy = "through_lesion_v1";
    through_config.validate();
    std::vector<CellInit> through_cells;
    for (const int x : {0, 1, 2, 4, 5}) {
        through_cells.push_back(scheduled_cell({x, 0, 0}, 5000 + x));
    }
    VasculatureState3D through_state;
    VesselNodeInit3D through_root;
    through_root.position = {-1, 0, 0};
    through_root.uid = 1;
    through_root.vessel_id = 1;
    through_root.role = VesselBranchRole::root;
    through_root.diameter_voxels = 1.0F;
    through_state.nodes.push_back(through_root);
    VesselTipInit3D through_tip;
    through_tip.position = through_root.position;
    through_tip.bias_axis = {1, 0, 0};
    through_tip.target = {2, 0, 0};
    through_tip.uid = 1;
    through_tip.vessel_id = 1;
    through_tip.current_node_uid = 1;
    through_tip.current_node_slot = 0;
    through_tip.role = VesselBranchRole::inward;
    through_tip.status = VesselTipStatus::active;
    through_tip.last_direction = 6;
    through_tip.pending_direction = 6;
    through_tip.diameter_voxels = 1.0F;
    through_tip.speed_voxels_per_hour = 1.0F;
    through_tip.max_length_voxels = 10.0F;
    through_tip.next_growth_time = 1.0;
    through_tip.schedule_generation = 1;
    through_state.tips.push_back(through_tip);
    through_state.next_vessel_id = 2;
    through_state.next_node_uid = 2;
    through_state.next_tip_uid = 2;
    Simulation3D through(through_config);
    through.restore(through_cells, 5006, {}, {}, {}, through_state);
    assert(through.step());
    assert(through.step());
    assert(through.step());
    assert(through.vessel_tips().position(0) == (Vec3i{2, 0, 0}));
    assert(through.vessel_tips().status(0) == VesselTipStatus::transiting);
    through.run();
    assert(through.vessel_tips().status(0) == VesselTipStatus::max_length);
    assert(through.vessel_tips().position(0).x > 5);
    assert(through.cells().alive_count() == 0);

    // An outward placeholder that reaches a disconnected tumour is converted
    // into an inward, cell-displacing branch targeted at that lesion. It does
    // not pass through the tumour as an outward non-replacing segment.
    Model3DConfig contact_config = through_config;
    contact_config.domain_min = {0, 0, 0};
    contact_config.max_events = 1;
    contact_config.angiogenesis.outward_other_lesion_contact_policy =
        "convert_to_inward";
    contact_config.validate();
    std::vector<CellInit> contact_cells{
        scheduled_cell({5, 0, 0}, 6000),
        scheduled_cell({6, 0, 0}, 6001)};
    VasculatureState3D contact_state;
    VesselNodeInit3D contact_root;
    contact_root.position = {4, 0, 0};
    contact_root.uid = 1;
    contact_root.vessel_id = 1;
    contact_root.role = VesselBranchRole::root;
    contact_root.diameter_voxels = 1.0F;
    contact_state.nodes.push_back(contact_root);
    VesselTipInit3D contact_tip = through_tip;
    contact_tip.position = contact_root.position;
    contact_tip.target = {8, 0, 0};
    contact_tip.role = VesselBranchRole::outward;
    contact_tip.max_length_voxels = 20.0F;
    contact_state.tips.push_back(contact_tip);
    contact_state.next_vessel_id = 2;
    contact_state.next_node_uid = 2;
    contact_state.next_tip_uid = 2;
    Simulation3D contact(contact_config);
    contact.restore(contact_cells, 6002, {}, {}, {}, contact_state);
    assert(contact.step());
    assert(contact.vessel_tips().role(0) == VesselBranchRole::inward);
    assert(contact.vessel_nodes().role(1) == VesselBranchRole::inward);
    assert(contact.vessel_tips().speed_voxels_per_hour(0) ==
           static_cast<float>(
               contact_config.angiogenesis.inward_speed_voxels_per_hour));
    assert(contact.vessel_tips().max_length_voxels(0) >=
           contact_config.angiogenesis.inward_max_length_voxels);
    assert(contact.vessel_tips().position(0) == (Vec3i{5, 0, 0}));
    assert(contact.cells().alive_count() == 1);
    assert(contact.stats().vascular_displacements == 1);

    // A temporarily blocked branch remains event-driven: its wakeup performs
    // no fake growth attempt, rechecks geometry, and schedules exactly one
    // later retry instead of becoming permanently terminal.
    Model3DConfig retry_config = test_config(1);
    retry_config.initial_r_cells = 0;
    retry_config.initial_K_cells = 0;
    retry_config.thin_layer = true;
    retry_config.bounded_domain = true;
    retry_config.domain_policy = "bounded";
    retry_config.domain_min = {0, 0, 0};
    retry_config.domain_max = {0, 0, 0};
    retry_config.end_time_hours = 1.5;
    retry_config.max_events = 1;
    retry_config.angiogenesis.trigger_activation_volume_voxels3 = 1.0e9;
    retry_config.angiogenesis.vessel_blocked_policy = "retry";
    retry_config.angiogenesis.vessel_blocked_retry_interval_hours = 1.0;
    retry_config.validate();
    VasculatureState3D retry_state;
    VesselNodeInit3D retry_root;
    retry_root.position = {0, 0, 0};
    retry_root.uid = 1;
    retry_root.vessel_id = 1;
    retry_root.role = VesselBranchRole::root;
    retry_root.diameter_voxels = 1.0F;
    retry_state.nodes.push_back(retry_root);
    VesselTipInit3D retry_tip = through_tip;
    retry_tip.position = retry_root.position;
    retry_tip.target = {8, 0, 0};
    retry_tip.role = VesselBranchRole::outward;
    retry_tip.pending_direction = kStayDirection;
    retry_tip.next_growth_time = 1.0;
    retry_state.tips.push_back(retry_tip);
    retry_state.next_vessel_id = 2;
    retry_state.next_node_uid = 2;
    retry_state.next_tip_uid = 2;
    Simulation3D retry(retry_config);
    retry.restore({}, 1, {}, {}, {}, retry_state);
    assert(retry.step());
    assert(retry.stats().vessel_growth_attempts == 0);
    assert(retry.vessel_tips().status(0) == VesselTipStatus::active);
    assert(retry.vessel_tips().pending_direction(0) == kStayDirection);
    assert(retry.vessel_tips().next_growth_time(0) == 2.0);

    return 0;
}
