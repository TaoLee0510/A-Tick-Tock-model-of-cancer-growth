#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "config/model_config.hpp"
#include "core/stateless_rng.hpp"
#include "engine/simulation.hpp"
#include "geometry/directions.hpp"
#include "rules/lifecycle.hpp"
#include "rules/migration.hpp"

int main() {
    using namespace atcg3d;
    Model3DConfig config;
    config.output_enabled = false;
    config.initial_r_cells = 6;
    config.initial_K_cells = 6;
    config.initial_radius = 8;
    config.end_time_hours = 8.0;
    config.max_events = 100000;
    config.density_block_edge = 2;
    config.threads = 1;

    Simulation3D first(config);
    first.run();
    assert(first.cells().alive_count() > 0);
    const auto checksum = first.state_checksum();

    config.threads = 4;
    config.parallel_min_events_per_thread = 1;
    config.parallel_thread_thresholds = {{0, 1.0}};
    Simulation3D second(config);
    second.run();
    assert(second.state_checksum() == checksum);
    assert(second.cells().alive_count() == first.cells().alive_count());

    // Rescheduling a blocked migration invalidates only the migration event;
    // unrelated division and death queue entries retain their generations.
    Model3DConfig isolated_config;
    isolated_config.output_enabled = false;
    isolated_config.domain_policy = "bounded";
    isolated_config.bounded_domain = true;
    isolated_config.domain_min = {0, 0, 0};
    isolated_config.domain_max = {0, 0, 0};
    isolated_config.density_block_edge = 1;
    isolated_config.migration_activation_enabled = false;
    isolated_config.end_time_hours = 2.0;
    CellInit isolated;
    isolated.uid = 1;
    isolated.anchor = {0, 0, 0};
    isolated.stage = CellStage::small;
    isolated.migration_rate = 1.0F;
    isolated.next_migration_time = 1.0;
    isolated.next_division_time = 10.0;
    Simulation3D generations(isolated_config);
    generations.restore({isolated}, 2, {}, {}, {});
    const Slot isolated_slot = generations.cells().alive_slots().front();
    assert(generations.cells().migration_schedule_generation(isolated_slot) == 1);
    assert(generations.cells().division_schedule_generation(isolated_slot) == 1);
    assert(generations.cells().death_schedule_generation(isolated_slot) == 1);
    assert(generations.step());
    assert(generations.cells().migration_schedule_generation(isolated_slot) == 2);
    assert(generations.cells().division_schedule_generation(isolated_slot) == 1);
    assert(generations.cells().death_schedule_generation(isolated_slot) == 1);

    // With activation enabled a normal cell still migrates. Active state must
    // additionally carry a finite end time, and any positive effective rate
    // must have a queued migration event.
    Model3DConfig gated_config = isolated_config;
    gated_config.migration_activation_enabled = true;
    CellInit gated_cell = isolated;
    gated_cell.normal_migration_rate = 1.0F;
    Simulation3D gated(gated_config);
    gated.restore({gated_cell}, 2, {}, {}, {});
    assert(gated.pending_event_count() == 2);  // migration and division
    assert(gated.step());
    assert(gated.stats().migration_attempts == 1);
    const auto expect_restore_rejected = [&](CellInit invalid) {
        bool rejected = false;
        try {
            Simulation3D candidate(gated_config);
            candidate.restore({invalid}, 2, {}, {}, {});
        } catch (const std::runtime_error&) {
            rejected = true;
        }
        assert(rejected);
    };
    CellInit missing_normal_migration = gated_cell;
    missing_normal_migration.next_migration_time = 0.0;
    expect_restore_rejected(missing_normal_migration);
    CellInit missing_migration = gated_cell;
    missing_migration.flags |= kMigrationActive;
    expect_restore_rejected(missing_migration);  // active without finite end
    CellInit inactive_with_end = gated_cell;
    inactive_with_end.migration_activation_end_time = 2.0;
    expect_restore_rejected(inactive_with_end);
    CellInit missing_division = gated_cell;
    missing_division.migration_rate = 0.0F;
    missing_division.next_division_time = 0.0;
    expect_restore_rejected(missing_division);  // active growth without division

    // A division-created local density threshold crossing activates migration
    // for both cells and schedules their first migration event.
    Model3DConfig crossing_config;
    crossing_config.output_enabled = false;
    crossing_config.migration_activation_enabled = true;
    crossing_config.migration_activation_window_edge = 3;
    crossing_config.migration_activation_block_edge = 1;
    crossing_config.migration_activation_threshold = 0.07;
    crossing_config.density_block_edge = 1;
    crossing_config.end_time_hours = 2.1;
    crossing_config.max_events = 1000;
    CellInit crossing_cell;
    crossing_cell.uid = 10;
    crossing_cell.type = CellType::K;
    crossing_cell.migration_rate = 1.0F;
    crossing_cell.normal_migration_rate = 1.0F;
    crossing_cell.division_work_remaining = 0.0F;
    crossing_cell.next_migration_time = 1.0;
    crossing_cell.next_division_time = 1.0;
    Simulation3D crossing(crossing_config);
    crossing.restore({crossing_cell}, 11, {}, {}, {});
    crossing.run();
    assert(crossing.stats().divisions == 1);
    assert(crossing.cells().alive_count() == 2);
    for (const Slot slot : crossing.cells().alive_slots()) {
        assert(migration_allowed_for_cell(slot, crossing.cells(), crossing_config));
        assert(crossing.cells().next_migration_time(slot) > 0.0);
    }
    assert(crossing.migration_activation_class_recomputes() > 0);

    // End-to-end activation lifecycle: falling density does not terminate an
    // active interval. Its exact end event returns the survivor to ordinary
    // migration, and a later high-density refresh may activate it again.
    Model3DConfig toggle_config;
    toggle_config.output_enabled = false;
    toggle_config.domain_policy = "bounded";
    toggle_config.bounded_domain = true;
    toggle_config.domain_min = {0, 0, 0};
    toggle_config.domain_max = {3, 3, 0};
    toggle_config.thin_layer = true;
    toggle_config.density_block_edge = 1;
    toggle_config.migration_activation_enabled = true;
    toggle_config.migration_activation_window_edge = 3;
    toggle_config.migration_activation_block_edge = 1;
    toggle_config.migration_activation_threshold = 0.15;
    toggle_config.end_time_hours = 3.1;
    toggle_config.max_events = 100;
    CellInit toggled;
    toggled.uid = 30;
    toggled.type = CellType::K;
    toggled.anchor = {1, 1, 0};
    toggled.flags = kDirtyDensity | kMigrationActive;
    toggled.migration_rate = 1.0F;
    toggled.normal_migration_rate = 1.0F;
    toggled.migration_activation_end_time = 1.0;
    toggled.next_migration_time = 1.0;
    toggled.density_growth_rate = 1.0F;
    toggled.division_work_remaining = 2.0F;
    toggled.next_division_time = 2.0;
    CellInit dying;
    dying.uid = 31;
    dying.type = CellType::K;
    dying.anchor = {2, 1, 0};
    dying.inherent_growth_rate = 0.001F;
    dying.density_growth_rate = 0.0F;
    dying.next_division_time = 0.0;
    dying.death_deadline = 0.5;
    Simulation3D toggling(toggle_config);
    toggling.restore({toggled, dying}, 32, {}, {}, {});
    assert(toggling.step());
    const Slot toggle_mother = toggling.cells().alive_slots().front();
    assert(toggling.clock().time_hours == 0.5);
    assert(toggling.cells().alive_count() == 1);
    assert(migration_allowed_for_cell(toggle_mother, toggling.cells(), toggle_config));
    assert((toggling.cells().flags(toggle_mother) & kMigrationActive) != 0);
    assert(toggling.cells().next_migration_time(toggle_mother) == 1.0);
    assert(toggling.stats().migration_attempts == 0);
    assert(toggling.migration_activation_bulk_slot_visits() > 0);
    assert(toggling.step());
    assert(toggling.clock().time_hours == 1.0);
    assert((toggling.cells().flags(toggle_mother) & kMigrationActive) == 0);
    assert(toggling.cells().migration_activation_end_time(toggle_mother) == 0.0);
    assert(toggling.cells().next_migration_time(toggle_mother) >
           toggling.clock().time_hours);
    assert(toggling.step());
    assert(toggling.stats().divisions == 1);
    assert(toggling.cells().alive_count() == 2);

    // Moving from an active query block into an inactive one must refresh the
    // moved cell even though neither block crosses its aggregate threshold.
    // More than 400 source-block residents prove the unchanged class does not
    // trigger the former whole-block scan.
    Model3DConfig frontier_config;
    frontier_config.output_enabled = false;
    frontier_config.thin_layer = true;
    frontier_config.density_block_edge = 4;
    frontier_config.growth_density_window_edge = 2;
    frontier_config.migration_activation_enabled = true;
    frontier_config.migration_activation_window_edge = 70;
    frontier_config.migration_activation_block_edge = 32;
    frontier_config.migration_activation_threshold = 0.001;
    frontier_config.end_time_hours = 1.1;
    frontier_config.max_events = 10;
    std::vector<CellInit> frontier_cells;
    CellUid support_uid = 1000;
    for (std::int32_t x = 0; x <= 12; ++x) {
        for (std::int32_t y = 0; y <= 30; ++y) {
            CellInit support;
            support.uid = support_uid++;
            support.type = CellType::K;
            support.anchor = {x, y, 0};
            support.flags = kDirtyDensity | kMigrationActive;
            support.migration_rate = 0.0F;
            support.migration_activation_end_time = 1000.0;
            support.density_growth_rate = 1.0F;
            support.division_work_remaining = 1000.0F;
            support.next_division_time = 1000.0;
            frontier_cells.push_back(support);
        }
    }
    CellUid moving_uid = support_uid;
    while (true) {
        const auto index = rng_bounded(
            frontier_config.seed, moving_uid,
            static_cast<std::uint64_t>(RngEventKind::migration_direction),
            0, 8);
        const DirectionId direction = static_cast<DirectionId>(index + 1);
        if (direction_vector(direction).x == 1) break;
        ++moving_uid;
    }
    CellInit frontier_moving;
    frontier_moving.uid = moving_uid;
    frontier_moving.type = CellType::K;
    frontier_moving.stage = CellStage::large;
    frontier_moving.anchor = {31, 0, 0};
    frontier_moving.flags = kDirtyDensity | kMigrationActive;
    frontier_moving.migration_rate = 1.0F;
    frontier_moving.migration_activation_end_time = 1000.0;
    frontier_moving.next_migration_time = 1.0;
    frontier_moving.density_growth_rate = 1.0F;
    frontier_moving.division_work_remaining = 1000.0F;
    frontier_moving.next_division_time = 1000.0;
    frontier_cells.push_back(frontier_moving);
    Simulation3D frontier(frontier_config);
    frontier.restore(frontier_cells, moving_uid + 1, {}, {}, {});
    assert(frontier.migration_activation_class_cache_size() == 1);
    assert(frontier.migration_activation_bulk_slot_visits() == 0);
    assert(frontier.step());
    Slot moved_slot = kEmptySlot;
    for (const Slot slot : frontier.cells().alive_slots()) {
        if (frontier.cells().uid(slot) == moving_uid) moved_slot = slot;
    }
    assert(moved_slot != kEmptySlot);
    assert(frontier.cells().anchor(moved_slot).x == 32);
    assert(migration_allowed_for_cell(
        moved_slot, frontier.cells(), frontier_config));
    assert((frontier.cells().flags(moved_slot) & kMigrationActive) != 0);
    assert(frontier.cells().next_migration_time(moved_slot) >
           frontier.clock().time_hours);
    assert(frontier.migration_activation_class_cache_size() == 2);
    assert(frontier.migration_activation_bulk_slot_visits() == 0);
    assert(frontier.migration_activation_direct_slot_visits() == 1);
    assert(frontier.migration_activation_class_recomputes() <= 2);

    // A stage-2 survivor can recover to a large anchor in a different query
    // block after its co-located neighbor dies. The recovered before/after
    // anchors are appended to dirty sites, materializing and directly visiting
    // the new block even when the original death site is now empty.
    Model3DConfig recovery_config;
    recovery_config.output_enabled = false;
    recovery_config.domain_policy = "bounded";
    recovery_config.bounded_domain = true;
    recovery_config.domain_min = {0, 0, 0};
    recovery_config.domain_max = {1, 1, 1};
    recovery_config.density_block_edge = 1;
    recovery_config.migration_activation_enabled = true;
    recovery_config.migration_activation_window_edge = 3;
    recovery_config.migration_activation_block_edge = 1;
    recovery_config.migration_activation_threshold = 0.07;
    recovery_config.end_time_hours = 1.0;
    recovery_config.max_events = 10;
    CellInit survivor;
    survivor.uid = 2000;
    survivor.type = CellType::K;
    survivor.stage = CellStage::ultrasmall;
    survivor.anchor = {1, 1, 1};
    survivor.flags = kDirtyDensity | kMigrationActive;
    survivor.migration_rate = 0.0F;
    survivor.migration_activation_end_time = 1000.0;
    survivor.density_growth_rate = 1.0F;
    survivor.division_work_remaining = 10.0F;
    survivor.next_division_time = 10.0;
    CellInit recovery_dying = survivor;
    recovery_dying.uid = 2001;
    recovery_dying.inherent_growth_rate = 0.001F;
    recovery_dying.density_growth_rate = 0.0F;
    recovery_dying.division_work_remaining = 0.0F;
    recovery_dying.next_division_time = 0.0;
    recovery_dying.death_deadline = 0.5;
    Simulation3D recovery(recovery_config);
    recovery.restore({survivor, recovery_dying}, 2002, {}, {}, {});
    assert(recovery.migration_activation_class_cache_size() == 1);
    assert(recovery.step());
    assert(recovery.cells().alive_count() == 1);
    const Slot recovered_slot = recovery.cells().alive_slots().front();
    assert(recovery.cells().stage(recovered_slot) == CellStage::large);
    assert(recovery.cells().anchor(recovered_slot) == Vec3i(0, 0, 0));
    assert(recovery.migration_activation_class_cache_size() == 2);
    assert(recovery.migration_activation_direct_slot_visits() >= 1);
    assert(migration_allowed_for_cell(
        recovered_slot, recovery.cells(), recovery_config));

    // Frequent migration refreshes integrate one persistent cycle; they do not
    // restart it indefinitely or accumulate an unbounded queue of stale growth
    // events.
    Model3DConfig progress_config;
    progress_config.output_enabled = false;
    progress_config.migration_activation_enabled = false;
    progress_config.density_block_edge = 1;
    progress_config.end_time_hours = 3.0;
    progress_config.max_events = 10000;
    progress_config.division_timing.base_cycle_hours = 2.0;
    progress_config.division_timing.minimum_fraction = 1.0;
    progress_config.division_timing.stochastic_tail_fraction = 0.0;
    CellInit moving;
    moving.uid = 20;
    moving.type = CellType::K;
    moving.migration_rate = 4.0F;
    moving.density_growth_rate = 1.0F;
    moving.division_work_remaining = 2.0F;
    moving.next_migration_time = 0.25;
    moving.next_division_time = 2.0;
    Simulation3D progressing(progress_config);
    progressing.restore({moving}, 21, {}, {}, {});
    progressing.run();
    assert(progressing.stats().migration_attempts >= 4);
    assert(progressing.stats().divisions >= 1);
    assert(progressing.pending_event_count() < 20);

    // Repeated local density changes can move many future division times. Old
    // generation entries are inert, and deterministic threshold compaction
    // keeps their heap storage bounded during a long migration-heavy run.
    Model3DConfig queue_config;
    queue_config.output_enabled = false;
    queue_config.domain_policy = "bounded";
    queue_config.bounded_domain = true;
    queue_config.domain_min = {0, 0, 0};
    queue_config.domain_max = {5, 5, 5};
    queue_config.migration_activation_enabled = false;
    queue_config.density_block_edge = 1;
    queue_config.growth_density_window_edge = 2;
    queue_config.K_limit = 4.0;
    queue_config.carrying_capacity_K = 4.0;
    queue_config.legacy_mapping.source_K_limit = 4.0 / 6.0;
    queue_config.legacy_mapping.source_carrying_capacity_K = 4.0 / 6.0;
    queue_config.K_death_delay_hours = 1000000.0;
    queue_config.end_time_hours = 250.0;
    queue_config.max_events = 30000;
    std::vector<CellInit> queue_cells(100);
    for (std::size_t index = 0; index < queue_cells.size(); ++index) {
        CellInit& cell = queue_cells[index];
        cell.uid = 1000 + index;
        cell.type = CellType::K;
        cell.anchor = {static_cast<std::int32_t>(index % 5),
                       static_cast<std::int32_t>((index / 5) % 5),
                       static_cast<std::int32_t>(index / 25)};
        cell.migration_rate = 1.0F;
        cell.next_migration_time = 1.0;
        cell.density_growth_rate = 1.0F;
        cell.division_work_remaining = 1000.0F;
        cell.next_division_time = 1000.0;
    }
    Simulation3D bounded_queue(queue_config);
    bounded_queue.restore(queue_cells, 1100, {}, {}, {});
    bounded_queue.run();
    assert(bounded_queue.clock().completed_events >= 20000);
    assert(bounded_queue.event_queue_rebuild_count() > 0);
    assert(bounded_queue.pending_event_count() < 500);

    // Same-time divisions propose against one occupancy snapshot. Both K
    // mothers have exactly one possible daughter voxel, so the higher stable
    // conflict hash must win even when that UID is numerically larger (the old
    // event-heap order would always let the smaller UID win). The loser only
    // retries and does not fall through to co-location.
    Model3DConfig division_conflict_config;
    division_conflict_config.output_enabled = false;
    division_conflict_config.domain_policy = "bounded";
    division_conflict_config.bounded_domain = true;
    division_conflict_config.domain_min = {0, 0, 0};
    division_conflict_config.domain_max = {2, 0, 0};
    division_conflict_config.thin_layer = true;
    division_conflict_config.density_block_edge = 1;
    division_conflict_config.migration_activation_enabled = false;
    division_conflict_config.ultrasmall_enabled = false;
    division_conflict_config.end_time_hours = 1.0;
    division_conflict_config.max_events = 2;
    CellUid lower_uid = 100;
    while (division_conflict_priority(
               division_conflict_config, 1.0, lower_uid + 1) <=
           division_conflict_priority(
               division_conflict_config, 1.0, lower_uid)) {
        lower_uid += 2;
    }
    const CellUid higher_uid = lower_uid + 1;
    assert(division_conflict_priority(
               division_conflict_config, 1.0, higher_uid) >
           division_conflict_priority(
               division_conflict_config, 1.0, lower_uid));

    const auto run_division_conflict = [&](int threads, bool reverse_uid_positions) {
        Model3DConfig local = division_conflict_config;
        local.threads = threads;
        CellInit left;
        left.uid = reverse_uid_positions ? higher_uid : lower_uid;
        left.type = CellType::K;
        left.stage = CellStage::small;
        left.anchor = {0, 0, 0};
        left.density_growth_rate = 1.0F;
        left.division_work_remaining = 0.0F;
        left.next_division_time = 1.0;
        CellInit right = left;
        right.uid = reverse_uid_positions ? lower_uid : higher_uid;
        right.anchor = {2, 0, 0};
        Simulation3D candidate(local);
        candidate.restore({left, right}, higher_uid + 1, {}, {}, {});
        assert(candidate.step());
        assert(candidate.stats().divisions == 1);
        assert(candidate.stats().deaths == 0);
        assert(candidate.stats().conflict_rejections == 1);
        assert(candidate.cells().alive_count() == 3);
        assert(candidate.lineage().size() == 1);
        assert(candidate.lineage().front().parent_uid == higher_uid);
        Slot loser = kEmptySlot;
        for (const Slot slot : candidate.cells().alive_slots()) {
            if (candidate.cells().uid(slot) == lower_uid) loser = slot;
        }
        assert(loser != kEmptySlot);
        assert(candidate.cells().next_division_time(loser) == 2.0);
        assert(candidate.cells().stage(loser) == CellStage::small);
        assert(candidate.cells().event_sequence(loser) == 0);
        return candidate.state_checksum();
    };
    const std::uint64_t division_conflict_checksum =
        run_division_conflict(1, false);
    assert(run_division_conflict(4, false) == division_conflict_checksum);
    (void)run_division_conflict(1, true);

    // Two scheduled stage-2 cells in one co-location group deliberately choose
    // different free targets. Their DivisionProposals still carry the same
    // source-group lock, so exactly one recovery commits and the other retries.
    Model3DConfig group_conflict_config = division_conflict_config;
    group_conflict_config.domain_min = {-1, -1, 0};
    group_conflict_config.domain_max = {1, 1, 0};
    group_conflict_config.ultrasmall_enabled = true;
    CellUid first_group_uid = 1000;
    CellUid second_group_uid = 1001;
    while (rng_bounded(
               group_conflict_config.seed, first_group_uid,
               static_cast<std::uint64_t>(RngEventKind::division_location),
               0, 8) ==
           rng_bounded(
               group_conflict_config.seed, second_group_uid,
               static_cast<std::uint64_t>(RngEventKind::division_location),
               0, 8)) {
        ++second_group_uid;
    }
    const CellUid expected_group_winner =
        division_conflict_priority(group_conflict_config, 1.0, first_group_uid) >
                division_conflict_priority(
                    group_conflict_config, 1.0, second_group_uid)
            ? first_group_uid
            : second_group_uid;
    const CellUid expected_group_loser = expected_group_winner == first_group_uid
        ? second_group_uid : first_group_uid;
    const auto run_group_conflict = [&](int threads) {
        Model3DConfig local = group_conflict_config;
        local.threads = threads;
        CellInit first_group;
        first_group.uid = first_group_uid;
        first_group.type = CellType::K;
        first_group.stage = CellStage::ultrasmall;
        first_group.anchor = {0, 0, 0};
        first_group.density_growth_rate = 1.0F;
        first_group.division_work_remaining = 0.0F;
        first_group.next_division_time = 1.0;
        CellInit second_group = first_group;
        second_group.uid = second_group_uid;
        Simulation3D candidate(local);
        candidate.restore(
            {first_group, second_group}, second_group_uid + 1, {}, {}, {});
        assert(candidate.step());
        assert(candidate.stats().divisions == 0);
        assert(candidate.stats().deaths == 0);
        assert(candidate.stats().conflict_rejections == 1);
        assert(candidate.cells().alive_count() == 2);
        Slot winner = kEmptySlot;
        Slot loser = kEmptySlot;
        for (const Slot slot : candidate.cells().alive_slots()) {
            if (candidate.cells().uid(slot) == expected_group_winner) winner = slot;
            if (candidate.cells().uid(slot) == expected_group_loser) loser = slot;
        }
        assert(winner != kEmptySlot && loser != kEmptySlot);
        assert(candidate.cells().anchor(winner) != Vec3i(0, 0, 0));
        assert(candidate.cells().anchor(loser) == Vec3i(0, 0, 0));
        assert(candidate.cells().event_sequence(winner) == 1);
        assert(candidate.cells().event_sequence(loser) == 0);
        assert(candidate.cells().next_division_time(loser) == 2.0);
        assert(candidate.cells().stage(winner) == CellStage::small);
        assert(candidate.cells().stage(loser) == CellStage::small);
        return candidate.state_checksum();
    };
    const std::uint64_t group_conflict_checksum = run_group_conflict(1);
    assert(run_group_conflict(4) == group_conflict_checksum);

    // A same-time death exposes two overlapping 2x2x2 stage-recovery
    // footprints. Both proposals are built after all deaths commit but before
    // either recovery mutates occupancy; the stable recovery hash selects one
    // complete footprint and rejects the other as a unit.
    Model3DConfig recovery_conflict_config;
    recovery_conflict_config.output_enabled = false;
    recovery_conflict_config.domain_policy = "bounded";
    recovery_conflict_config.bounded_domain = true;
    recovery_conflict_config.domain_min = {0, 0, 0};
    recovery_conflict_config.domain_max = {2, 1, 1};
    recovery_conflict_config.density_block_edge = 1;
    recovery_conflict_config.migration_activation_enabled = false;
    recovery_conflict_config.end_time_hours = 1.0;
    recovery_conflict_config.max_events = 1;
    const CellUid left_recovery_uid = 3000;
    const CellUid right_recovery_uid = 3001;
    const CellUid expected_recovery_winner =
        stage_recovery_conflict_priority(
            recovery_conflict_config, 1.0, left_recovery_uid) >
                stage_recovery_conflict_priority(
                    recovery_conflict_config, 1.0, right_recovery_uid)
            ? left_recovery_uid
            : right_recovery_uid;
    const auto run_recovery_conflict = [&](int threads) {
        Model3DConfig local = recovery_conflict_config;
        local.threads = threads;
        CellInit left;
        left.uid = left_recovery_uid;
        left.type = CellType::K;
        left.stage = CellStage::small;
        left.anchor = {0, 0, 0};
        left.density_growth_rate = 1.0F;
        left.division_work_remaining = 100.0F;
        left.next_division_time = 100.0;
        CellInit right = left;
        right.uid = right_recovery_uid;
        right.anchor = {2, 0, 0};
        CellInit dying;
        dying.uid = 3002;
        dying.type = CellType::K;
        dying.stage = CellStage::small;
        dying.anchor = {1, 0, 0};
        dying.inherent_growth_rate = 0.001F;
        dying.density_growth_rate = 0.0F;
        dying.division_work_remaining = 0.0F;
        dying.next_division_time = 0.0;
        dying.death_deadline = 1.0;
        Simulation3D candidate(local);
        candidate.restore({left, right, dying}, 3003, {}, {}, {});
        assert(candidate.step());
        assert(candidate.stats().deaths == 1);
        assert(candidate.stats().conflict_rejections == 1);
        assert(candidate.cells().alive_count() == 2);
        Slot left_slot = kEmptySlot;
        Slot right_slot = kEmptySlot;
        for (const Slot slot : candidate.cells().alive_slots()) {
            if (candidate.cells().uid(slot) == left_recovery_uid) left_slot = slot;
            if (candidate.cells().uid(slot) == right_recovery_uid) right_slot = slot;
        }
        assert(left_slot != kEmptySlot && right_slot != kEmptySlot);
        assert((candidate.cells().stage(left_slot) == CellStage::large) ==
               (expected_recovery_winner == left_recovery_uid));
        assert((candidate.cells().stage(right_slot) == CellStage::large) ==
               (expected_recovery_winner == right_recovery_uid));
        return candidate.state_checksum();
    };
    const std::uint64_t recovery_conflict_checksum =
        run_recovery_conflict(1);
    assert(run_recovery_conflict(4) == recovery_conflict_checksum);

    // A same-time proposal batch is atomic with respect to max_events: if the
    // complete conflict set does not fit, step stops before the batch rather
    // than exceeding the hard cap or changing conflict semantics.
    Model3DConfig capped_config;
    capped_config.output_enabled = false;
    capped_config.migration_activation_enabled = false;
    capped_config.end_time_hours = 2.0;
    capped_config.max_events = 2;
    std::vector<CellInit> capped_cells(3);
    for (std::size_t index = 0; index < capped_cells.size(); ++index) {
        capped_cells[index].uid = 100 + index;
        capped_cells[index].anchor = {static_cast<std::int32_t>(index * 10), 0, 0};
        capped_cells[index].migration_rate = 1.0F;
        capped_cells[index].next_migration_time = 1.0;
        capped_cells[index].next_division_time = 10.0;
    }
    Simulation3D capped(capped_config);
    capped.restore(capped_cells, 103, {}, {}, {});
    assert(!capped.step());
    assert(capped.clock().completed_events == 0);
    assert(capped.stats().migration_attempts == 0);

    capped_config.max_events = 3;
    Simulation3D exact_batch(capped_config);
    exact_batch.restore(capped_cells, 103, {}, {}, {});
    assert(exact_batch.step());
    assert(exact_batch.clock().completed_events == 3);
}
