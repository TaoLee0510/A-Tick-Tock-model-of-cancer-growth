#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "config/model_config.hpp"
#include "core/cell_store.hpp"
#include "space/chunk_grid.hpp"
#include "vasculature/vessel_grid.hpp"

int main() {
    using namespace atcg3d;
    static_assert(CellStore3D::logical_bytes_per_slot() <= 100);
    static_assert(CellStore3D::logical_bytes_per_slot() == 99);

    CellStore3D cells;
    cells.reserve(10);
    CellInit first;
    first.uid = 1;
    first.anchor = {-1, -1, -1};
    const Slot first_slot = cells.create(first);
    assert(cells.alive_count() == 1);
    assert(cells.stage_count(CellStage::small) == 1);
    assert(cells.stage_count(CellStage::large) == 0);
    cells.set_stage(first_slot, CellStage::large);
    cells.set_stage(first_slot, CellStage::large);
    assert(cells.stage_count(CellStage::small) == 0);
    assert(cells.stage_count(CellStage::large) == 1);
    cells.erase(first_slot);
    assert(cells.stage_count(CellStage::large) == 0);
    CellInit second;
    second.uid = 2;
    second.anchor = {4, 4, 4};
    second.division_work_remaining = 12.5F;
    const Slot second_slot = cells.create(second);
    assert(second_slot == first_slot);
    assert(cells.uid(second_slot) == 2);
    assert(cells.division_work_remaining(second_slot) == 12.5F);
    assert(cells.stage_count(CellStage::small) == 1);
    assert(cells.bump_migration_schedule_generation(second_slot) == 1);
    assert(cells.migration_schedule_generation(second_slot) == 1);
    assert(cells.division_schedule_generation(second_slot) == 0);
    assert(cells.death_schedule_generation(second_slot) == 0);
    cells.set_division_schedule_generation(second_slot, 7);
    cells.set_death_schedule_generation(second_slot, 11);
    const CellInit generations = cells.snapshot(second_slot);
    assert(generations.migration_schedule_generation == 1);
    assert(generations.division_schedule_generation == 7);
    assert(generations.death_schedule_generation == 11);
    assert(generations.schedule_generation == 11);

    const float lower_time = 1.0F;
    const float upper_time = std::nextafter(
        lower_time, std::numeric_limits<float>::infinity());
    const double non_float_time = static_cast<double>(lower_time) +
        0.75 * (static_cast<double>(upper_time) - lower_time);
    assert(cells.set_next_migration_time(second_slot, non_float_time) ==
           static_cast<double>(upper_time));
    assert(cells.next_migration_time(second_slot) >= non_float_time);
    assert(cells.set_next_division_time(second_slot, non_float_time) ==
           static_cast<double>(upper_time));
    assert(cells.set_migration_activation_end_time(
               second_slot, non_float_time) == static_cast<double>(upper_time));
    // Observation times use nearest rounding, and the returned value is the
    // actual stored value used by lazy progress integration.
    assert(cells.set_last_update_time(second_slot, non_float_time) ==
           static_cast<double>(upper_time));
    assert(cells.snapshot(second_slot).last_update_time ==
           cells.last_update_time(second_slot));
    const double integration_start = cells.last_update_time(second_slot);
    double integrated_elapsed = 0.0;
    double previous_stored_time = integration_start;
    for (int step = 1; step <= 100; ++step) {
        const double requested = non_float_time +
            static_cast<double>(step) * 1.0e-9;
        const double stored =
            cells.set_last_update_time(second_slot, requested);
        integrated_elapsed += std::max(0.0, stored - previous_stored_time);
        previous_stored_time = stored;
    }
    assert(integrated_elapsed ==
           cells.last_update_time(second_slot) - integration_start);

    const auto expect_invalid_time = [](auto&& operation) {
        bool rejected = false;
        try {
            operation();
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        assert(rejected);
    };
    const double too_large =
        static_cast<double>(std::numeric_limits<float>::max()) * 2.0;
    expect_invalid_time([&] {
        cells.set_next_migration_time(
            second_slot, std::numeric_limits<double>::infinity());
    });
    expect_invalid_time([&] {
        cells.set_next_division_time(second_slot, too_large);
    });
    expect_invalid_time([&] {
        cells.set_last_update_time(
            second_slot, std::numeric_limits<double>::quiet_NaN());
    });

    CellStore3D invalid_assignment;
    CellInit invalid_cell;
    invalid_cell.uid = 99;
    invalid_cell.migration_activation_end_time = too_large;
    expect_invalid_time([&] {
        invalid_assignment.restore_layout(1, {0}, {invalid_cell}, {});
    });
    assert(invalid_assignment.slot_count() == 0);

    CellStore3D journal_store;
    CellInit journal_a;
    journal_a.uid = 100;
    const Slot journal_slot_a = journal_store.create(journal_a);
    CellInit journal_b;
    journal_b.uid = 101;
    const Slot journal_slot_b = journal_store.create(journal_b);
    journal_store.reset_checkpoint_journal();
    journal_store.set_anchor(journal_slot_a, {2, 3, 4});
    journal_store.set_density_growth_rate(journal_slot_a, 0.5F);
    journal_store.erase(journal_slot_b);
    CheckpointCellJournal3D first_journal =
        journal_store.take_checkpoint_journal();
    assert(first_journal.slot_count == 2);
    assert(first_journal.mutations.size() == 2);
    assert(first_journal.free_list_mutations.size() == 1);
    assert(first_journal.free_list_mutations.front().kind ==
           FreeListMutationKind3D::push);
    assert(first_journal.free_list_mutations.front().slot ==
           journal_slot_b);
    assert(journal_store.take_checkpoint_journal().mutations.empty());

    CellInit journal_reuse;
    journal_reuse.uid = 102;
    assert(journal_store.create(journal_reuse) == journal_slot_b);
    CheckpointCellJournal3D reuse_journal =
        journal_store.take_checkpoint_journal();
    assert(reuse_journal.mutations.size() == 1);
    assert(reuse_journal.mutations.front().alive);
    assert(reuse_journal.mutations.front().cell.uid == 102);
    assert(reuse_journal.free_list_mutations.size() == 1);
    assert(reuse_journal.free_list_mutations.front().kind ==
           FreeListMutationKind3D::pop);

    CellStore3D stable_rejection;
    CellInit valid_cell;
    valid_cell.uid = 90;
    const Slot rejected_slot = stable_rejection.create(valid_cell);
    stable_rejection.erase(rejected_slot);
    const std::vector<Slot> free_before_rejection = stable_rejection.free_slots();
    invalid_cell.next_migration_time = too_large;
    invalid_cell.migration_activation_end_time = 0.0;
    expect_invalid_time([&] { (void)stable_rejection.create(invalid_cell); });
    assert(stable_rejection.free_slots() == free_before_rejection);

    Model3DConfig config;
    config.chunk_edge = 4;
    DomainPolicy domain(config);
    SparseChunkGrid3D grid(config.chunk_edge, domain);
    assert(grid.place_single({-1, -1, -1}, second_slot));
    assert(grid.owner({-1, -1, -1}) == second_slot);

    CellInit third;
    third.uid = 3;
    third.anchor = {-1, -1, -1};
    const Slot third_slot = cells.create(third);
    assert(cells.stage_count(CellStage::small) == 2);
    assert(grid.add_colocated({-1, -1, -1}, third_slot));
    assert(grid.occupants({-1, -1, -1}).size() == 2);
    assert(grid.remove({-1, -1, -1}, second_slot));
    assert(grid.owner({-1, -1, -1}) == third_slot);

    CellInit large;
    large.uid = 4;
    large.stage = CellStage::large;
    large.anchor = {7, 7, 7};
    const Slot large_slot = cells.create(large);
    assert(cells.stage_count(CellStage::large) == 1);
    assert(grid.place_large(large.anchor, large_slot));
    assert(grid.can_move_large(large.anchor, {1, 0, 0}));
    assert(!grid.can_place_large(large.anchor));

    SparseVesselGrid3D vessels(config.chunk_edge, DomainPolicy(config));
    SparseChunkGrid3D cells_with_vessels(config.chunk_edge, DomainPolicy(config));
    cells_with_vessels.attach_vessel_grid(&vessels);
    assert(vessels.add({-5, -1, 2}, VesselBranchRole::outward).placed);
    assert(cells_with_vessels.empty({-5, -1, 2}));
    assert(cells_with_vessels.blocked_by_vessel({-5, -1, 2}));
    assert(!cells_with_vessels.available({-5, -1, 2}));
    assert(!cells_with_vessels.place_single({-5, -1, 2}, 99));
    assert(!cells_with_vessels.add_colocated({-5, -1, 2}, 99));

    // A restored checkpoint must preserve stable slots and the exact LIFO
    // free-list order, otherwise future creates diverge after resume.
    CellStore3D fragmented;
    std::vector<Slot> original_slots;
    for (std::uint64_t uid = 10; uid < 15; ++uid) {
        CellInit cell;
        cell.uid = uid;
        cell.anchor = {static_cast<std::int32_t>(uid), 0, 0};
        cell.stage = uid == 12 ? CellStage::large : CellStage::small;
        original_slots.push_back(fragmented.create(cell));
    }
    fragmented.erase(original_slots[1]);
    fragmented.erase(original_slots[3]);

    const std::vector<Slot> live_layout_slots = fragmented.alive_slots();
    std::vector<CellInit> live_layout_cells;
    live_layout_cells.reserve(live_layout_slots.size());
    for (const Slot slot : live_layout_slots) {
        live_layout_cells.push_back(fragmented.snapshot(slot));
    }
    const std::vector<Slot> free_layout_slots = fragmented.free_slots();

    CellStore3D restored;
    restored.restore_layout(fragmented.slot_count(), live_layout_slots,
                            live_layout_cells, free_layout_slots);
    assert(restored.slot_count() == fragmented.slot_count());
    assert(restored.alive_slots() == live_layout_slots);
    assert(restored.free_slots() == free_layout_slots);
    assert(restored.stage_count(CellStage::large) == 1);
    assert(restored.stage_count(CellStage::small) == 2);
    for (const Slot slot : live_layout_slots) {
        assert(restored.uid(slot) == fragmented.uid(slot));
        assert(restored.anchor(slot) == fragmented.anchor(slot));
    }

    CellInit replacement;
    replacement.uid = 20;
    const Slot original_reuse_1 = fragmented.create(replacement);
    const Slot restored_reuse_1 = restored.create(replacement);
    assert(original_reuse_1 == original_slots[3]);
    assert(restored_reuse_1 == original_reuse_1);
    replacement.uid = 21;
    const Slot original_reuse_2 = fragmented.create(replacement);
    const Slot restored_reuse_2 = restored.create(replacement);
    assert(original_reuse_2 == original_slots[1]);
    assert(restored_reuse_2 == original_reuse_2);
}
