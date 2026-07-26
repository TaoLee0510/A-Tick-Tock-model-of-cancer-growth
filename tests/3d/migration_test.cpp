#include <cassert>
#include <cstdint>

#include "config/model_config.hpp"
#include "core/cell_store.hpp"
#include "rules/migration.hpp"
#include "space/chunk_grid.hpp"
#include "space/density_index.hpp"

int main() {
    using namespace atcg3d;
    Model3DConfig config;
    config.thin_layer = true;
    config.direction_density_threshold = 1.0;
    config.continue_probability = 0.90;
    config.density_block_edge = 1;
    DomainPolicy domain(config);
    SparseChunkGrid3D grid(8, domain);
    BlockDensityIndex3D density(1);
    CellStore3D cells;
    CellInit cell;
    cell.uid = 11;
    cell.type = CellType::r;
    cell.stage = CellStage::small;
    cell.anchor = {0, 0, 0};
    cell.last_direction = 1;
    cell.flags |= kMigrationActive;
    const Slot slot = cells.create(cell);
    assert(grid.place_single(cell.anchor, slot));
    density.add(cell.anchor, cell.type, slot);

    std::uint64_t forward = 0;
    std::uint64_t side_2 = 0;
    std::uint64_t side_8 = 0;
    for (std::uint64_t sequence = 0; sequence < 20000; ++sequence) {
        const DirectionId direction = select_migration_direction(slot, cells, grid, density, config, sequence);
        if (direction == 1) ++forward;
        else if (direction == 2) ++side_2;
        else if (direction == 8) ++side_8;
        else assert(false);
    }
    assert(forward > 17500 && forward < 18500);
    assert(side_2 > 700 && side_2 < 1300);
    assert(side_8 > 700 && side_8 < 1300);

    // Blocking one of the two 45-degree turns gives the 90/10 thin-layer case.
    assert(grid.place_single(direction_vector(2), 98));
    forward = side_2 = side_8 = 0;
    for (std::uint64_t sequence = 0; sequence < 10000; ++sequence) {
        const DirectionId direction = select_migration_direction(slot, cells, grid, density, config, sequence);
        if (direction == 1) ++forward;
        else if (direction == 8) ++side_8;
        else assert(false);
    }
    assert(forward > 8700 && forward < 9300);
    assert(side_8 > 700 && side_8 < 1300);
    assert(grid.remove(direction_vector(2), 98));

    assert(grid.place_single(direction_vector(1), 99));
    forward = side_2 = side_8 = 0;
    for (std::uint64_t sequence = 0; sequence < 10000; ++sequence) {
        const DirectionId direction = select_migration_direction(slot, cells, grid, density, config, sequence);
        if (direction == 2) ++side_2;
        else if (direction == 8) ++side_8;
        else assert(false);
    }
    assert(side_2 > 4500 && side_2 < 5500);
    assert(side_8 > 4500 && side_8 < 5500);
    assert(grid.remove(direction_vector(1), 99));

    // The same r cell in normal state ignores its persistence history and is
    // uniform over all feasible lattice directions.
    cells.set_flags(slot, cells.flags(slot) &
        static_cast<std::uint8_t>(~kMigrationActive));
    std::uint64_t normal_counts[9]{};
    for (std::uint64_t sequence = 0; sequence < 8000; ++sequence) {
        ++normal_counts[select_migration_direction(
            slot, cells, grid, density, config, sequence)];
    }
    assert(normal_counts[0] == 0);
    for (DirectionId direction = 1; direction <= 8; ++direction) {
        assert(normal_counts[direction] > 800 && normal_counts[direction] < 1200);
    }

    // A K cell is uniform over the 26 lattice directions in a symmetric 3D
    // environment and identical keys reproduce identical choices.
    Model3DConfig symmetric_config;
    symmetric_config.output_enabled = false;
    symmetric_config.direction_density_threshold = 1.0;
    symmetric_config.density_block_edge = 1;
    DomainPolicy symmetric_domain(symmetric_config);
    SparseChunkGrid3D symmetric_grid(8, symmetric_domain);
    BlockDensityIndex3D symmetric_density(1);
    CellStore3D symmetric_cells;
    CellInit K_cell;
    K_cell.uid = 1234;
    K_cell.type = CellType::K;
    const Slot K_slot = symmetric_cells.create(K_cell);
    assert(symmetric_grid.place_single(K_cell.anchor, K_slot));
    symmetric_density.add(K_cell.anchor, K_cell.type, K_slot);
    std::uint64_t counts[27]{};
    for (std::uint64_t sequence = 0; sequence < 26000; ++sequence) {
        const DirectionId direction = select_migration_direction(
            K_slot, symmetric_cells, symmetric_grid, symmetric_density, symmetric_config, sequence);
        ++counts[direction];
        assert(direction == select_migration_direction(
            K_slot, symmetric_cells, symmetric_grid, symmetric_density, symmetric_config, sequence));
    }
    assert(counts[0] == 0);
    for (DirectionId direction = 1; direction <= 26; ++direction) {
        assert(counts[direction] > 800 && counts[direction] < 1200);
    }

    // Normal K migration remains enabled and uniformly random; activation is a
    // direction/rate state, not an on/off permission gate.
    symmetric_config.migration_activation_enabled = true;
    const MoveProposal inactive = make_move_proposal(
        K_slot, symmetric_cells, symmetric_grid, symmetric_density,
        symmetric_config, 1, 1);
    assert(inactive.slot == K_slot && inactive.direction != kStayDirection);
    symmetric_cells.set_flags(
        K_slot, symmetric_cells.flags(K_slot) | kMigrationActive);
    const MoveProposal active = make_move_proposal(
        K_slot, symmetric_cells, symmetric_grid, symmetric_density,
        symmetric_config, 1, 1);
    assert(active.slot == K_slot && active.direction != kStayDirection);
    // Conflict priority is a property of seed/time/UID/event kind. Biological
    // direction draws may use event_sequence, but that sequence must not bias
    // which same-time actor wins a reservation conflict.
    assert(active.priority == inactive.priority);
    const MoveProposal later_sequence = make_move_proposal(
        K_slot, symmetric_cells, symmetric_grid, symmetric_density,
        symmetric_config, 999, 1);
    assert(later_sequence.priority == active.priority);
    const MoveProposal later_time_bucket = make_move_proposal(
        K_slot, symmetric_cells, symmetric_grid, symmetric_density,
        symmetric_config, 1, 2);
    assert(later_time_bucket.priority != active.priority);

    // Production density filtering can reject all otherwise feasible initial
    // directions; the committed stay resets persistence.
    CellInit r_cell;
    r_cell.uid = 4321;
    r_cell.type = CellType::r;
    r_cell.last_direction = kStayDirection;
    r_cell.flags |= kMigrationActive;
    r_cell.anchor = {20, 20, 20};
    const Slot r_slot = symmetric_cells.create(r_cell);
    assert(symmetric_grid.place_single(r_cell.anchor, r_slot));
    symmetric_density.add(r_cell.anchor, r_cell.type, r_slot);
    symmetric_config.direction_density_threshold = 0.0;
    for (DirectionId direction = 1; direction <= 26; ++direction) {
        const Vec3i delta = direction_vector(direction);
        symmetric_density.add(r_cell.anchor + Vec3i{2 * delta.x, 2 * delta.y, 2 * delta.z},
                              CellType::K, static_cast<Slot>(100 + direction));
    }
    const MoveProposal blocked = make_move_proposal(
        r_slot, symmetric_cells, symmetric_grid, symmetric_density, symmetric_config, 0, 0);
    assert(blocked.direction == kStayDirection);
    assert(!commit_move(blocked, symmetric_cells, symmetric_grid, symmetric_density));
    assert(symmetric_cells.last_direction(r_slot) == kStayDirection);

    // A crowding exchange is an atomic two-cell transaction. It is available
    // only for two singleton stage-1 cells and swaps both the spatial grid and
    // density-index anchors without changing either stable UID/slot.
    Model3DConfig swap_config;
    swap_config.output_enabled = false;
    swap_config.direction_density_threshold = 1.0;
    DomainPolicy swap_domain(swap_config);
    SparseChunkGrid3D swap_grid(8, swap_domain);
    BlockDensityIndex3D swap_density(2);
    CellStore3D swap_cells;
    CellInit swap_actor;
    swap_actor.uid = 5001;
    swap_actor.type = CellType::r;
    swap_actor.stage = CellStage::small;
    swap_actor.anchor = {0, 0, 0};
    CellInit swap_partner = swap_actor;
    swap_partner.uid = 5002;
    swap_partner.type = CellType::K;
    swap_partner.anchor = {1, 0, 0};
    const Slot swap_actor_slot = swap_cells.create(swap_actor);
    const Slot swap_partner_slot = swap_cells.create(swap_partner);
    assert(swap_grid.place_single(swap_actor.anchor, swap_actor_slot));
    assert(swap_grid.place_single(swap_partner.anchor, swap_partner_slot));
    swap_density.add(swap_actor.anchor, swap_actor.type, swap_actor_slot);
    swap_density.add(swap_partner.anchor, swap_partner.type, swap_partner_slot);
    const auto swap_directions =
        feasible_crowding_swap_directions(
            swap_actor_slot, swap_cells, swap_grid, false);
    assert(swap_directions.size() == 1 && swap_directions.front() == 6);
    const MoveProposal exchange = make_crowding_swap_proposal(
        swap_actor_slot, 6, swap_cells, swap_grid, swap_config, 7);
    assert(exchange.swaps_anchors);
    assert(exchange.swap_partner == swap_partner_slot);
    assert(commit_crowding_swap(
        exchange, swap_cells, swap_grid, swap_density));
    assert(swap_cells.anchor(swap_actor_slot) == Vec3i(1, 0, 0));
    assert(swap_cells.anchor(swap_partner_slot) == Vec3i(0, 0, 0));
    assert(swap_grid.owner({1, 0, 0}) == swap_actor_slot);
    assert(swap_grid.owner({0, 0, 0}) == swap_partner_slot);
    assert(swap_cells.last_direction(swap_actor_slot) == 6);
    assert(swap_cells.last_direction(swap_partner_slot) == 2);
    assert(!commit_crowding_swap(
        exchange, swap_cells, swap_grid, swap_density));

    // Co-location groups are explicitly excluded: exchanging one member would
    // otherwise silently split a many-cell occupancy transaction.
    CellInit colocated = swap_partner;
    colocated.uid = 5003;
    colocated.anchor = {0, 0, 0};
    colocated.stage = CellStage::ultrasmall;
    const Slot colocated_slot = swap_cells.create(colocated);
    assert(swap_grid.add_colocated({0, 0, 0}, colocated_slot));
    swap_cells.set_stage(swap_partner_slot, CellStage::ultrasmall);
    assert(feasible_crowding_swap_directions(
               swap_actor_slot, swap_cells, swap_grid, false).empty());
}
