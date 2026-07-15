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
}
