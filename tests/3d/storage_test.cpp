#include <cassert>

#include "config/model_config.hpp"
#include "core/cell_store.hpp"
#include "space/chunk_grid.hpp"

int main() {
    using namespace atcg3d;
    static_assert(CellStore3D::logical_bytes_per_slot() <= 100);

    CellStore3D cells;
    cells.reserve(10);
    CellInit first;
    first.uid = 1;
    first.anchor = {-1, -1, -1};
    const Slot first_slot = cells.create(first);
    assert(cells.alive_count() == 1);
    cells.erase(first_slot);
    CellInit second;
    second.uid = 2;
    second.anchor = {4, 4, 4};
    const Slot second_slot = cells.create(second);
    assert(second_slot == first_slot);
    assert(cells.uid(second_slot) == 2);

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
    assert(grid.add_colocated({-1, -1, -1}, third_slot));
    assert(grid.occupants({-1, -1, -1}).size() == 2);
    assert(grid.remove({-1, -1, -1}, second_slot));
    assert(grid.owner({-1, -1, -1}) == third_slot);

    CellInit large;
    large.uid = 4;
    large.stage = CellStage::large;
    large.anchor = {7, 7, 7};
    const Slot large_slot = cells.create(large);
    assert(grid.place_large(large.anchor, large_slot));
    assert(grid.can_move_large(large.anchor, {1, 0, 0}));
    assert(!grid.can_place_large(large.anchor));
}
