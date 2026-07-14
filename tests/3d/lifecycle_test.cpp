#include <cassert>
#include <vector>

#include "config/model_config.hpp"
#include "core/cell_store.hpp"
#include "geometry/footprint.hpp"
#include "rules/initialization.hpp"
#include "rules/lifecycle.hpp"
#include "space/chunk_grid.hpp"
#include "space/density_index.hpp"

namespace {

atcg3d::Model3DConfig bounded_single_voxel_config() {
    atcg3d::Model3DConfig config;
    config.output_enabled = false;
    config.bounded_domain = true;
    config.domain_policy = "bounded";
    config.domain_min = {0, 0, 0};
    config.domain_max = {0, 0, 0};
    config.density_block_edge = 1;
    return config;
}

}  // namespace

int main() {
    using namespace atcg3d;

    // Large division first uses a complete daughter footprint on radius-2 shell.
    {
        Model3DConfig config;
        config.output_enabled = false;
        DomainPolicy domain(config);
        SparseChunkGrid3D grid(8, domain);
        BlockDensityIndex3D density(1);
        CellStore3D cells;
        CellInit mother_cell;
        mother_cell.uid = 1;
        mother_cell.stage = CellStage::large;
        mother_cell.anchor = {0, 0, 0};
        const Slot mother = cells.create(mother_cell);
        assert(grid.place_large(mother_cell.anchor, mother));
        density.add(mother_cell.anchor, mother_cell.type);
        CellUid next_uid = 2;
        std::vector<LineageEdge> lineage;
        const DivisionResult result = divide_cell(mother, 1.0, next_uid, cells, grid,
                                                  density, config, lineage);
        assert(result.changed && result.daughter != kEmptySlot);
        assert(cells.stage(result.daughter) == CellStage::large);
        assert(chebyshev_distance(cells.anchor(mother), cells.anchor(result.daughter)) == 2);
        assert(lineage.size() == 1 && lineage[0].parent_uid == 1);
    }

    // A bounded 4x4x4 region rejects every large shell candidate and exercises
    // the exact [-1,2]^3 shape-reduction fallback.
    {
        Model3DConfig config;
        config.output_enabled = false;
        config.bounded_domain = true;
        config.domain_policy = "bounded";
        config.domain_min = {-1, -1, -1};
        config.domain_max = {2, 2, 2};
        config.density_block_edge = 1;
        DomainPolicy domain(config);
        SparseChunkGrid3D grid(4, domain);
        BlockDensityIndex3D density(1);
        CellStore3D cells;
        CellInit initial;
        initial.uid = 10;
        initial.stage = CellStage::large;
        const Slot mother = cells.create(initial);
        assert(grid.place_large(initial.anchor, mother));
        density.add(initial.anchor, initial.type);
        CellUid next_uid = 11;
        std::vector<LineageEdge> lineage;
        const DivisionResult result = divide_cell(mother, 2.0, next_uid, cells, grid,
                                                  density, config, lineage);
        assert(result.changed && result.daughter != kEmptySlot);
        assert(cells.stage(mother) == CellStage::small);
        assert(cells.stage(result.daughter) == CellStage::small);
        assert(cells.anchor(mother) != cells.anchor(result.daughter));
    }

    // A blocked r cell follows the legacy non-survival rule.
    {
        Model3DConfig config = bounded_single_voxel_config();
        DomainPolicy domain(config);
        SparseChunkGrid3D grid(4, domain);
        BlockDensityIndex3D density(1);
        CellStore3D cells;
        CellInit initial;
        initial.uid = 20;
        initial.type = CellType::r;
        const Slot mother = cells.create(initial);
        assert(grid.place_single(initial.anchor, mother));
        density.add(initial.anchor, initial.type);
        CellUid next_uid = 21;
        std::vector<LineageEdge> lineage;
        const DivisionResult result = divide_cell(mother, 3.0, next_uid, cells, grid,
                                                  density, config, lineage);
        assert(result.mother_removed && !cells.valid(mother));
        assert(grid.occupants({0, 0, 0}).empty());
    }

    // A blocked K cell forms an explicit co-location group. Removing one member
    // leaves exactly one stage-1 cell; removing the last releases the voxel.
    {
        Model3DConfig config = bounded_single_voxel_config();
        DomainPolicy domain(config);
        SparseChunkGrid3D grid(4, domain);
        BlockDensityIndex3D density(1);
        CellStore3D cells;
        CellInit initial;
        initial.uid = 30;
        initial.type = CellType::K;
        const Slot mother = cells.create(initial);
        assert(grid.place_single(initial.anchor, mother));
        density.add(initial.anchor, initial.type);
        CellUid next_uid = 31;
        std::vector<LineageEdge> lineage;
        const DivisionResult result = divide_cell(mother, 4.0, next_uid, cells, grid,
                                                  density, config, lineage);
        assert(result.changed && grid.occupants({0, 0, 0}).size() == 2);
        assert(cells.stage(mother) == CellStage::ultrasmall);
        assert(cells.stage(result.daughter) == CellStage::ultrasmall);
        assert(remove_cell(result.daughter, cells, grid, density));
        assert(grid.occupants({0, 0, 0}).size() == 1);
        assert(cells.stage(mother) == CellStage::small);
        assert(remove_cell(mother, cells, grid, density));
        assert(grid.occupants({0, 0, 0}).empty());
    }

    // An isolated stage-1 cell enumerates the eight containing cubes and grows.
    {
        Model3DConfig config;
        config.output_enabled = false;
        DomainPolicy domain(config);
        SparseChunkGrid3D grid(8, domain);
        CellStore3D cells;
        CellInit initial;
        initial.uid = 40;
        const Slot slot = cells.create(initial);
        assert(grid.place_single(initial.anchor, slot));
        assert(try_stage_recovery(slot, cells, grid, config, 1));
        assert(cells.stage(slot) == CellStage::large);
        for (const Vec3i voxel : large_footprint(cells.anchor(slot))) {
            assert(grid.owner(voxel) == slot);
        }
    }

    // The initializer creates true 3D sphere/shell anchors without overlap.
    {
        Model3DConfig config;
        config.output_enabled = false;
        config.initial_r_cells = 8;
        config.initial_K_cells = 8;
        config.initial_radius = 7;
        config.initial_shell_thickness = 2;
        DomainPolicy domain(config);
        SparseChunkGrid3D grid(8, domain);
        BlockDensityIndex3D density(1);
        CellStore3D cells;
        const InitializationResult result = initialize_sphere_and_shell(cells, grid, density, config);
        assert(cells.alive_count() == 16 && result.next_uid == 17);
        bool saw_nonzero_z = false;
        for (const Slot slot : cells.alive_slots()) {
            saw_nonzero_z = saw_nonzero_z || cells.anchor(slot).z != 0;
            assert(grid.owner(cells.anchor(slot)) != kEmptySlot);
        }
        assert(saw_nonzero_z);
    }
}
