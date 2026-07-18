#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "config/model_config.hpp"
#include "core/cell_store.hpp"
#include "geometry/footprint.hpp"
#include "space/chunk_grid.hpp"
#include "vasculature/lesion_index.hpp"
#include "vasculature/surface_index.hpp"

namespace {

using namespace atcg3d;

void add_site_cell(LesionIndex3D& index, Vec3i site, double volume = 1.0) {
    index.add_cell_anchor(site, volume);
    index.add_occupied_site(site);
}

void fill_block(LesionIndex3D& index, Vec3i block, int edge) {
    for (int z = 0; z < edge; ++z) {
        for (int y = 0; y < edge; ++y) {
            for (int x = 0; x < edge; ++x) {
                const Vec3i site{block.x * edge + x, block.y * edge + y,
                                 block.z * edge + z};
                add_site_cell(index, site);
            }
        }
    }
}

LesionIndexConfig3D dense_block_config(int edge = 2) {
    LesionIndexConfig3D config;
    config.block_edge = edge;
    config.connectivity = LesionConnectivity3D::full_26;
    config.core_activation_occupied_fraction = 0.5;
    config.core_deactivation_occupied_fraction = 0.25;
    config.minimum_cells_per_core_block = 1;
    config.halo_blocks = 0;
    return config;
}

}  // namespace

int main() {
    using namespace atcg3d;

    // Strict construction validation.
    {
        bool rejected = false;
        try {
            LesionIndexConfig3D invalid;
            invalid.block_edge = 0;
            (void)LesionIndex3D(invalid);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        assert(rejected);
    }

    // Negative coordinates must use mathematical floor division.
    {
        LesionIndexConfig3D config = dense_block_config(8);
        LesionIndex3D index(config);
        assert(index.block_coordinate({-1, -1, -1}) == (Vec3i{-1, -1, -1}));
        assert(index.block_coordinate({-8, -8, -8}) == (Vec3i{-1, -1, -1}));
        assert(index.block_coordinate({-9, -9, -9}) == (Vec3i{-2, -2, -2}));
        assert(index.block_coordinate({0, 7, 8}) == (Vec3i{0, 0, 1}));
    }

    // Six-neighbor and 26-neighbor connectivity are distinct and configurable.
    {
        LesionIndexConfig3D six_config = dense_block_config(1);
        six_config.connectivity = LesionConnectivity3D::face_6;
        LesionIndex3D six(six_config);
        add_site_cell(six, {0, 0, 0});
        add_site_cell(six, {1, 1, 1});
        six.refresh_topology();
        assert(six.lesions().size() == 2);

        LesionIndexConfig3D full_config = six_config;
        full_config.connectivity = LesionConnectivity3D::full_26;
        LesionIndex3D full(full_config);
        add_site_cell(full, {1, 1, 1});  // reverse insertion must not matter
        add_site_cell(full, {0, 0, 0});
        full.refresh_topology();
        assert(full.lesions().size() == 1);
        assert(full.lesions().front().core_blocks.size() == 2);
    }

    // Sparse topology never scans the enormous empty bounding box between
    // distant observations, including near signed-coordinate limits.
    {
        LesionIndexConfig3D config = dense_block_config(1);
        LesionIndex3D index(config);
        const Vec3i low{std::numeric_limits<std::int32_t>::min(), 0, 0};
        const Vec3i high{std::numeric_limits<std::int32_t>::max(), 0, 0};
        add_site_cell(index, low);
        add_site_cell(index, high);
        index.refresh_topology();
        assert(index.observed_block_count() == 2);
        assert(index.core_block_count() == 2);
        assert(index.lesions().size() == 2);
    }

    // Core threshold hysteresis survives a full observation rebuild.
    {
        LesionIndex3D index(dense_block_config());
        assert(!index.refresh_needed());
        for (int x = 0; x < 4; ++x) add_site_cell(index, {x % 2, x / 2, 0});
        assert(index.topology_dirty());
        assert(index.statistics_dirty());
        index.refresh_topology();
        assert(!index.refresh_needed());
        assert(index.lesions().size() == 1);
        const LesionId stable_id = index.lesions().front().id;

        index.begin_rebuild();
        assert(index.refresh_needed());
        add_site_cell(index, {0, 0, 0});
        add_site_cell(index, {1, 0, 0});  // 2/8 == deactivation threshold
        index.refresh_topology();
        assert(index.lesions().size() == 1);
        assert(index.lesions().front().id == stable_id);

        index.begin_rebuild();
        add_site_cell(index, {0, 0, 0});  // below deactivation threshold
        index.refresh_topology();
        assert(index.lesions().empty());
    }

    // A sparse migration trail must not bridge two dense lesions. Halo blocks
    // may be attributed for surface lookup but never enter connectivity.
    {
        LesionIndexConfig3D config = dense_block_config();
        config.halo_blocks = 1;
        LesionIndex3D index(config);
        fill_block(index, {0, 0, 0}, 2);
        fill_block(index, {4, 0, 0}, 2);
        add_site_cell(index, {2, 0, 0});  // block 1: sparse, assigned to left
        add_site_cell(index, {4, 0, 0});  // block 2: sparse, no core in halo 1
        add_site_cell(index, {6, 0, 0});  // block 3: sparse, assigned to right
        const LesionTopologyDelta3D initial_delta = index.refresh_topology();
        assert(index.lesions().size() == 2);
        assert(initial_delta.created.size() == 2);
        assert(index.lesion_at_site({2, 0, 0}).has_value());
        assert(!index.lesion_at_site({4, 0, 0}).has_value());
        assert(index.lesion_at_site({6, 0, 0}).has_value());
        assert(index.lesion_at_site({2, 0, 0}) != index.lesion_at_site({6, 0, 0}));
    }

    // Largest-overlap continuation is deterministic through growth, merge and
    // split. On a split only the lexicographically first equal-size child can
    // retain the old ID; the other receives a fresh monotonic ID.
    {
        LesionIndex3D index(dense_block_config(1));
        add_site_cell(index, {0, 0, 0});
        add_site_cell(index, {4, 0, 0});
        index.refresh_topology();
        assert(index.lesions().size() == 2);
        const LesionId left_id = *index.lesion_at_site({0, 0, 0});
        const LesionId right_id = *index.lesion_at_site({4, 0, 0});
        assert(left_id < right_id);

        index.begin_rebuild();
        for (int x = 0; x <= 4; ++x) add_site_cell(index, {x, 0, 0});
        const LesionTopologyDelta3D merge_delta = index.refresh_topology();
        assert(index.lesions().size() == 1);
        assert(index.lesions().front().id == left_id);  // equal overlap merge tie
        assert(merge_delta.merges.size() == 1);
        assert(merge_delta.merges.front().result == left_id);
        assert((merge_delta.merges.front().predecessors ==
                std::vector<LesionId>{left_id, right_id}));
        assert(merge_delta.removed == std::vector<LesionId>{right_id});

        index.begin_rebuild();
        add_site_cell(index, {0, 0, 0});
        add_site_cell(index, {1, 0, 0});
        add_site_cell(index, {3, 0, 0});
        add_site_cell(index, {4, 0, 0});
        const LesionTopologyDelta3D split_delta = index.refresh_topology();
        assert(index.lesions().size() == 2);
        assert(*index.lesion_at_site({0, 0, 0}) == left_id);
        assert(*index.lesion_at_site({3, 0, 0}) != left_id);
        assert(*index.lesion_at_site({3, 0, 0}) > right_id);
        assert(split_delta.splits.size() == 1);
        assert(split_delta.splits.front().predecessor == left_id);
        assert(split_delta.splits.front().retained_child == left_id);
        assert(split_delta.splits.front().children.size() == 2);
        assert(split_delta.created.size() == 1);
    }

    // Checkpoint identity is sorted, strict, and restores stable IDs without
    // serializing derived halo attribution.
    {
        LesionIndex3D original(dense_block_config(1));
        add_site_cell(original, {-3, 0, 0});
        add_site_cell(original, {7, 0, 0});
        original.refresh_topology();
        const std::vector<LesionCoreIdentity3D> identity =
            original.snapshot_core_identity();
        assert(identity.size() == 2);
        assert(identity[0].block == (Vec3i{-3, 0, 0}));
        assert(identity[1].block == (Vec3i{7, 0, 0}));

        LesionIndex3D restored(dense_block_config(1));
        add_site_cell(restored, {7, 0, 0});
        add_site_cell(restored, {-3, 0, 0});
        restored.refresh_topology();
        restored.restore_core_identity(identity, 100);
        assert(restored.snapshot_core_identity() == identity);
        assert(restored.next_lesion_id() == 100);

        std::vector<LesionCoreIdentity3D> duplicate = identity;
        duplicate.push_back(identity.front());
        bool rejected_duplicate = false;
        try {
            restored.restore_core_identity(duplicate, 100);
        } catch (const std::invalid_argument&) {
            rejected_duplicate = true;
        }
        assert(rejected_duplicate);

        std::vector<LesionCoreIdentity3D> zero_id = identity;
        zero_id.front().lesion_id = kNoLesionId;
        bool rejected_zero = false;
        try {
            restored.restore_core_identity(zero_id, 100);
        } catch (const std::invalid_argument&) {
            rejected_zero = true;
        }
        assert(rejected_zero);
    }

    // Statistics count biological cells once while occupancy can reflect a
    // multi-voxel footprint. Face ownership follows the inside surface site.
    {
        LesionIndexConfig3D config = dense_block_config(2);
        config.core_activation_occupied_fraction = 0.25;
        config.core_deactivation_occupied_fraction = 0.25;
        config.halo_blocks = 1;
        LesionIndex3D index(config);
        index.add_cell_anchor({-1, 0, 0}, 8.0);  // one stage-0-like cell
        index.add_occupied_site({-1, 0, 0});
        index.add_occupied_site({-2, 0, 0});
        index.add_cell_anchor({0, 0, 0}, 1.0);
        index.add_occupied_site({0, 0, 0});
        index.add_occupied_site({1, 0, 0});
        index.refresh_topology();
        assert(index.lesions().size() == 1);
        const LesionSummary3D& lesion = index.lesions().front();
        assert(lesion.cell_count == 2);
        assert(lesion.occupied_voxel_count == 4);
        assert(std::abs(lesion.biological_volume - 9.0) < 1e-12);
        assert(std::abs(lesion.centroid.x + 0.5) < 1e-12);
        assert(lesion.minimum_site.x == -2);
        assert(lesion.maximum_site.x == 1);
        const ExposedFace3D face{{-1, 0, 0}, {-1, 0, 0}};
        assert(index.lesion_for_face(face) == lesion.id);
        assert(index.lesion_owns_face(lesion.id, face));
        assert(index.find_lesion(lesion.id) != nullptr);
        assert(index.allocated_bytes() > 0);
    }

    // Incremental add/remove API updates observations without rebuilding a
    // dense domain and rejects count underflow.
    {
        LesionIndex3D index(dense_block_config(1));
        add_site_cell(index, {-1, 0, 0});
        index.refresh_topology();
        assert(index.lesions().size() == 1);
        index.remove_cell_anchor({-1, 0, 0});
        index.remove_occupied_site({-1, 0, 0});
        index.refresh_topology();
        assert(index.lesions().empty());
        assert(index.observed_block_count() == 0);
        bool rejected = false;
        try {
            index.remove_occupied_site({-1, 0, 0});
        } catch (const std::logic_error&) {
            rejected = true;
        }
        assert(rejected);
    }

    // Full CellStore/Grid rebuild counts a large cell once while obtaining its
    // eight unique occupied voxels directly from the sparse space index.
    {
        Model3DConfig domain_config;
        domain_config.chunk_edge = 4;
        CellStore3D cells;
        SparseChunkGrid3D grid(domain_config.chunk_edge, DomainPolicy(domain_config));

        CellInit large;
        large.uid = 1;
        large.anchor = {0, 0, 0};
        large.stage = CellStage::large;
        const Slot large_slot = cells.create(large);
        assert(grid.place_large(large.anchor, large_slot));

        CellInit isolated;
        isolated.uid = 2;
        isolated.anchor = {20, 0, 0};
        isolated.stage = CellStage::small;
        const Slot isolated_slot = cells.create(isolated);
        assert(grid.place_single(isolated.anchor, isolated_slot));

        LesionIndexConfig3D lesion_config = dense_block_config(2);
        lesion_config.halo_blocks = 0;
        LesionIndex3D index(lesion_config);
        index.rebuild_from(cells, grid);
        assert(index.lesions().size() == 1);
        assert(index.lesions().front().cell_count == 1);
        assert(index.lesions().front().occupied_voxel_count == 8);
        assert(index.lesions().front().biological_volume == 8.0);
        assert(index.lesion_for_site({0, 0, 0}) == index.lesions().front().id);
        assert(!index.lesion_for_site(isolated.anchor).has_value());

        // A colocated addition changes biological cell count without changing
        // voxel occupancy. Marking that site refreshes only its coarse block
        // and deduplicates primary/overflow occupants by stable slot.
        CellInit colocated = isolated;
        colocated.uid = 3;
        colocated.stage = CellStage::ultrasmall;
        cells.set_stage(isolated_slot, CellStage::ultrasmall);
        const Slot colocated_slot = cells.create(colocated);
        assert(grid.add_colocated(isolated.anchor, colocated_slot));
        const std::array<Vec3i, 1> changed{{isolated.anchor}};
        assert(index.mark_dirty_sites(changed) == 1);
        assert(index.observations_dirty());
        bool rejected_stale_topology = false;
        try {
            index.refresh_topology();
        } catch (const std::logic_error&) {
            rejected_stale_topology = true;
        }
        assert(rejected_stale_topology);
        assert(index.refresh_dirty_blocks_from(cells, grid) == 1);
        assert(!index.observations_dirty());
        const auto isolated_block = index.block_stats(index.block_coordinate(isolated.anchor));
        assert(isolated_block.has_value());
        assert(isolated_block->cell_count == 2);
        assert(isolated_block->occupied_voxel_count == 1);
        assert(!index.topology_dirty());  // still below the occupied-fraction threshold
        index.refresh_topology();
        assert(index.lesions().size() == 1);
    }

    // A checkpoint taken between lesion refreshes must preserve the old
    // observations/topology while the CellStore/Grid already contain newer
    // occupancy. Resume replays the same later refresh and obtains the same
    // stable IDs and summaries as an uninterrupted index.
    {
        Model3DConfig domain_config;
        domain_config.chunk_edge = 4;
        CellStore3D cells;
        SparseChunkGrid3D grid(domain_config.chunk_edge, DomainPolicy(domain_config));

        CellInit initial;
        initial.uid = 10;
        initial.anchor = {-2, 0, 0};
        initial.stage = CellStage::large;
        const Slot initial_slot = cells.create(initial);
        assert(grid.place_large(initial.anchor, initial_slot));

        LesionIndexConfig3D config = dense_block_config(2);
        config.halo_blocks = 0;
        LesionIndex3D uninterrupted(config);
        uninterrupted.rebuild_from(cells, grid);
        const std::vector<LesionCoreIdentity3D> checkpoint_identity =
            uninterrupted.snapshot_core_identity();
        const LesionId checkpoint_next_id = uninterrupted.next_lesion_id();
        assert(checkpoint_identity.size() == 1);
        const LesionSummary3D checkpoint_summary = uninterrupted.lesions().front();

        const auto old_footprint = large_footprint(initial.anchor);
        grid.remove_large(initial.anchor, initial_slot);
        cells.erase(initial_slot);
        CellInit moved = initial;
        moved.uid = 11;
        moved.anchor = {4, 0, 0};
        const Slot moved_slot = cells.create(moved);
        assert(grid.place_large(moved.anchor, moved_slot));
        const auto new_footprint = large_footprint(moved.anchor);
        std::vector<Vec3i> changed(old_footprint.begin(), old_footprint.end());
        changed.insert(changed.end(), new_footprint.begin(), new_footprint.end());
        uninterrupted.mark_dirty_sites(changed);

        const std::vector<LesionDirtyBlockState3D> checkpoint_dirty =
            uninterrupted.snapshot_dirty_block_state();
        assert(checkpoint_dirty.size() == 2);
        assert(checkpoint_dirty[0].block == (Vec3i{-1, 0, 0}));
        assert(checkpoint_dirty[0].exists);
        assert(checkpoint_dirty[0].cell_count == 1);
        assert(checkpoint_dirty[0].occupied_voxel_count == 8);
        assert(checkpoint_dirty[0].biological_volume == 8.0);
        assert(checkpoint_dirty[0].cell_coordinate_sum_x == -2);
        assert(checkpoint_dirty[0].occupied_coordinate_sum_x == -12);
        assert(checkpoint_dirty[0].occupied_coordinate_sum_y == 4);
        assert(checkpoint_dirty[0].occupied_coordinate_sum_z == 4);
        assert(checkpoint_dirty[1].block == (Vec3i{2, 0, 0}));
        assert(!checkpoint_dirty[1].exists);

        LesionIndex3D resumed(config);
        resumed.rebuild_from(cells, grid);  // observes the moved cell too early
        assert(resumed.lesion_for_anchor(moved.anchor).has_value());
        resumed.restore_checkpoint_state(checkpoint_identity, checkpoint_next_id,
                                         checkpoint_dirty);
        assert(resumed.observations_dirty());
        assert(resumed.dirty_block_count() == 2);
        assert(resumed.snapshot_dirty_block_state() == checkpoint_dirty);
        assert(resumed.snapshot_core_identity() == checkpoint_identity);
        assert(resumed.next_lesion_id() == checkpoint_next_id);
        assert(resumed.lesions().size() == 1);
        assert(resumed.lesions().front().id == checkpoint_summary.id);
        assert(resumed.lesions().front().cell_count == checkpoint_summary.cell_count);
        assert(resumed.lesions().front().occupied_voxel_count ==
               checkpoint_summary.occupied_voxel_count);
        assert(resumed.lesions().front().biological_volume ==
               checkpoint_summary.biological_volume);
        assert(resumed.lesions().front().centroid == checkpoint_summary.centroid);
        assert(resumed.lesion_for_anchor(initial.anchor) == checkpoint_summary.id);
        assert(!resumed.lesion_for_anchor(moved.anchor).has_value());

        assert(uninterrupted.refresh_dirty_blocks_from(cells, grid) == 2);
        assert(resumed.refresh_dirty_blocks_from(cells, grid) == 2);
        const LesionTopologyDelta3D uninterrupted_delta =
            uninterrupted.refresh_topology();
        const LesionTopologyDelta3D resumed_delta = resumed.refresh_topology();
        assert(uninterrupted.snapshot_core_identity() ==
               resumed.snapshot_core_identity());
        assert(uninterrupted.next_lesion_id() == resumed.next_lesion_id());
        assert(uninterrupted_delta.created == resumed_delta.created);
        assert(uninterrupted_delta.removed == resumed_delta.removed);
        assert(uninterrupted.lesions().size() == resumed.lesions().size());
        assert(uninterrupted.lesions().front().id == resumed.lesions().front().id);
        assert(uninterrupted.lesions().front().centroid ==
               resumed.lesions().front().centroid);

        const auto expect_invalid_dirty = [&](std::vector<LesionDirtyBlockState3D> bad) {
            LesionIndex3D candidate(config);
            candidate.rebuild_from(cells, grid);
            bool rejected = false;
            try {
                candidate.restore_checkpoint_state(
                    checkpoint_identity, checkpoint_next_id, bad);
            } catch (const std::invalid_argument&) {
                rejected = true;
            }
            assert(rejected);
        };

        std::vector<LesionDirtyBlockState3D> duplicate = checkpoint_dirty;
        duplicate.push_back(checkpoint_dirty.front());
        expect_invalid_dirty(duplicate);

        std::vector<LesionDirtyBlockState3D> nan_volume = checkpoint_dirty;
        nan_volume[0].biological_volume =
            std::numeric_limits<double>::quiet_NaN();
        expect_invalid_dirty(nan_volume);

        std::vector<LesionDirtyBlockState3D> negative_volume = checkpoint_dirty;
        negative_volume[0].biological_volume = -1.0;
        expect_invalid_dirty(negative_volume);

        std::vector<LesionDirtyBlockState3D> excessive_occupancy = checkpoint_dirty;
        excessive_occupancy[0].occupied_voxel_count = 9;
        expect_invalid_dirty(excessive_occupancy);

        std::vector<LesionDirtyBlockState3D> absent_payload = checkpoint_dirty;
        absent_payload[1].cell_count = 1;
        expect_invalid_dirty(absent_payload);

        std::vector<LesionDirtyBlockState3D> empty_present = checkpoint_dirty;
        empty_present[0] = {.block = {-1, 0, 0}, .exists = true};
        expect_invalid_dirty(empty_present);

        std::vector<LesionDirtyBlockState3D> invalid_sum = checkpoint_dirty;
        invalid_sum[0].cell_coordinate_sum_x = 0;
        expect_invalid_dirty(invalid_sum);

        std::vector<LesionDirtyBlockState3D> excessive_cells = checkpoint_dirty;
        excessive_cells[0].cell_count =
            static_cast<std::uint64_t>(std::numeric_limits<Slot>::max()) + 1U;
        expect_invalid_dirty(excessive_cells);
    }

    return 0;
}
