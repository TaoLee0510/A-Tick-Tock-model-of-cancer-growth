#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "config/model_config.hpp"
#include "core/cell_store.hpp"
#include "geometry/footprint.hpp"
#include "rules/initialization.hpp"
#include "rules/lifecycle.hpp"
#include "space/chunk_grid.hpp"
#include "space/density_index.hpp"
#include "vasculature/vessel_grid.hpp"

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

    {
        DivisionTimingConfig timing;
        timing.base_cycle_hours = 30.0;
        timing.minimum_fraction = 1.0;
        timing.stochastic_tail_fraction = 0.0;
        timing.stochastic_time_quantum_hours = 0.5;
        assert(sample_division_delay(2.0, timing, 1, 1, 1) == 15.5);
        assert(sample_death_delay(48.0, 7, 99, 3) >= 1.0);
        assert(sample_death_delay(48.0, 7, 99, 3) ==
               sample_death_delay(48.0, 7, 99, 3));
    }

    // The conversion decision is thresholded and deterministic without
    // consuming mutable RNG state. The mapped density estimator remains
    // stage-aware: the same anchors occupy 2/8 of a 2^3 small-cell window and
    // 2/(8/8) of a large-cell-capacity window.
    {
        RToKConversionConfig conversion;
        conversion.enabled = true;
        conversion.density_window_edge = 2;
        conversion.query_block_edge = 1;
        conversion.density_threshold = 0.25;
        conversion.probability_per_division = 1.0;
        BlockDensityIndex3D density(1);
        density.add({0, 0, 0}, CellType::r, 0);
        density.add({1, 0, 0}, CellType::K, 1);
        assert(r_to_K_division_density(
                   density, {0, 0, 0}, CellStage::small, conversion) == 0.25);
        assert(r_to_K_division_density(
                   density, {0, 0, 0}, CellStage::large, conversion) == 2.0);
        assert(!should_convert_r_daughter(
            CellType::r, 0.249, conversion, 9, 10, 11));
        assert(should_convert_r_daughter(
            CellType::r, 0.25, conversion, 9, 10, 11));
        assert(!should_convert_r_daughter(
            CellType::K, 1.0, conversion, 9, 10, 11));
        conversion.probability_per_division = 0.37;
        assert(should_convert_r_daughter(
                   CellType::r, 0.25, conversion, 9, 10, 11) ==
               should_convert_r_daughter(
                   CellType::r, 0.25, conversion, 9, 10, 11));
        conversion.probability_per_division = 0.0;
        assert(!should_convert_r_daughter(
            CellType::r, 1.0, conversion, 9, 10, 11));
    }

    // Density refresh integrates the old rate over elapsed time and preserves
    // the one sampled cycle instead of consuming RNG and drawing a new delay.
    {
        Model3DConfig config;
        config.output_enabled = false;
        config.migration_activation_enabled = false;
        BlockDensityIndex3D density(1);
        CellStore3D cells;
        CellInit initial;
        initial.uid = 1000;
        initial.density_growth_rate = 2.0F;
        initial.division_work_remaining = 10.0F;
        initial.next_division_time = 5.0;
        const Slot slot = cells.create(initial);
        density.add(initial.anchor, initial.type, slot);
        const std::uint64_t sequence = cells.event_sequence(slot);
        const GrowthRefreshResult first =
            refresh_growth_state(slot, 1.0, cells, density, config);
        assert(cells.division_work_remaining(slot) == 8.0F);
        assert(cells.event_sequence(slot) == sequence);
        const double division_time_after_first = cells.next_division_time(slot);
        const float rate_after_first = cells.density_growth_rate(slot);
        const float work_after_first = cells.division_work_remaining(slot);
        (void)refresh_growth_state(slot, 1.0, cells, density, config);
        assert(cells.division_work_remaining(slot) == work_after_first);
        assert(cells.event_sequence(slot) == sequence);
        assert(cells.next_division_time(slot) == division_time_after_first);
        (void)refresh_growth_state(slot, 2.0, cells, density, config);
        const float expected = std::max(0.0F, work_after_first - rate_after_first);
        assert(std::abs(cells.division_work_remaining(slot) - expected) < 1.0e-5F);
        assert(cells.event_sequence(slot) == sequence);
        assert(first.division_time_changed);

        initialize_division_cycle(slot, 2.0, cells, config);
        assert(cells.division_work_remaining(slot) > 0.0F);
        assert(cells.event_sequence(slot) == sequence + 1);
    }

    // Density triggers a finite active interval rather than gating migration.
    // Low-density cells keep an ordinary rate, falling density does not end an
    // interval early, expiry resamples the normal r rate, and later density can
    // activate the same cell again.
    {
        Model3DConfig config;
        config.output_enabled = false;
        config.migration_activation_enabled = true;
        config.migration_activation_window_edge = 2;
        config.migration_activation_block_edge = 1;
        config.migration_activation_threshold = 0.20;
        BlockDensityIndex3D density(1);
        CellStore3D cells;
        CellInit first;
        first.uid = 1100;
        first.migration_rate = 2.0F;
        first.normal_migration_rate = 0.25F;
        first.density_growth_rate = 1.0F;
        first.division_work_remaining = 10.0F;
        const Slot first_slot = cells.create(first);
        density.add({0, 0, 0}, first.type, first_slot);
        assert(!refresh_migration_activation_state(
            first_slot, 0.0, cells, density, config));
        assert(migration_allowed_for_cell(first_slot, cells, config));
        assert(effective_migration_rate(first_slot, cells, config) == 0.25);

        CellInit second;
        second.uid = 1101;
        second.anchor = {1, 0, 0};
        const Slot second_slot = cells.create(second);
        density.add(second.anchor, second.type, second_slot);
        assert(refresh_migration_activation_state(
            first_slot, 0.0, cells, density, config));
        const double first_end = cells.migration_activation_end_time(first_slot);
        assert(first_end > 0.0 && first_end < 10.0);
        assert(effective_migration_rate(first_slot, cells, config) == 2.0);
        density.remove(second.anchor, second.type, second_slot);
        assert(!refresh_migration_activation_state(
            first_slot, first_end * 0.5, cells, density, config));
        assert((cells.flags(first_slot) & kMigrationActive) != 0);
        assert(expire_migration_activation_state(
            first_slot, first_end, cells, config));
        assert((cells.flags(first_slot) & kMigrationActive) == 0);
        assert(cells.migration_activation_end_time(first_slot) == 0.0);
        assert(effective_migration_rate(first_slot, cells, config) > 0.0);

        density.add(second.anchor, second.type, second_slot);
        assert(refresh_migration_activation_state(
            first_slot, first_end, cells, density, config));
        assert(cells.migration_activation_end_time(first_slot) > first_end);

        // If the remaining cycle is shorter than one float ULP at the current
        // clock, no strictly-future in-cycle end event is representable. The
        // cell safely remains inactive and does not consume an RNG sequence.
        CellInit tiny_interval;
        tiny_interval.uid = 1102;
        tiny_interval.anchor = {0, 1, 0};
        tiny_interval.density_growth_rate = 1.0F;
        tiny_interval.division_work_remaining = 1.0e-8F;
        const double tiny_now = 1024.0 + 1.0e-8;
        tiny_interval.last_update_time = tiny_now;
        const Slot tiny_slot = cells.create(tiny_interval);
        density.add(tiny_interval.anchor, tiny_interval.type, tiny_slot);
        config.migration_activation_threshold = 0.0;
        const std::uint64_t tiny_sequence = cells.event_sequence(tiny_slot);
        assert(!refresh_migration_activation_state(
            tiny_slot, tiny_now, cells, density, config));
        assert((cells.flags(tiny_slot) & kMigrationActive) == 0);
        assert(cells.migration_activation_end_time(tiny_slot) == 0.0);
        assert(cells.event_sequence(tiny_slot) == tiny_sequence);

        config.migration_activation_enabled = false;
        assert(migration_allowed_for_cell(first_slot, cells, config));
    }

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
        mother_cell.flags |= kMigrationActive;
        mother_cell.migration_rate = 2.0F;
        mother_cell.normal_migration_rate = 0.2F;
        mother_cell.migration_activation_end_time = 100.0;
        const Slot mother = cells.create(mother_cell);
        assert(grid.place_large(mother_cell.anchor, mother));
        density.add(mother_cell.anchor, mother_cell.type, mother);
        CellUid next_uid = 2;
        std::vector<LineageEdge> lineage;
        const DivisionResult result = divide_cell(mother, 1.0, next_uid, cells, grid,
                                                  density, config, lineage);
        assert(result.changed && result.daughter != kEmptySlot);
        assert(cells.stage(result.daughter) == CellStage::large);
        assert(chebyshev_distance(cells.anchor(mother), cells.anchor(result.daughter)) == 2);
        assert(lineage.size() == 1 && lineage[0].parent_uid == 1);
        assert((cells.flags(mother) & kMigrationActive) == 0);
        assert((cells.flags(result.daughter) & kMigrationActive) == 0);
        assert(cells.migration_activation_end_time(mother) == 0.0);
        assert(cells.migration_activation_end_time(result.daughter) == 0.0);
        assert(cells.migration_rate(mother) == config.initial_r_migration_rate);
        assert(cells.migration_rate(result.daughter) ==
               config.initial_r_migration_rate);
        assert(cells.normal_migration_rate(mother) > 0.0F);
        assert(cells.normal_migration_rate(result.daughter) > 0.0F);
    }

    // A successful dense r division can create a K daughter. Converted rates
    // are freshly sampled from the configured K model, then capped; the r
    // mother's inherited rate is independently capped at the r limit.
    {
        Model3DConfig config;
        config.output_enabled = false;
        config.migration_activation_enabled = false;
        config.density_block_edge = 1;
        config.r_to_K_conversion.enabled = true;
        config.r_to_K_conversion.density_window_edge = 2;
        config.r_to_K_conversion.query_block_edge = 1;
        config.r_to_K_conversion.density_threshold = 0.25;
        config.r_to_K_conversion.probability_per_division = 1.0;
        config.division_timing.inherited_growth_multiplier_min = 2.0;
        config.division_timing.inherited_growth_multiplier_max = 2.0;
        config.division_timing.r_max_inherent_growth_rate = 1.10;
        config.division_timing.K_max_inherent_growth_rate = 0.90;
        config.initial_growth_rate_model = "fixed";
        config.initial_K_growth_rate = 2.5;
        config.initial_migration_rate_model = "fixed";
        config.initial_K_migration_rate = 0.125;
        DomainPolicy domain(config);
        SparseChunkGrid3D grid(8, domain);
        BlockDensityIndex3D density(1);
        CellStore3D cells;
        CellInit mother_cell;
        mother_cell.uid = 200;
        mother_cell.type = CellType::r;
        mother_cell.inherent_growth_rate = 1.0F;
        const Slot mother = cells.create(mother_cell);
        assert(grid.place_single(mother_cell.anchor, mother));
        density.add(mother_cell.anchor, mother_cell.type, mother);
        CellInit neighbor;
        neighbor.uid = 201;
        neighbor.type = CellType::K;
        neighbor.anchor = {0, 1, 0};
        const Slot neighbor_slot = cells.create(neighbor);
        assert(grid.place_single(neighbor.anchor, neighbor_slot));
        density.add(neighbor.anchor, neighbor.type, neighbor_slot);

        CellUid next_uid = 202;
        std::vector<LineageEdge> lineage;
        const DivisionResult result = divide_cell(
            mother, 1.0, next_uid, cells, grid, density, config, lineage);
        assert(result.changed && result.daughter != kEmptySlot);
        assert(cells.type(mother) == CellType::r);
        assert(cells.type(result.daughter) == CellType::K);
        assert(std::abs(cells.inherent_growth_rate(mother) - 1.10F) < 1.0e-6F);
        assert(std::abs(cells.inherent_growth_rate(result.daughter) - 0.90F) < 1.0e-6F);
        assert(std::abs(cells.migration_rate(result.daughter) - 0.125F) < 1.0e-6F);
        assert(std::abs(cells.normal_migration_rate(result.daughter) - 0.125F) < 1.0e-6F);
        assert(lineage.size() == 1 && lineage.front().type == CellType::K);
        const DensityCounts3D counts = density.estimate_box({-2, -2, -2}, {2, 2, 2});
        assert(counts.r == 1 && counts.K == 2);
    }

    // Inherited K rates are capped for both mother and daughter when no type
    // conversion applies.
    {
        Model3DConfig config;
        config.output_enabled = false;
        config.migration_activation_enabled = false;
        config.division_timing.inherited_growth_multiplier_min = 2.0;
        config.division_timing.inherited_growth_multiplier_max = 2.0;
        config.division_timing.K_max_inherent_growth_rate = 0.75;
        DomainPolicy domain(config);
        SparseChunkGrid3D grid(8, domain);
        BlockDensityIndex3D density(1);
        CellStore3D cells;
        CellInit initial;
        initial.uid = 210;
        initial.type = CellType::K;
        initial.inherent_growth_rate = 1.0F;
        const Slot mother = cells.create(initial);
        assert(grid.place_single(initial.anchor, mother));
        density.add(initial.anchor, initial.type, mother);
        CellUid next_uid = 211;
        std::vector<LineageEdge> lineage;
        const DivisionResult result = divide_cell(
            mother, 1.0, next_uid, cells, grid, density, config, lineage);
        assert(result.changed && result.daughter != kEmptySlot);
        assert(cells.inherent_growth_rate(mother) == 0.75F);
        assert(cells.inherent_growth_rate(result.daughter) == 0.75F);
    }

    // A failed r division does not commit its conversion draw, inherited-rate
    // cap, or event sequence. This keeps an eventual retry bit-reproducible.
    {
        Model3DConfig config;
        config.output_enabled = false;
        config.bounded_domain = true;
        config.domain_policy = "bounded";
        config.domain_min = {0, 0, 0};
        config.domain_max = {1, 1, 1};
        config.density_block_edge = 1;
        config.allow_shape_reduction = false;
        config.r_to_K_conversion.enabled = true;
        config.r_to_K_conversion.density_window_edge = 2;
        config.r_to_K_conversion.query_block_edge = 1;
        config.r_to_K_conversion.density_threshold = 0.0;
        config.r_to_K_conversion.probability_per_division = 1.0;
        config.division_timing.inherited_growth_multiplier_min = 2.0;
        config.division_timing.inherited_growth_multiplier_max = 2.0;
        config.division_timing.r_max_inherent_growth_rate = 0.5;
        DomainPolicy domain(config);
        SparseChunkGrid3D grid(4, domain);
        BlockDensityIndex3D density(1);
        CellStore3D cells;
        CellInit initial;
        initial.uid = 220;
        initial.type = CellType::r;
        initial.stage = CellStage::large;
        initial.inherent_growth_rate = 1.25F;
        initial.migration_rate = 0.4F;
        initial.event_sequence = 7;
        const Slot mother = cells.create(initial);
        assert(grid.place_large(initial.anchor, mother));
        density.add(initial.anchor, initial.type, mother);
        CellUid next_uid = 221;
        std::vector<LineageEdge> lineage;
        const DivisionResult result = divide_cell(
            mother, 1.0, next_uid, cells, grid, density, config, lineage);
        assert(!result.changed && result.daughter == kEmptySlot);
        assert(cells.type(mother) == CellType::r);
        assert(cells.inherent_growth_rate(mother) == 1.25F);
        assert(std::abs(cells.migration_rate(mother) - 0.4F) < 1.0e-6F);
        assert(cells.event_sequence(mother) == 7);
        assert(next_uid == 221 && lineage.empty());
        for (const Vec3i voxel : large_footprint(initial.anchor)) {
            assert(grid.owner(voxel) == mother);
        }
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
        SparseVesselGrid3D vessels(4, domain);
        SparseChunkGrid3D grid(4, domain);
        grid.attach_vessel_grid(&vessels);
        BlockDensityIndex3D density(1);
        CellStore3D cells;
        CellInit initial;
        initial.uid = 10;
        initial.stage = CellStage::large;
        const Slot mother = cells.create(initial);
        assert(grid.place_large(initial.anchor, mother));
        density.add(initial.anchor, initial.type, mother);
        // Reserve every shape-reduction site outside the mother's current
        // footprint. The fallback must choose only vessel-free former mother
        // voxels and commit both single-voxel placements atomically.
        for (const Vec3i site : shape_reduction_sites(initial.anchor)) {
            if (grid.owner(site) == kEmptySlot) {
                assert(vessels.add(site, VesselBranchRole::outward).placed);
            }
        }
        CellUid next_uid = 11;
        std::vector<LineageEdge> lineage;
        const DivisionResult result = divide_cell(mother, 2.0, next_uid, cells, grid,
                                                  density, config, lineage);
        assert(result.changed && result.daughter != kEmptySlot);
        assert(cells.stage(mother) == CellStage::small);
        assert(cells.stage(result.daughter) == CellStage::small);
        assert(cells.anchor(mother) != cells.anchor(result.daughter));
        assert(!grid.blocked_by_vessel(cells.anchor(mother)));
        assert(!grid.blocked_by_vessel(cells.anchor(result.daughter)));
        std::size_t occupied_old_footprint = 0;
        for (const Vec3i site : large_footprint(initial.anchor)) {
            if (grid.owner(site) != kEmptySlot) ++occupied_old_footprint;
        }
        assert(occupied_old_footprint == 2);
        assert(density.estimate_box({-1, -1, -1}, {2, 2, 2}).total() == 2);
    }

    // Thin-layer shape reduction keeps biological anchors on z=0 even though
    // z=1 remains addressable as the upper footprint half of a large cell.
    {
        Model3DConfig config;
        config.output_enabled = false;
        config.thin_layer = true;
        config.bounded_domain = true;
        config.domain_policy = "bounded";
        config.domain_min = {0, 0, 0};
        config.domain_max = {1, 1, 1};
        config.density_block_edge = 1;
        DomainPolicy domain(config);
        SparseChunkGrid3D grid(4, domain);
        BlockDensityIndex3D density(1);
        CellStore3D cells;
        CellInit initial;
        initial.uid = 15;
        initial.stage = CellStage::large;
        const Slot mother = cells.create(initial);
        assert(grid.place_large(initial.anchor, mother));
        density.add(initial.anchor, initial.type, mother);
        CellUid next_uid = 16;
        std::vector<LineageEdge> lineage;
        const DivisionResult result = divide_cell(
            mother, 2.0, next_uid, cells, grid, density, config, lineage);
        assert(result.changed && result.daughter != kEmptySlot);
        assert(cells.stage(mother) == CellStage::small);
        assert(cells.stage(result.daughter) == CellStage::small);
        assert(cells.anchor(mother).z == 0);
        assert(cells.anchor(result.daughter).z == 0);
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
        density.add(initial.anchor, initial.type, mother);
        CellUid next_uid = 21;
        std::vector<LineageEdge> lineage;
        const DivisionResult result = divide_cell(mother, 3.0, next_uid, cells, grid,
                                                  density, config, lineage);
        assert(result.mother_removed && !cells.valid(mother));
        assert(grid.occupants({0, 0, 0}).empty());
    }

    // A failed biological division is transactional for inherited state. A K
    // cell with ultrasmall fallback disabled keeps its original inherent rate,
    // UID stream, lineage, stage, and occupancy when no daughter can be placed.
    {
        Model3DConfig config = bounded_single_voxel_config();
        config.ultrasmall_enabled = false;
        config.division_timing.inherited_growth_multiplier_min = 0.5;
        config.division_timing.inherited_growth_multiplier_max = 0.5;
        DomainPolicy domain(config);
        SparseChunkGrid3D grid(4, domain);
        BlockDensityIndex3D density(1);
        CellStore3D cells;
        CellInit initial;
        initial.uid = 25;
        initial.type = CellType::K;
        initial.inherent_growth_rate = 1.25F;
        const Slot mother = cells.create(initial);
        assert(grid.place_single(initial.anchor, mother));
        density.add(initial.anchor, initial.type, mother);
        CellUid next_uid = 26;
        std::vector<LineageEdge> lineage;
        const std::uint64_t sequence_before = cells.event_sequence(mother);
        const DivisionResult result = divide_cell(mother, 3.5, next_uid, cells, grid,
                                                  density, config, lineage);
        assert(!result.changed && !result.mother_removed && result.daughter == kEmptySlot);
        assert(cells.valid(mother));
        assert(cells.inherent_growth_rate(mother) == 1.25F);
        assert(cells.stage(mother) == CellStage::small);
        assert(grid.owner(initial.anchor) == mother);
        assert(cells.event_sequence(mother) == sequence_before);
        assert(next_uid == 26 && lineage.empty());
    }

    // A blocked K cell forms an explicit co-location group. Removing one member
    // leaves exactly one stage-1 cell; removing the last releases the voxel.
    {
        Model3DConfig config = bounded_single_voxel_config();
        config.ultrasmall_enabled = true;
        DomainPolicy domain(config);
        SparseChunkGrid3D grid(4, domain);
        BlockDensityIndex3D density(1);
        CellStore3D cells;
        CellInit initial;
        initial.uid = 30;
        initial.type = CellType::K;
        const Slot mother = cells.create(initial);
        assert(grid.place_single(initial.anchor, mother));
        density.add(initial.anchor, initial.type, mother);
        CellUid next_uid = 31;
        std::vector<LineageEdge> lineage;
        const DivisionResult result = divide_cell(mother, 4.0, next_uid, cells, grid,
                                                  density, config, lineage);
        assert(result.changed && grid.occupants({0, 0, 0}).size() == 2);
        assert(cells.stage(mother) == CellStage::ultrasmall);
        assert(cells.stage(result.daughter) == CellStage::ultrasmall);

        // With no free neighboring voxel, a later stage-2 recovery attempt is
        // retryable: both co-located cells survive and the failed attempt must
        // not permanently mutate the mother's inherited growth property.
        const float inherent_before_retry = cells.inherent_growth_rate(mother);
        const std::uint64_t sequence_before_retry = cells.event_sequence(mother);
        const CellUid uid_before_retry = next_uid;
        const std::size_t lineage_before_retry = lineage.size();
        const DivisionResult retry = divide_cell(mother, 5.0, next_uid, cells, grid,
                                                 density, config, lineage);
        assert(!retry.changed && !retry.mother_removed && retry.daughter == kEmptySlot);
        assert(cells.valid(mother) && cells.valid(result.daughter));
        assert(grid.occupants({0, 0, 0}).size() == 2);
        assert(cells.stage(mother) == CellStage::ultrasmall);
        assert(cells.stage(result.daughter) == CellStage::ultrasmall);
        assert(cells.inherent_growth_rate(mother) == inherent_before_retry);
        assert(cells.event_sequence(mother) == sequence_before_retry);
        assert(next_uid == uid_before_retry);
        assert(lineage.size() == lineage_before_retry);

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
