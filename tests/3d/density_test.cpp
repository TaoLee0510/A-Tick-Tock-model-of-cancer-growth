#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <thread>
#include <vector>

#include "config/model_config.hpp"
#include "core/cell_store.hpp"
#include "rules/density.hpp"
#include "space/density_index.hpp"
#include "space/quantized_box_count_index.hpp"

int main() {
    using namespace atcg3d;
    CellStore3D cells;
    BlockDensityIndex3D exact_index(1);
    for (int index = 0; index < 12; ++index) {
        CellInit cell;
        cell.uid = static_cast<CellUid>(index + 1);
        cell.anchor = {index - 6, index % 3 - 1, index % 5 - 2};
        cell.type = index % 2 == 0 ? CellType::r : CellType::K;
        const Slot slot = cells.create(cell);
        exact_index.add(cells.anchor(slot), cells.type(slot), slot);
    }
    const DensityCounts3D exact_box = exact_box_counts(cells, {-3, -3, -3}, {3, 3, 3});
    const DensityCounts3D indexed_box = exact_index.estimate_box({-3, -3, -3}, {3, 3, 3});
    assert(exact_box.r == indexed_box.r);
    assert(exact_box.K == indexed_box.K);

    const double exact_cone = exact_directional_density(cells, {0, 0, 0}, 6, 5, 45.0);
    const double indexed_cone = exact_index.estimate_directional_density({0, 0, 0}, 6, 5, 45.0);
    assert(std::abs(exact_cone - indexed_cone) < 1e-12);

    BlockDensityIndex3D production(4);
    const double empty = production.estimate_directional_density({0, 0, 0}, 6, 5, 45.0);
    production.add({2, 0, 0}, CellType::r, 20);
    const double one = production.estimate_directional_density({0, 0, 0}, 6, 5, 45.0);
    production.add({3, 0, 0}, CellType::K, 21);
    const double two = production.estimate_directional_density({0, 0, 0}, 6, 5, 45.0);
    // This anchor shares a production block with the cone, but is geometrically
    // behind the queried direction and must not inflate its density.
    production.add({-1, 0, 0}, CellType::K, 22);
    const double outside_same_block =
        production.estimate_directional_density({0, 0, 0}, 6, 5, 45.0);
    assert(empty == 0.0);
    assert(one > empty);
    assert(two > one);
    assert(outside_same_block == two);

    // The production block index scans only anchors in nearby blocks; its
    // geometric test and cached cone volume agree with the brute-force oracle.
    CellStore3D cone_cells;
    BlockDensityIndex3D cone_blocks(4);
    for (int x = -5; x <= 5; ++x) {
        for (int y = -5; y <= 5; ++y) {
            for (int z = -5; z <= 5; ++z) {
                if ((x * 17 + y * 11 + z * 7) % 13 != 0) continue;
                CellInit cell;
                cell.uid = static_cast<CellUid>(cone_cells.slot_count() + 1000);
                cell.anchor = {x, y, z};
                cell.type = (x + y + z) % 2 == 0 ? CellType::r : CellType::K;
                const Slot slot = cone_cells.create(cell);
                cone_blocks.add(cell.anchor, cell.type, slot);
            }
        }
    }
    for (DirectionId direction = 1; direction <= 26; ++direction) {
        const double oracle =
            exact_directional_density(cone_cells, {0, 0, 0}, direction, 5, 45.0);
        const double indexed =
            cone_blocks.estimate_directional_density({0, 0, 0}, direction, 5, 45.0);
        assert(std::abs(oracle - indexed) < 1e-12);
    }
    const auto all_directional =
        cone_blocks.estimate_all_directional_densities(
            {0, 0, 0}, 5, 45.0);
    for (DirectionId direction = 1; direction <= 26; ++direction) {
        const double indexed =
            cone_blocks.estimate_directional_density(
                {0, 0, 0}, direction, 5, 45.0);
        assert(std::abs(indexed - all_directional[direction]) < 1e-12);
    }
    const auto all_thin_directional =
        cone_blocks.estimate_all_directional_densities(
            {0, 0, 0}, 5, 45.0, true);
    for (DirectionId direction = 1; direction <= 26; ++direction) {
        const double indexed =
            cone_blocks.estimate_directional_density(
                {0, 0, 0}, direction, 5, 45.0, true);
        assert(std::abs(indexed - all_thin_directional[direction]) < 1e-12);
    }

    // Complete blocks are aggregated and boundary blocks inspect their compact
    // anchor entries, so a large activation box does not scan cells or 70^3
    // voxels and still excludes anchors just outside its boundary.
    BlockDensityIndex3D activation(4);
    activation.add({-34, 0, 0}, CellType::r, 30);
    activation.add({35, 0, 0}, CellType::K, 31);
    activation.add({-35, 0, 0}, CellType::r, 32);
    activation.add({36, 0, 0}, CellType::K, 33);
    const DensityCounts3D activation_box = activation.estimate_box(
        {-34, -34, -34}, {35, 35, 35});
    assert(activation_box.r == 1 && activation_box.K == 1);
    const double small_activation = migration_activation_density(
        activation, {0, 0, 0}, CellStage::small, 70, 1);
    const double large_activation = migration_activation_density(
        activation, {0, 0, 0}, CellStage::large, 70, 1);
    assert(std::abs(large_activation - 8.0 * small_activation) < 1e-15);

    // Thin-layer mode is an exact 2D rule reduction: directional cones and
    // activation/growth windows inspect only the anchor z plane. Off-plane
    // anchors must not change either numerator or cache invalidation result.
    CellStore3D thin_cells;
    BlockDensityIndex3D thin_density(4);
    const auto add_thin_cell = [&](CellUid uid, Vec3i anchor, CellType type) {
        CellInit cell;
        cell.uid = uid;
        cell.anchor = anchor;
        cell.type = type;
        cell.inherent_growth_rate = 1.2F;
        const Slot slot = thin_cells.create(cell);
        thin_density.add(anchor, type, slot);
        return slot;
    };
    const Slot thin_target = add_thin_cell(100, {0, 0, 0}, CellType::r);
    add_thin_cell(101, {2, 0, 0}, CellType::K);
    add_thin_cell(102, {3, 1, 0}, CellType::r);
    add_thin_cell(103, {2, 0, 1}, CellType::K);
    add_thin_cell(104, {3, 1, -1}, CellType::r);
    const double thin_oracle =
        exact_directional_density(thin_cells, {0, 0, 0}, 6, 5, 45.0, true);
    const double thin_indexed =
        thin_density.estimate_directional_density(
            {0, 0, 0}, 6, 5, 45.0, true);
    assert(std::abs(thin_oracle - thin_indexed) < 1e-12);

    const DensityGrowthCounts thin_growth =
        growth_counts(thin_density, {0, 0, 0}, 6, true);
    assert(thin_growth.rc == 2);
    assert(thin_growth.kc == 1);
    assert(thin_growth.cells_number == 3);
    const double thin_small_activation = migration_activation_density(
        thin_density, {0, 0, 0}, CellStage::small, 70, 32, true);
    const double thin_large_activation = migration_activation_density(
        thin_density, {0, 0, 0}, CellStage::large, 70, 32, true);
    assert(std::abs(thin_small_activation - 3.0 / (70.0 * 70.0)) < 1e-15);
    assert(std::abs(thin_large_activation - 4.0 * thin_small_activation) < 1e-15);

    Model3DConfig thin_config;
    thin_config.thin_layer = true;
    thin_config.growth_density_window_edge = 6;
    const double thin_growth_rate = density_growth_rate_for_cell(
        thin_cells, thin_target, thin_density, thin_config);
    const double legacy_growth_rate = calculate_density_growth_rate(
        static_cast<int>(CellType::r),
        static_cast<double>(thin_cells.inherent_growth_rate(thin_target)),
        thin_growth,
        thin_config.legacy_mapping.source_r_limit,
        thin_config.legacy_mapping.source_K_limit,
        thin_config.alpha,
        thin_config.beta,
        thin_config.legacy_mapping.source_carrying_capacity_r,
        thin_config.legacy_mapping.source_carrying_capacity_K);
    assert(std::abs(thin_growth_rate - legacy_growth_rate) < 1e-12);

    // A large cell contributes one biological anchor, not eight footprint voxels.
    BlockDensityIndex3D unique_anchor(1);
    unique_anchor.add({0, 0, 0}, CellType::r, 40);
    const DensityGrowthCounts counts = growth_counts(unique_anchor, {0, 0, 0}, 6);
    assert(counts.rc == 1);
    assert(counts.cells_number == 1);

    // Slot visitors expose local cells without allocating an intermediate
    // result inside the index, including across negative block coordinates.
    BlockDensityIndex3D visitors(4);
    visitors.add({-1, -1, -1}, CellType::r, 50);
    visitors.add({-2, -2, -2}, CellType::K, 51);
    visitors.add({0, 0, 0}, CellType::r, 52);
    std::vector<Slot> visited_block;
    visitors.for_each_slot_in_block({-1, -1, -1},
                                    [&](Slot slot) { visited_block.push_back(slot); });
    std::sort(visited_block.begin(), visited_block.end());
    assert((visited_block == std::vector<Slot>{50, 51}));
    std::vector<Slot> visited_box;
    visitors.for_each_slot_in_box({-1, -1, -1}, {0, 0, 0},
                                  [&](Slot slot) { visited_box.push_back(slot); });
    std::sort(visited_box.begin(), visited_box.end());
    assert((visited_box == std::vector<Slot>{50, 52}));
    visitors.move({-1, -1, -1}, {5, 5, 5}, CellType::r, 50);
    visited_block.clear();
    visitors.for_each_slot_in_block({1, 1, 1},
                                    [&](Slot slot) { visited_block.push_back(slot); });
    assert((visited_block == std::vector<Slot>{50}));
    visitors.remove({5, 5, 5}, CellType::r, 50);

    // Cached activation boxes are invalidated only when the changed anchor can
    // contribute to that box. Cache identity includes both query dimensions.
    BlockDensityIndex3D cached(1);
    cached.add({0, 0, 0}, CellType::r, 60);
    assert(cached.estimate_quantized_box({0, 0, 0}, 8, 4).total() == 1);
    assert(cached.quantized_cache_size() == 1);
    cached.add({100, 100, 100}, CellType::K, 61);
    assert(cached.quantized_cache_size() == 1);
    cached.add({1, 0, 0}, CellType::K, 62);
    assert(cached.quantized_cache_size() == 0);
    assert(cached.estimate_quantized_box({0, 0, 0}, 8, 4).total() == 2);
    cached.add({2, 2, 2}, CellType::r, 63);
    assert(cached.estimate_quantized_box({0, 0, 0}, 1, 4).total() == 1);
    assert(cached.estimate_quantized_box({0, 0, 0}, 8, 4).total() == 3);
    cached.move({1, 0, 0}, {101, 100, 100}, CellType::K, 62);
    assert(cached.estimate_quantized_box({0, 0, 0}, 8, 4).total() == 2);
    cached.remove({2, 2, 2}, CellType::r, 63);
    assert(cached.estimate_quantized_box({0, 0, 0}, 1, 4).total() == 0);

    BlockDensityIndex3D bounded_cache(1);
    for (std::size_t index = 0;
         index <= BlockDensityIndex3D::quantized_cache_capacity(); ++index) {
        bounded_cache.estimate_quantized_box(
            {static_cast<std::int32_t>(index), 0, 0}, 1, 1);
    }
    assert(bounded_cache.quantized_cache_size() <=
           BlockDensityIndex3D::quantized_cache_capacity());

    // Neighborhood refresh workers share the quantized migration-activation
    // cache. Concurrent cache misses must not race while registering layouts,
    // inserting keys, or clearing at capacity.
    BlockDensityIndex3D concurrent_cache(4);
    concurrent_cache.add({0, 0, 0}, CellType::r, 70);
    std::atomic<bool> start{false};
    std::atomic<int> failures{0};
    std::vector<std::thread> readers;
    for (int worker = 0; worker < 18; ++worker) {
        readers.emplace_back([&, worker] {
            while (!start.load(std::memory_order_acquire)) {
                std::this_thread::yield();
            }
            for (int query = 0; query < 2048; ++query) {
                const Vec3i anchor{query - 1024, worker * 3 - 27,
                                   (query + worker) % 31 - 15};
                const DensityCounts3D value =
                    concurrent_cache.estimate_quantized_box(anchor, 8, 4);
                if (value.total() > 1) {
                    failures.fetch_add(1, std::memory_order_relaxed);
                }
            }
        });
    }
    start.store(true, std::memory_order_release);
    for (std::thread& reader : readers) reader.join();
    assert(failures.load(std::memory_order_relaxed) == 0);
    assert(concurrent_cache.quantized_cache_size() <=
           BlockDensityIndex3D::quantized_cache_capacity());

    // The runtime migration-activation index maintains the same exact
    // quantized box counts incrementally, including negative coordinates and
    // movements across query-block boundaries.
    BlockDensityIndex3D incremental_oracle(4);
    QuantizedBoxCountIndex3D incremental(70, 32, false);
    const std::vector<Vec3i> incremental_anchors{
        {-65, -33, -1}, {-34, 0, 0}, {-1, -1, -1}, {0, 0, 0},
        {31, 31, 31}, {32, 32, 32}, {35, 0, 0}, {70, 2, -40}};
    for (std::size_t index = 0; index < incremental_anchors.size(); ++index) {
        incremental_oracle.add(incremental_anchors[index], CellType::r,
                               static_cast<Slot>(100 + index));
        incremental.add(incremental_anchors[index]);
    }
    for (int x = -3; x <= 3; ++x) {
        for (int y = -2; y <= 2; ++y) {
            for (int z = -2; z <= 2; ++z) {
                const Vec3i block{x, y, z};
                const Vec3i representative{x * 32, y * 32, z * 32};
                assert(incremental.count(block) ==
                       incremental_oracle
                           .estimate_quantized_box(representative, 70, 32)
                           .total());
            }
        }
    }
    incremental_oracle.move({31, 31, 31}, {64, -32, 1}, CellType::r, 104);
    incremental.move({31, 31, 31}, {64, -32, 1});
    incremental_oracle.remove({-34, 0, 0}, CellType::r, 101);
    incremental.remove({-34, 0, 0});
    for (const Vec3i block : incremental.resident_blocks()) {
        const Vec3i representative{block.x * 32, block.y * 32, block.z * 32};
        assert(incremental.count(block) ==
               incremental_oracle
                   .estimate_quantized_box(representative, 70, 32)
                   .total());
        assert(incremental.resident_count(block) > 0);
    }

    // Growth windows use an exact per-slot incremental cache in production.
    // Every mutation must remain identical to the brute block estimator,
    // including sparse slots, negative coordinates, moves, and removals.
    BlockDensityIndex3D local_windows(4);
    local_windows.configure_local_window_counts(6, false);
    std::vector<std::pair<Slot, Vec3i>> local_cells;
    const auto add_local = [&](Slot slot, Vec3i anchor, CellType type) {
        local_windows.add(anchor, type, slot);
        local_cells.push_back({slot, anchor});
    };
    add_local(2, {-4, -3, -2}, CellType::r);
    add_local(5, {-1, -1, 0}, CellType::K);
    add_local(7, {0, 0, 0}, CellType::r);
    add_local(11, {3, 2, 1}, CellType::K);
    const auto verify_local = [&] {
        for (const auto& [slot, anchor] : local_cells) {
            const DensityCounts3D cached_counts =
                local_windows.local_window_counts(slot);
            const DensityCounts3D exact_counts = local_windows.estimate_box(
                anchor - Vec3i{2, 2, 2}, anchor + Vec3i{3, 3, 3});
            assert(cached_counts.r == exact_counts.r);
            assert(cached_counts.K == exact_counts.K);
        }
    };
    verify_local();
    local_windows.move({-1, -1, 0}, {4, -2, 3}, CellType::K, 5);
    local_cells[1].second = {4, -2, 3};
    verify_local();
    local_windows.remove({-4, -3, -2}, CellType::r, 2);
    local_cells.erase(local_cells.begin());
    verify_local();

    BlockDensityIndex3D local_thin(4);
    local_thin.configure_local_window_counts(6, true);
    local_thin.add({0, 0, 0}, CellType::r, 1);
    local_thin.add({1, 1, 1}, CellType::K, 3);
    assert(local_thin.local_window_counts(1).total() == 1);
    assert(local_thin.local_window_counts(3).total() == 1);

    // An axial one-voxel move updates only the symmetric differences of the
    // old/new six-wide windows. Each difference is two 6x6 slabs (72 lattice
    // sites), rather than rescanning two complete 6x6x6 windows (432 sites).
    // The moving cell's translated observation window is updated by the same
    // bounded method, and every cached count remains exact.
    BlockDensityIndex3D axial_move(4);
    axial_move.configure_local_window_counts(6, false);
    struct AxialCell {
        Slot slot{};
        Vec3i anchor{};
        CellType type{CellType::r};
    };
    std::vector<AxialCell> axial_cells;
    axial_cells.push_back({0, {0, 0, 0}, CellType::r});
    axial_move.add({0, 0, 0}, CellType::r, 0);
    Slot next_axial_slot = 1;
    for (int x = -4; x <= 4; ++x) {
        for (int y = -3; y <= 3; ++y) {
            for (int z = -3; z <= 3; ++z) {
                const Vec3i anchor{x, y, z};
                if (anchor == Vec3i{0, 0, 0} ||
                    anchor == Vec3i{1, 0, 0}) {
                    continue;
                }
                const CellType type =
                    (x + 2 * y + 3 * z) % 2 == 0
                        ? CellType::r : CellType::K;
                axial_move.add(anchor, type, next_axial_slot);
                axial_cells.push_back({next_axial_slot, anchor, type});
                ++next_axial_slot;
            }
        }
    }
    axial_move.move({0, 0, 0}, {1, 0, 0}, CellType::r, 0);
    axial_cells.front().anchor = {1, 0, 0};
    const LocalWindowMoveDiagnostics3D axial_diagnostics =
        axial_move.last_local_window_move_diagnostics();
    assert(axial_diagnostics.affected_window_query_sites == 72);
    assert(axial_diagnostics.moving_window_query_sites == 72);
    assert(axial_diagnostics.affected_window_query_sites <
           2ULL * 6ULL * 6ULL * 6ULL);
    assert(axial_diagnostics.affected_slot_visits > 0);
    const auto verify_axial_cells = [&] {
        for (const AxialCell& cell : axial_cells) {
            const DensityCounts3D cached_counts =
                axial_move.local_window_counts(cell.slot);
            const DensityCounts3D exact_counts = axial_move.estimate_box(
                cell.anchor - Vec3i{2, 2, 2},
                cell.anchor + Vec3i{3, 3, 3});
            assert(cached_counts.r == exact_counts.r);
            assert(cached_counts.K == exact_counts.K);
        }
    };
    verify_axial_cells();
    std::vector<Slot> former_anchor_slots;
    axial_move.for_each_slot_in_box(
        {0, 0, 0}, {0, 0, 0},
        [&](Slot slot) { former_anchor_slots.push_back(slot); });
    assert(std::find(former_anchor_slots.begin(), former_anchor_slots.end(), 0) ==
           former_anchor_slots.end());
    std::vector<Slot> target_anchor_slots;
    axial_move.for_each_slot_in_box(
        {1, 0, 0}, {1, 0, 0},
        [&](Slot slot) { target_anchor_slots.push_back(slot); });
    assert((target_anchor_slots == std::vector<Slot>{0}));

    // Plane- and space-diagonal unit moves exercise all disjoint slab axes.
    // Their symmetric-difference volumes are likewise much smaller than two
    // complete windows and retain exact cached counts even with co-location at
    // the destination anchor.
    axial_move.move({1, 0, 0}, {2, 1, 0}, CellType::r, 0);
    axial_cells.front().anchor = {2, 1, 0};
    const LocalWindowMoveDiagnostics3D plane_diagnostics =
        axial_move.last_local_window_move_diagnostics();
    assert(plane_diagnostics.affected_window_query_sites == 132);
    assert(plane_diagnostics.moving_window_query_sites == 132);
    verify_axial_cells();

    axial_move.move({2, 1, 0}, {3, 2, 1}, CellType::r, 0);
    axial_cells.front().anchor = {3, 2, 1};
    const LocalWindowMoveDiagnostics3D spatial_diagnostics =
        axial_move.last_local_window_move_diagnostics();
    assert(spatial_diagnostics.affected_window_query_sites == 182);
    assert(spatial_diagnostics.moving_window_query_sites == 182);
    verify_axial_cells();

    local_thin.move({0, 0, 0}, {1, 0, 0}, CellType::r, 1);
    const LocalWindowMoveDiagnostics3D thin_move_diagnostics =
        local_thin.last_local_window_move_diagnostics();
    assert(thin_move_diagnostics.affected_window_query_sites == 12);
    assert(thin_move_diagnostics.moving_window_query_sites == 12);
    assert(local_thin.local_window_counts(1).total() == 1);
    assert(local_thin.local_window_counts(3).total() == 1);
}
