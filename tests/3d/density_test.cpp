#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "config/model_config.hpp"
#include "core/cell_store.hpp"
#include "rules/density.hpp"
#include "space/density_index.hpp"

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
}
