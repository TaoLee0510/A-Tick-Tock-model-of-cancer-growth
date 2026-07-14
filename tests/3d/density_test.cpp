#include <cassert>
#include <cmath>

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
        exact_index.add(cells.anchor(slot), cells.type(slot));
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
    production.add({2, 0, 0}, CellType::r);
    const double one = production.estimate_directional_density({0, 0, 0}, 6, 5, 45.0);
    production.add({3, 0, 0}, CellType::K);
    const double two = production.estimate_directional_density({0, 0, 0}, 6, 5, 45.0);
    assert(empty == 0.0);
    assert(one > empty);
    assert(two > one);

    // Complete blocks are aggregated and boundary blocks inspect their compact
    // anchor entries, so a large activation box does not scan cells or 70^3
    // voxels and still excludes anchors just outside its boundary.
    BlockDensityIndex3D activation(4);
    activation.add({-34, 0, 0}, CellType::r);
    activation.add({35, 0, 0}, CellType::K);
    activation.add({-35, 0, 0}, CellType::r);
    activation.add({36, 0, 0}, CellType::K);
    const DensityCounts3D activation_box = activation.estimate_box(
        {-34, -34, -34}, {35, 35, 35});
    assert(activation_box.r == 1 && activation_box.K == 1);
    const double small_activation = migration_activation_density(
        activation, {0, 0, 0}, CellStage::small, 70, 1);
    const double large_activation = migration_activation_density(
        activation, {0, 0, 0}, CellStage::large, 70, 1);
    assert(std::abs(large_activation - 8.0 * small_activation) < 1e-15);

    // A large cell contributes one biological anchor, not eight footprint voxels.
    BlockDensityIndex3D unique_anchor(1);
    unique_anchor.add({0, 0, 0}, CellType::r);
    const DensityGrowthCounts counts = growth_counts(unique_anchor, {0, 0, 0}, 6);
    assert(counts.r_count == 1);
    assert(counts.all_count == 1);
}
