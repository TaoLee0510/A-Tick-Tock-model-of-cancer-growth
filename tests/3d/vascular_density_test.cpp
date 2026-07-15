#include <cassert>
#include <cmath>

#include "config/model_config.hpp"
#include "core/cell_store.hpp"
#include "rules/density.hpp"
#include "space/density_index.hpp"
#include "vasculature/influence_field.hpp"

int main() {
    using namespace atcg3d;

    DensityGrowthCounts integer_counts{17, 23, 40};
    const double legacy = calculate_density_growth_rate(
        1, 1.25, integer_counts, 15.5, 18.0, 2.2, 0.0, 31.0, 36.0);
    const double continuous = calculate_density_growth_rate_continuous(
        1, 1.25, 17.0, 23.0, 40.0, 15.5, 18.0, 2.2, 0.0, 31.0, 36.0);
    assert(legacy == continuous);

    Model3DConfig config;
    CellStore3D cells;
    BlockDensityIndex3D density(4);
    CellUid uid = 1;
    CellInit target_cell;
    target_cell.anchor = {0, 0, 0};
    target_cell.uid = uid++;
    target_cell.clone_id = 1;
    target_cell.type = CellType::r;
    target_cell.stage = CellStage::small;
    target_cell.inherent_growth_rate = 1.0F;
    const Slot target = cells.create(target_cell);
    density.add(target_cell.anchor, target_cell.type, target);

    int added = 1;
    for (int z = -2; z <= 3 && added < 120; ++z) {
        for (int y = -2; y <= 3 && added < 120; ++y) {
            for (int x = -2; x <= 3 && added < 120; ++x) {
                if (x == 0 && y == 0 && z == 0) continue;
                CellInit cell;
                cell.anchor = {x, y, z};
                cell.uid = uid++;
                cell.clone_id = static_cast<std::uint32_t>(cell.uid);
                cell.type = CellType::K;
                cell.stage = CellStage::small;
                cell.inherent_growth_rate = 1.0F;
                const Slot slot = cells.create(cell);
                density.add(cell.anchor, cell.type, slot);
                ++added;
            }
        }
    }
    const double raw = density_growth_rate_for_cell(cells, target, density, config);

    VascularInfluenceField3D influence(4, 2.0F, 0.8F);
    influence.add_source(cells.anchor(target));
    const double relieved = density_growth_rate_for_cell(cells, target, density, config, &influence);
    assert(relieved > raw);
    assert(std::abs(relieved - cells.inherent_growth_rate(target)) < 1e-6);

    VascularInfluenceField3D exponential(
        4, 5.0F, 0.8F, VascularInfluenceProfile3D::exponential, 1.0F);
    exponential.add_source({0, 0, 0});
    assert(std::abs(exponential.relief({0, 0, 0}) - 0.8F) < 1e-6F);
    assert(exponential.relief({1, 0, 0}) > exponential.relief({2, 0, 0}));
    assert(exponential.relief({5, 0, 0}) == 0.0F);

    return 0;
}
