#include "rules/density.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_set>

namespace atcg3d {

DensityCounts3D exact_box_counts(const CellStore3D& cells, Vec3i minimum, Vec3i maximum) {
    DensityCounts3D result;
    for (const Slot slot : cells.alive_slots()) {
        const Vec3i anchor = cells.anchor(slot);
        if (anchor.x < minimum.x || anchor.x > maximum.x ||
            anchor.y < minimum.y || anchor.y > maximum.y ||
            anchor.z < minimum.z || anchor.z > maximum.z) {
            continue;
        }
        if (cells.type(slot) == CellType::r) {
            ++result.r;
        } else {
            ++result.K;
        }
    }
    return result;
}

double exact_directional_density(const CellStore3D& cells,
                                 Vec3i anchor,
                                 DirectionId direction,
                                 int radius,
                                 double half_angle_degrees) {
    const auto offsets = directional_cone_offsets(direction, radius, half_angle_degrees);
    if (offsets.empty()) {
        return 1.0;
    }
    std::unordered_set<Vec3i, Vec3iHash> sites;
    sites.reserve(offsets.size());
    for (const Vec3i offset : offsets) {
        sites.insert(anchor + offset);
    }
    std::size_t count = 0;
    for (const Slot slot : cells.alive_slots()) {
        if (sites.contains(cells.anchor(slot))) {
            ++count;
        }
    }
    return static_cast<double>(count) / static_cast<double>(sites.size());
}

DensityGrowthCounts growth_counts(const BlockDensityIndex3D& density,
                                  Vec3i anchor,
                                  int window_edge) {
    const int lower = (window_edge - 1) / 2;
    const int upper = window_edge - lower - 1;
    const DensityCounts3D counts = density.estimate_box(
        anchor - Vec3i{lower, lower, lower}, anchor + Vec3i{upper, upper, upper});
    return {static_cast<long>(counts.r), static_cast<long>(counts.K),
            static_cast<long>(counts.total())};
}

double migration_activation_density(const BlockDensityIndex3D& density,
                                    Vec3i anchor,
                                    CellStage stage,
                                    int window_edge,
                                    int query_block_edge) {
    const DensityCounts3D counts = density.estimate_quantized_box(
        anchor, window_edge, query_block_edge);
    double capacity = static_cast<double>(window_edge) * window_edge * window_edge;
    if (stage == CellStage::large) capacity /= 8.0;
    return capacity > 0.0 ? static_cast<double>(counts.total()) / capacity : 0.0;
}

double density_growth_rate_for_cell(const CellStore3D& cells,
                                    Slot slot,
                                    const BlockDensityIndex3D& density,
                                    const Model3DConfig& config) {
    const DensityGrowthCounts counts = growth_counts(density, cells.anchor(slot), config.growth_density_window_edge);
    return calculate_density_growth_rate(static_cast<int>(cells.type(slot)),
                                         static_cast<double>(cells.inherent_growth_rate(slot)),
                                         counts,
                                         config.r_limit * config.carrying_capacity_scale_2d_to_3d,
                                         config.K_limit * config.carrying_capacity_scale_2d_to_3d,
                                         config.alpha,
                                         config.beta,
                                         config.carrying_capacity_r * config.carrying_capacity_scale_2d_to_3d,
                                         config.carrying_capacity_K * config.carrying_capacity_scale_2d_to_3d);
}

}  // namespace atcg3d
