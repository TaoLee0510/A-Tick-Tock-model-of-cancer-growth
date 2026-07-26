#include "rules/density.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_set>

#include "vasculature/influence_field.hpp"

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
                                 double half_angle_degrees,
                                 bool thin_layer) {
    const auto offsets = directional_cone_offsets(direction, radius, half_angle_degrees);
    if (offsets.empty()) {
        return 1.0;
    }
    std::unordered_set<Vec3i, Vec3iHash> sites;
    sites.reserve(offsets.size());
    for (const Vec3i offset : offsets) {
        if (thin_layer && offset.z != 0) continue;
        sites.insert(anchor + offset);
    }
    if (sites.empty()) {
        return 1.0;
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
                                  int window_edge,
                                  bool thin_layer) {
    const int lower = (window_edge - 1) / 2;
    const int upper = window_edge - lower - 1;
    const DensityCounts3D counts = density.estimate_box(
        anchor - Vec3i{lower, lower, thin_layer ? 0 : lower},
        anchor + Vec3i{upper, upper, thin_layer ? 0 : upper});
    return {static_cast<long>(counts.r), static_cast<long>(counts.K),
            static_cast<long>(counts.total())};
}

DensityGrowthCounts growth_counts(const BlockDensityIndex3D& density,
                                  Slot slot,
                                  Vec3i anchor,
                                  int window_edge,
                                  bool thin_layer) {
    if (!density.has_local_window_counts(window_edge, thin_layer)) {
        return growth_counts(density, anchor, window_edge, thin_layer);
    }
    const DensityCounts3D counts = density.local_window_counts(slot);
    return {static_cast<long>(counts.r), static_cast<long>(counts.K),
            static_cast<long>(counts.total())};
}

double migration_activation_density(const BlockDensityIndex3D& density,
                                    Vec3i anchor,
                                    CellStage stage,
                                    int window_edge,
                                    int query_block_edge,
                                    bool thin_layer) {
    const DensityCounts3D counts = density.estimate_quantized_box(
        anchor, window_edge, query_block_edge, thin_layer);
    double capacity = static_cast<double>(window_edge) * window_edge;
    if (!thin_layer) capacity *= window_edge;
    if (stage == CellStage::large) capacity /= thin_layer ? 4.0 : 8.0;
    return capacity > 0.0 ? static_cast<double>(counts.total()) / capacity : 0.0;
}

double density_growth_rate_for_cell(const CellStore3D& cells,
                                    Slot slot,
                                    const BlockDensityIndex3D& density,
                                    const Model3DConfig& config,
                                    const VascularInfluenceField3D* vascular_influence) {
    const DensityGrowthCounts counts = growth_counts(
        density, slot, cells.anchor(slot), config.growth_density_window_edge,
        config.thin_layer);
    const double retained_density = vascular_influence == nullptr
        ? 1.0
        : 1.0 - static_cast<double>(vascular_influence->relief(cells.anchor(slot)));
    const double r_limit = config.thin_layer
        ? config.legacy_mapping.source_r_limit
        : config.r_limit * config.carrying_capacity_scale_2d_to_3d;
    const double K_limit = config.thin_layer
        ? config.legacy_mapping.source_K_limit
        : config.K_limit * config.carrying_capacity_scale_2d_to_3d;
    const double carrying_capacity_r = config.thin_layer
        ? config.legacy_mapping.source_carrying_capacity_r
        : config.carrying_capacity_r * config.carrying_capacity_scale_2d_to_3d;
    const double carrying_capacity_K = config.thin_layer
        ? config.legacy_mapping.source_carrying_capacity_K
        : config.carrying_capacity_K * config.carrying_capacity_scale_2d_to_3d;
    return calculate_density_growth_rate_continuous(
        static_cast<int>(cells.type(slot)),
        static_cast<double>(cells.inherent_growth_rate(slot)),
        static_cast<double>(counts.rc) * retained_density,
        static_cast<double>(counts.kc) * retained_density,
        static_cast<double>(counts.cells_number) * retained_density,
        r_limit,
        K_limit,
        config.alpha,
        config.beta,
        carrying_capacity_r,
        carrying_capacity_K);
}

}  // namespace atcg3d
