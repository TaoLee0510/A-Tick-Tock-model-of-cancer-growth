#pragma once

#include "common/density_growth_rule.hpp"
#include "config/model_config.hpp"
#include "core/cell_store.hpp"
#include "engine/environment.hpp"
#include "geometry/footprint.hpp"
#include "space/density_index.hpp"

namespace atcg3d {

DensityCounts3D exact_box_counts(const CellStore3D& cells, Vec3i minimum, Vec3i maximum);
double exact_directional_density(const CellStore3D& cells,
                                 Vec3i anchor,
                                 DirectionId direction,
                                 int radius,
                                 double half_angle_degrees,
                                 bool thin_layer = false);
DensityGrowthCounts growth_counts(const BlockDensityIndex3D& density,
                                  Vec3i anchor,
                                  int window_edge,
                                  bool thin_layer = false);
DensityGrowthCounts growth_counts(const BlockDensityIndex3D& density,
                                  Slot slot,
                                  Vec3i anchor,
                                  int window_edge,
                                  bool thin_layer = false);
double migration_activation_density(const BlockDensityIndex3D& density,
                                    Vec3i anchor,
                                    CellStage stage,
                                    int window_edge,
                                    int query_block_edge,
                                    bool thin_layer = false);
double density_growth_rate_for_cell(const CellStore3D& cells,
                                    Slot slot,
                                    const BlockDensityIndex3D& density,
                                    const Model3DConfig& config,
                                    const LocalDensityModifier3D* density_modifier = nullptr);

}  // namespace atcg3d
