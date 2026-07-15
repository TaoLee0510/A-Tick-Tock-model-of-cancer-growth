#pragma once

#include <vector>

#include "core/types.hpp"

namespace atcg3d {

double squared_distance_to_segment(Vec3i point, Vec3i start, Vec3i end) noexcept;
double segment_length(Vec3i start, Vec3i end) noexcept;

// A voxel is part of the capsule when the distance from its integer-valued
// center to the centerline segment is no greater than diameter/2. The result is
// unique and lexicographically ordered, which makes proposals deterministic.
std::vector<Vec3i> rasterize_capsule(Vec3i start,
                                    Vec3i end,
                                    float diameter_voxels);

}  // namespace atcg3d
