#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "core/types.hpp"
#include "geometry/directions.hpp"

namespace atcg3d {

inline constexpr std::array<Vec3i, 8> kLargeFootprintOffsets{{
    {0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1},
    {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1},
}};

inline std::array<Vec3i, 8> large_footprint(Vec3i anchor) {
    std::array<Vec3i, 8> result{};
    for (std::size_t index = 0; index < result.size(); ++index) {
        result[index] = anchor + kLargeFootprintOffsets[index];
    }
    return result;
}

inline std::vector<Vec3i> entering_voxels(Vec3i anchor, Vec3i displacement) {
    const auto before = large_footprint(anchor);
    const auto after = large_footprint(anchor + displacement);
    std::vector<Vec3i> result;
    result.reserve(7);
    for (const Vec3i site : after) {
        if (std::find(before.begin(), before.end(), site) == before.end()) {
            result.push_back(site);
        }
    }
    return result;
}

inline std::vector<Vec3i> stage_recovery_anchors(Vec3i contained_site) {
    std::vector<Vec3i> result;
    result.reserve(8);
    for (int dx = -1; dx <= 0; ++dx) {
        for (int dy = -1; dy <= 0; ++dy) {
            for (int dz = -1; dz <= 0; ++dz) {
                result.push_back(contained_site + Vec3i{dx, dy, dz});
            }
        }
    }
    return result;
}

inline std::vector<Vec3i> chebyshev_shell(Vec3i center, int radius) {
    std::vector<Vec3i> result;
    if (radius <= 0) {
        return result;
    }
    const int edge = 2 * radius + 1;
    result.reserve(static_cast<std::size_t>(edge * edge * edge - (edge - 2) * (edge - 2) * (edge - 2)));
    for (int dx = -radius; dx <= radius; ++dx) {
        for (int dy = -radius; dy <= radius; ++dy) {
            for (int dz = -radius; dz <= radius; ++dz) {
                if (std::max({std::abs(dx), std::abs(dy), std::abs(dz)}) == radius) {
                    result.push_back(center + Vec3i{dx, dy, dz});
                }
            }
        }
    }
    return result;
}

inline std::vector<Vec3i> shape_reduction_sites(Vec3i anchor) {
    std::vector<Vec3i> result;
    result.reserve(64);
    for (int dx = -1; dx <= 2; ++dx) {
        for (int dy = -1; dy <= 2; ++dy) {
            for (int dz = -1; dz <= 2; ++dz) {
                result.push_back(anchor + Vec3i{dx, dy, dz});
            }
        }
    }
    return result;
}

inline int chebyshev_distance(Vec3i lhs, Vec3i rhs) {
    return std::max({std::abs(lhs.x - rhs.x), std::abs(lhs.y - rhs.y), std::abs(lhs.z - rhs.z)});
}

inline std::vector<Vec3i> directional_cone_offsets(DirectionId direction,
                                                    int radius,
                                                    double half_angle_degrees) {
    std::vector<Vec3i> result;
    if (direction == 0 || direction > 26 || radius <= 0) {
        return result;
    }
    const Vec3i forward = direction_vector(direction);
    const double forward_length = std::sqrt(static_cast<double>(squared_length(forward)));
    const double minimum_cosine = std::cos(half_angle_degrees * std::acos(-1.0) / 180.0);
    for (int dx = -radius; dx <= radius; ++dx) {
        for (int dy = -radius; dy <= radius; ++dy) {
            for (int dz = -radius; dz <= radius; ++dz) {
                const Vec3i offset{dx, dy, dz};
                const int distance = std::max({std::abs(dx), std::abs(dy), std::abs(dz)});
                if (distance == 0 || distance > radius) {
                    continue;
                }
                const double offset_length = std::sqrt(static_cast<double>(squared_length(offset)));
                const double cosine = static_cast<double>(dot(offset, forward)) /
                                      (offset_length * forward_length);
                if (cosine + 1e-12 >= minimum_cosine) {
                    result.push_back(offset);
                }
            }
        }
    }
    return result;
}

}  // namespace atcg3d
