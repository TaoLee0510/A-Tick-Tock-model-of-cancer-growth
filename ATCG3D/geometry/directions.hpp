#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "core/types.hpp"

namespace atcg3d {

inline constexpr std::array<Vec3i, 27> kDirectionVectors{{
    {0, 0, 0},
    {-1, -1, 0}, {-1, 0, 0}, {-1, 1, 0}, {0, 1, 0},
    {1, 1, 0}, {1, 0, 0}, {1, -1, 0}, {0, -1, 0},
    {-1, -1, -1}, {-1, 0, -1}, {-1, 1, -1}, {0, 1, -1},
    {1, 1, -1}, {1, 0, -1}, {1, -1, -1}, {0, -1, -1}, {0, 0, -1},
    {-1, -1, 1}, {-1, 0, 1}, {-1, 1, 1}, {0, 1, 1},
    {1, 1, 1}, {1, 0, 1}, {1, -1, 1}, {0, -1, 1}, {0, 0, 1},
}};

inline constexpr Vec3i direction_vector(DirectionId id) {
    return kDirectionVectors[id <= 26 ? id : 0];
}

inline constexpr int squared_length(Vec3i value) {
    return value.x * value.x + value.y * value.y + value.z * value.z;
}

inline constexpr int dot(Vec3i lhs, Vec3i rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

inline double direction_angle_degrees(DirectionId lhs, DirectionId rhs) {
    const Vec3i a = direction_vector(lhs);
    const Vec3i b = direction_vector(rhs);
    const int a2 = squared_length(a);
    const int b2 = squared_length(b);
    if (a2 == 0 || b2 == 0) {
        return 180.0;
    }
    const double cosine = std::clamp(
        static_cast<double>(dot(a, b)) / std::sqrt(static_cast<double>(a2 * b2)), -1.0, 1.0);
    return std::acos(cosine) * 180.0 / std::acos(-1.0);
}

inline DirectionId opposite_direction(DirectionId id) {
    if (id == 0 || id > 26) {
        return 0;
    }
    const Vec3i wanted{-direction_vector(id).x, -direction_vector(id).y, -direction_vector(id).z};
    for (DirectionId candidate = 1; candidate <= 26; ++candidate) {
        if (direction_vector(candidate) == wanted) {
            return candidate;
        }
    }
    return 0;
}

inline std::vector<DirectionId> turn_neighbors(DirectionId previous, double half_angle_degrees) {
    std::vector<DirectionId> result;
    if (previous == 0 || previous > 26) {
        return result;
    }
    constexpr double tolerance = 1e-10;
    for (DirectionId candidate = 1; candidate <= 26; ++candidate) {
        if (candidate != previous &&
            direction_angle_degrees(previous, candidate) <= half_angle_degrees + tolerance) {
            result.push_back(candidate);
        }
    }
    return result;
}

}  // namespace atcg3d
