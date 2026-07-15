#include "vasculature/geometry.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace atcg3d {
namespace {

std::int64_t checked_floor(double value) {
    const double rounded = std::floor(value);
    if (rounded < static_cast<double>(std::numeric_limits<std::int32_t>::min()) ||
        rounded > static_cast<double>(std::numeric_limits<std::int32_t>::max())) {
        throw std::overflow_error("vascular capsule exceeds int32 lattice coordinates");
    }
    return static_cast<std::int64_t>(rounded);
}

std::int64_t checked_ceil(double value) {
    const double rounded = std::ceil(value);
    if (rounded < static_cast<double>(std::numeric_limits<std::int32_t>::min()) ||
        rounded > static_cast<double>(std::numeric_limits<std::int32_t>::max())) {
        throw std::overflow_error("vascular capsule exceeds int32 lattice coordinates");
    }
    return static_cast<std::int64_t>(rounded);
}

}  // namespace

double squared_distance_to_segment(Vec3i point, Vec3i start, Vec3i end) noexcept {
    const double ab_x = static_cast<double>(end.x) - start.x;
    const double ab_y = static_cast<double>(end.y) - start.y;
    const double ab_z = static_cast<double>(end.z) - start.z;
    const double ap_x = static_cast<double>(point.x) - start.x;
    const double ap_y = static_cast<double>(point.y) - start.y;
    const double ap_z = static_cast<double>(point.z) - start.z;
    const double length_squared = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;
    const double projection = length_squared > 0.0
        ? std::clamp((ap_x * ab_x + ap_y * ab_y + ap_z * ab_z) / length_squared,
                     0.0, 1.0)
        : 0.0;
    const double dx = ap_x - projection * ab_x;
    const double dy = ap_y - projection * ab_y;
    const double dz = ap_z - projection * ab_z;
    return dx * dx + dy * dy + dz * dz;
}

double segment_length(Vec3i start, Vec3i end) noexcept {
    const double dx = static_cast<double>(end.x) - start.x;
    const double dy = static_cast<double>(end.y) - start.y;
    const double dz = static_cast<double>(end.z) - start.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

std::vector<Vec3i> rasterize_capsule(Vec3i start,
                                    Vec3i end,
                                    float diameter_voxels) {
    if (!(diameter_voxels > 0.0F) || !std::isfinite(diameter_voxels)) {
        throw std::invalid_argument("vascular diameter must be finite and positive");
    }
    const double radius = static_cast<double>(diameter_voxels) * 0.5;
    const double radius_squared = radius * radius;
    const double tolerance = 1e-12 * std::max(1.0, radius_squared);
    const std::int64_t min_x = checked_floor(static_cast<double>(std::min(start.x, end.x)) - radius);
    const std::int64_t min_y = checked_floor(static_cast<double>(std::min(start.y, end.y)) - radius);
    const std::int64_t min_z = checked_floor(static_cast<double>(std::min(start.z, end.z)) - radius);
    const std::int64_t max_x = checked_ceil(static_cast<double>(std::max(start.x, end.x)) + radius);
    const std::int64_t max_y = checked_ceil(static_cast<double>(std::max(start.y, end.y)) + radius);
    const std::int64_t max_z = checked_ceil(static_cast<double>(std::max(start.z, end.z)) + radius);

    std::vector<Vec3i> result;
    for (std::int64_t x = min_x; x <= max_x; ++x) {
        for (std::int64_t y = min_y; y <= max_y; ++y) {
            for (std::int64_t z = min_z; z <= max_z; ++z) {
                const Vec3i site{static_cast<std::int32_t>(x), static_cast<std::int32_t>(y),
                                 static_cast<std::int32_t>(z)};
                if (squared_distance_to_segment(site, start, end) <=
                    radius_squared + tolerance) {
                    result.push_back(site);
                }
            }
        }
    }
    return result;
}

}  // namespace atcg3d
