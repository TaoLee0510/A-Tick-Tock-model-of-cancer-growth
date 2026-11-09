#pragma once

#include <algorithm>
#include <cstdint>

namespace atcg3d::structured_pde {
inline constexpr double kFixed26DiffusionFactor = 9.0 / 26.0;

inline double beta_mean(const BetaRateConfig& config) noexcept {
    return config.scale * config.alpha / (config.alpha + config.beta);
}

inline double positive_part(double value) noexcept {
    return std::max(0.0, value);
}

inline int floor_div(int value, int divisor) noexcept {
    int quotient = value / divisor;
    const int remainder = value % divisor;
    if (remainder != 0 && ((remainder < 0) != (divisor < 0))) --quotient;
    return quotient;
}

inline std::uint64_t bounds_size(const StructuredActiveBounds3D& bounds) noexcept {
    if (!bounds.valid) return 0;
    return static_cast<std::uint64_t>(bounds.x1 - bounds.x0) *
        static_cast<std::uint64_t>(bounds.y1 - bounds.y0) *
        static_cast<std::uint64_t>(bounds.z1 - bounds.z0);
}
inline StructuredActiveBounds3D union_bounds(
    const StructuredActiveBounds3D& lhs,
    const StructuredActiveBounds3D& rhs) noexcept {
    if (!lhs.valid) return rhs;
    if (!rhs.valid) return lhs;
    return {
        std::min(lhs.x0, rhs.x0),
        std::min(lhs.y0, rhs.y0),
        std::min(lhs.z0, rhs.z0),
        std::max(lhs.x1, rhs.x1),
        std::max(lhs.y1, rhs.y1),
        std::max(lhs.z1, rhs.z1),
        true};
}

}  // namespace atcg3d::structured_pde
