#pragma once

#include <algorithm>
#include <cmath>

namespace atcg3d::continuum {

// Implicit local Michaelis-Menten uptake after an explicit diffusion step.
// This nonnegative root is shared by the transient ABM and both field models.
inline double resource_after_uptake(double diffused, double consumers,
                                    double rate, double half, double dt,
                                    double maximum) noexcept {
    const double demand = rate * consumers;
    const double b = half + dt * demand - diffused;
    const double discriminant = std::max(0.0, b * b + 4.0 * half * diffused);
    return std::clamp(0.5 * (-b + std::sqrt(discriminant)), 0.0, maximum);
}

}  // namespace atcg3d::continuum
