#pragma once

#include <cstdint>
#include <span>

#include "core/types.hpp"

namespace atcg3d {

struct VesselDirectionParameters3D {
    double forward_half_angle_degrees{45.0};
    double turn_half_angle_degrees{45.0};
    double persistence_probability{0.90};
    double forward_bias{1.0};
    double distance_weight_exponent{0.0};
};

// Selects one of the supplied feasible lattice directions. Candidates must be
// in the forward cone around bias_axis and, after the first step, in the turn
// cone around last_direction. Persistence is applied before weighted random
// selection. Empty candidate sets return DirectionId 0.
DirectionId select_vessel_growth_direction(
    std::span<const DirectionId> feasible,
    Vec3i bias_axis,
    DirectionId last_direction,
    const VesselDirectionParameters3D& parameters,
    std::uint64_t seed,
    std::uint64_t tip_uid,
    std::uint64_t event_sequence);

double vector_angle_degrees(Vec3i lhs, Vec3i rhs) noexcept;

}  // namespace atcg3d
