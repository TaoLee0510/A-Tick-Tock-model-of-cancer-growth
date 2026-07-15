#include "vasculature/growth_rules.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "core/stateless_rng.hpp"
#include "geometry/directions.hpp"

namespace atcg3d {
namespace {

constexpr std::uint64_t kVesselDirectionEvent = 0x56455353454c4449ULL;
constexpr double kAngleTolerance = 1e-10;

double direction_weight(DirectionId direction,
                        Vec3i bias_axis,
                        const VesselDirectionParameters3D& parameters) {
    const Vec3i vector = direction_vector(direction);
    const double vector_length = std::sqrt(static_cast<double>(squared_length(vector)));
    const double bias_length = std::sqrt(static_cast<double>(squared_length(bias_axis)));
    const double cosine = bias_length > 0.0
        ? std::clamp(static_cast<double>(dot(vector, bias_axis)) /
                         (vector_length * bias_length),
                     -1.0, 1.0)
        : 0.0;
    return std::exp(parameters.forward_bias * cosine) *
           std::pow(vector_length, -parameters.distance_weight_exponent);
}

DirectionId weighted_choice(std::span<const DirectionId> candidates,
                            Vec3i bias_axis,
                            const VesselDirectionParameters3D& parameters,
                            std::uint64_t seed,
                            std::uint64_t tip_uid,
                            std::uint64_t event_sequence,
                            std::uint64_t draw) {
    if (candidates.empty()) return kStayDirection;

    double total = 0.0;
    std::vector<double> cumulative;
    cumulative.reserve(candidates.size());
    for (const DirectionId direction : candidates) {
        total += direction_weight(direction, bias_axis, parameters);
        cumulative.push_back(total);
    }
    const double target = rng_unit(seed, tip_uid, kVesselDirectionEvent,
                                   event_sequence, draw) * total;
    const auto iterator = std::lower_bound(cumulative.begin(), cumulative.end(), target);
    const std::size_t index = std::min<std::size_t>(
        static_cast<std::size_t>(std::distance(cumulative.begin(), iterator)),
        candidates.size() - 1U);
    return candidates[index];
}

void validate_parameters(const VesselDirectionParameters3D& parameters) {
    const double values[] = {
        parameters.forward_half_angle_degrees,
        parameters.turn_half_angle_degrees,
        parameters.persistence_probability,
        parameters.forward_bias,
        parameters.distance_weight_exponent,
    };
    for (const double value : values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("vessel direction parameters must be finite");
        }
    }
    if (parameters.forward_half_angle_degrees <= 0.0 ||
        parameters.forward_half_angle_degrees > 180.0 ||
        parameters.turn_half_angle_degrees <= 0.0 ||
        parameters.turn_half_angle_degrees > 180.0 ||
        parameters.persistence_probability < 0.0 ||
        parameters.persistence_probability > 1.0 ||
        parameters.forward_bias < 0.0 || parameters.distance_weight_exponent < 0.0) {
        throw std::invalid_argument("invalid vessel direction parameters");
    }
}

}  // namespace

double vector_angle_degrees(Vec3i lhs, Vec3i rhs) noexcept {
    const int lhs_squared = squared_length(lhs);
    const int rhs_squared = squared_length(rhs);
    if (lhs_squared == 0 || rhs_squared == 0) return 180.0;
    const double cosine = std::clamp(
        static_cast<double>(dot(lhs, rhs)) /
            std::sqrt(static_cast<double>(lhs_squared) * rhs_squared),
        -1.0, 1.0);
    return std::acos(cosine) * 180.0 / std::acos(-1.0);
}

DirectionId select_vessel_growth_direction(
    std::span<const DirectionId> feasible,
    Vec3i bias_axis,
    DirectionId last_direction,
    const VesselDirectionParameters3D& parameters,
    std::uint64_t seed,
    std::uint64_t tip_uid,
    std::uint64_t event_sequence) {
    validate_parameters(parameters);
    if (squared_length(bias_axis) == 0) return kStayDirection;

    std::vector<DirectionId> candidates;
    candidates.reserve(feasible.size());
    for (const DirectionId direction : feasible) {
        if (direction == kStayDirection || direction > 26) continue;
        const Vec3i candidate = direction_vector(direction);
        if (vector_angle_degrees(candidate, bias_axis) >
            parameters.forward_half_angle_degrees + kAngleTolerance) {
            continue;
        }
        if (last_direction != kStayDirection && last_direction <= 26 &&
            direction_angle_degrees(direction, last_direction) >
                parameters.turn_half_angle_degrees + kAngleTolerance) {
            continue;
        }
        candidates.push_back(direction);
    }
    if (candidates.empty()) return kStayDirection;

    const auto previous = std::find(candidates.begin(), candidates.end(), last_direction);
    if (previous != candidates.end()) {
        if (candidates.size() == 1U ||
            rng_unit(seed, tip_uid, kVesselDirectionEvent, event_sequence, 0) <
                parameters.persistence_probability) {
            return last_direction;
        }
        candidates.erase(previous);
        return weighted_choice(candidates, bias_axis, parameters, seed, tip_uid,
                               event_sequence, 1);
    }
    return weighted_choice(candidates, bias_axis, parameters, seed, tip_uid,
                           event_sequence, 0);
}

}  // namespace atcg3d
