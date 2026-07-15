#include <array>
#include <cassert>
#include <cmath>

#include "geometry/directions.hpp"
#include "vasculature/growth_rules.hpp"

int main() {
    using namespace atcg3d;

    assert(std::abs(vector_angle_degrees({1, 0, 0}, {2, 0, 0})) < 1e-12);
    assert(std::abs(vector_angle_degrees({1, 0, 0}, {-1, 0, 0}) - 180.0) < 1e-12);

    std::array<DirectionId, 26> all{};
    for (DirectionId id = 1; id <= 26; ++id) all[id - 1] = id;

    VesselDirectionParameters3D parameters;
    parameters.forward_half_angle_degrees = 45.0;
    parameters.turn_half_angle_degrees = 180.0;
    parameters.persistence_probability = 0.0;
    for (std::uint64_t sequence = 0; sequence < 1000; ++sequence) {
        const DirectionId selected = select_vessel_growth_direction(
            all, {10, 0, 0}, kStayDirection, parameters, 9, 17, sequence);
        assert(selected != kStayDirection);
        assert(vector_angle_degrees(direction_vector(selected), {1, 0, 0}) <= 45.0 + 1e-10);
    }

    const DirectionId positive_x = 6;
    parameters.persistence_probability = 1.0;
    for (std::uint64_t sequence = 0; sequence < 100; ++sequence) {
        assert(select_vessel_growth_direction(all, {1, 0, 0}, positive_x,
                                              parameters, 11, 99, sequence) == positive_x);
    }

    const std::array<DirectionId, 1> backward{{2}};
    assert(select_vessel_growth_direction(backward, {1, 0, 0}, kStayDirection,
                                          parameters, 1, 1, 1) == kStayDirection);
    assert(select_vessel_growth_direction(all, {}, kStayDirection,
                                          parameters, 1, 1, 1) == kStayDirection);

    const DirectionId first = select_vessel_growth_direction(
        all, {1, 1, 1}, kStayDirection, parameters, 123, 456, 789);
    const DirectionId second = select_vessel_growth_direction(
        all, {1, 1, 1}, kStayDirection, parameters, 123, 456, 789);
    assert(first == second);

    return 0;
}
