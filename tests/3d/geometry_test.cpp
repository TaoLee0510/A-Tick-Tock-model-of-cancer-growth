#include <cassert>
#include <set>

#include "geometry/directions.hpp"
#include "geometry/footprint.hpp"

int main() {
    using namespace atcg3d;
    assert(direction_vector(1) == (Vec3i{-1, -1, 0}));
    assert(direction_vector(8) == (Vec3i{0, -1, 0}));
    assert(direction_vector(17) == (Vec3i{0, 0, -1}));
    assert(direction_vector(26) == (Vec3i{0, 0, 1}));

    std::set<Vec3i> unique;
    for (DirectionId id = 1; id <= 26; ++id) {
        unique.insert(direction_vector(id));
        assert(opposite_direction(opposite_direction(id)) == id);
        if (id <= 8) assert(direction_vector(id).z == 0);
    }
    assert(unique.size() == 26);

    assert(large_footprint({0, 0, 0}).size() == 8);
    assert(entering_voxels({0, 0, 0}, {1, 0, 0}).size() == 4);
    assert(entering_voxels({0, 0, 0}, {1, 1, 0}).size() == 6);
    assert(entering_voxels({0, 0, 0}, {1, 1, 1}).size() == 7);
    assert(stage_recovery_anchors({0, 0, 0}).size() == 8);
    assert(chebyshev_shell({0, 0, 0}, 2).size() == 98);
    assert(shape_reduction_sites({0, 0, 0}).size() == 64);

    const auto turns = turn_neighbors(1, 45.0);
    assert(std::find(turns.begin(), turns.end(), 2) != turns.end());
    assert(std::find(turns.begin(), turns.end(), 8) != turns.end());
    assert(!directional_cone_offsets(1, 5, 45.0).empty());
}
