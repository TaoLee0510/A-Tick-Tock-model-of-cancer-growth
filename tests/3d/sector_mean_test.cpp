#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "geometry/directions.hpp"
#include "model/sector_mean_field.hpp"

namespace {
using namespace atcg3d;

std::pair<double, std::uint32_t> direct_mean(
    const std::vector<double>& field, std::array<int, 3> shape,
    Vec3i point, DirectionId direction, int edge, bool thin) {
    const auto forward = direction_vector(direction);
    const int lower = (edge - 1) / 2, upper = edge - lower - 1;
    double sum = 0.0;
    std::uint32_t count = 0;
    for (int z = thin ? 0 : -lower; z <= (thin ? 0 : upper); ++z) {
        for (int y = -lower; y <= upper; ++y) {
            for (int x = -lower; x <= upper; ++x) {
                const Vec3i offset{x, y, z};
                if (squared_length(offset) == 0) continue;
                const double cosine = double(dot(offset, forward)) /
                    std::sqrt(double(squared_length(offset) * squared_length(forward)));
                if (cosine + 1e-12 < std::cos(std::acos(-1.0) / 4.0)) continue;
                const auto target = point + offset;
                if (target.x < 0 || target.x >= shape[0] || target.y < 0 ||
                    target.y >= shape[1] || target.z < 0 || target.z >= shape[2]) continue;
                sum += field[(std::size_t(target.z) * shape[1] + target.y) * shape[0] + target.x];
                ++count;
            }
        }
    }
    return {count ? sum / count : 0.0, count};
}

double check_geometry(bool thin, int edge, int extent) {
    const std::array<int, 3> shape{extent, extent, thin ? 1 : extent};
    const std::size_t size = std::size_t(shape[0]) * shape[1] * shape[2];
    std::vector<double> field(size);
    for (std::size_t i = 0; i < size; ++i)
        field[i] = double((i * 131 + 17) % 1009) / 1009.0;
    continuum::SectorMeanField3D cache(shape, thin, edge, 45.0, 1);
    const int tile = std::min(32, extent);
    const std::array<continuum::SectorQueryBox3D, 2> boxes{{
        {{0, 0, 0}, {tile, tile, thin ? 1 : tile}},
        {{extent - tile, extent - tile, thin ? 0 : extent - tile}, shape}}};
    cache.prepare(field, 1.0, boxes, 1);
    const std::array<Vec3i, 4> points{{{0, 0, 0}, {1, 1, thin ? 0 : 1},
        {extent - 1, extent - 1, thin ? 0 : extent - 1},
        {extent - 2, extent - 3, thin ? 0 : extent - 2}}};
    std::vector<double> serial;
    double maximum_error = 0.0;
    for (const auto point : points) {
        const auto i = (std::size_t(point.z) * extent + point.y) * extent + point.x;
        for (DirectionId direction = 1; direction <= 26; ++direction) {
            if (thin && direction_vector(direction).z != 0) {
                assert(cache.count(i, direction) == 0);
                continue;
            }
            const auto [expected, count] = direct_mean(field, shape, point, direction, edge, thin);
            assert(cache.count(i, direction) == count);
            const double actual = cache.mean(i, direction);
            maximum_error = std::max(maximum_error, std::abs(actual - expected));
            assert(std::abs(actual - expected) < 5e-4);
            serial.push_back(actual);
        }
    }
    cache.set_threads(8);
    cache.prepare(field, 1.0, boxes, 2);
    std::size_t query = 0;
    for (const auto point : points) {
        const auto i = (std::size_t(point.z) * extent + point.y) * extent + point.x;
        for (DirectionId direction = 1; direction <= 26; ++direction) {
            if (!thin || direction_vector(direction).z == 0)
                assert(cache.mean(i, direction) == serial[query++]);
        }
    }
    // A new field epoch must replace means while reusing geometric counts.
    std::fill(field.begin(), field.end(), 0.25);
    cache.prepare(field, 1.0, boxes, 3);
    for (const auto point : points) {
        const auto i = (std::size_t(point.z) * extent + point.y) * extent + point.x;
        for (DirectionId direction = 1; direction <= 26; ++direction)
            if (cache.count(i, direction) > 0)
                assert(std::abs(cache.mean(i, direction) - 0.25) < 5e-4);
    }
    if (extent > 64) {
        bool rejected = false;
        try {
            cache.mean((std::size_t(extent / 2) * extent + extent / 2) * extent + extent / 2, 1);
        } catch (const std::logic_error&) {
            rejected = true;
        }
        assert(rejected);
    }
    std::cout << std::setprecision(17) << "thin=" << thin << " edge=" << edge
              << " extent=" << extent << " maximum_mean_error=" << maximum_error
              << " allocated_bytes=" << cache.allocated_bytes() << '\n';
    return maximum_error;
}

void check_orientation() {
    const std::array<int, 3> shape{17, 17, 17};
    std::vector<double> field(17 * 17 * 17);
    for (std::size_t i = 0; i < field.size(); ++i)
        field[i] = double(i % 17) / 16.0;
    continuum::SectorMeanField3D cache(shape, false, 8, 45.0, 4);
    const std::array<continuum::SectorQueryBox3D, 1> boxes{{{{0, 0, 0}, shape}}};
    cache.prepare(field, 1.0, boxes);
    const auto center = (8 * 17 + 8) * 17 + 8;
    for (DirectionId direction = 1; direction <= 26; ++direction) {
        const auto [expected, count] = direct_mean(field, shape, {8, 8, 8}, direction, 8, false);
        assert(cache.count(center, direction) == count);
        assert(std::abs(cache.mean(center, direction) - expected) < 5e-4);
        const auto forward = direction_vector(direction);
        if (forward.x > 0) assert(cache.mean(center, direction) > 0.5);
        if (forward.x < 0) assert(cache.mean(center, direction) < 0.5);
    }
}
}  // namespace

int main() {
    check_geometry(true, 8, 16);
    check_geometry(false, 8, 16);
    check_geometry(false, 70, 17);
    check_geometry(false, 70, 128);
    check_orientation();
}
