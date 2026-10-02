#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace atcg3d::continuum {

struct MovingTumorFrontWorkspace2D {
    std::vector<double> prefix;
    std::vector<std::uint8_t> candidate;
    std::vector<std::uint8_t> exterior;
    std::vector<std::size_t> queue;
    std::vector<std::size_t> largest_component;
};

struct MovingTumorFrontSummary2D {
    std::size_t tumour_voxels{};
    std::size_t front_voxels{};
};

// Builds one deterministic, hole-filled tumour body from a local 2D occupied
// fraction field. Eight-neighbour connectivity keeps diagonally continuous
// invasion paths together; four-neighbour exterior flooding prevents enclosed
// low-density holes from becoming artificial nutrient reservoirs.
MovingTumorFrontSummary2D build_moving_tumor_front_mask_2d(
    std::span<const double> occupied_fraction,
    int width,
    int height,
    int smoothing_radius_voxels,
    double density_threshold,
    MovingTumorFrontWorkspace2D& workspace,
    std::vector<std::uint8_t>& tumour_mask);

}  // namespace atcg3d::continuum
