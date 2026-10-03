#include "model/moving_tumor_front.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace atcg3d::continuum {

MovingTumorFrontSummary2D build_moving_tumor_front_mask_2d(
    std::span<const double> occupied_fraction,
    int width,
    int height,
    int smoothing_radius_voxels,
    double density_threshold,
    MovingTumorFrontWorkspace2D& workspace,
    std::vector<std::uint8_t>& tumour_mask) {
    if (width <= 0 || height <= 0 || smoothing_radius_voxels < 0 ||
        !(density_threshold > 0.0) || density_threshold > 1.0 ||
        occupied_fraction.size() !=
            static_cast<std::size_t>(width) * height) {
        throw std::invalid_argument("invalid moving tumour-front field");
    }

    const std::size_t sites = static_cast<std::size_t>(width) * height;
    const int pitch = width + 1;
    workspace.prefix.assign(
        static_cast<std::size_t>(pitch) * (height + 1), 0.0);
    workspace.candidate.assign(sites, 0U);
    workspace.exterior.assign(sites, 0U);
    tumour_mask.assign(sites, 0U);
    workspace.largest_component.clear();

    for (int y = 0; y < height; ++y) {
        double row_sum = 0.0;
        for (int x = 0; x < width; ++x) {
            row_sum += occupied_fraction[
                static_cast<std::size_t>(y) * width + x];
            workspace.prefix[
                static_cast<std::size_t>(y + 1) * pitch + x + 1] =
                workspace.prefix[static_cast<std::size_t>(y) * pitch + x + 1] +
                row_sum;
        }
    }

    const int radius = smoothing_radius_voxels;
    const double denominator = static_cast<double>(
        (2 * radius + 1) * (2 * radius + 1));
    const auto rectangle_sum = [&](int x0, int y0, int x1, int y1) {
        return workspace.prefix[static_cast<std::size_t>(y1) * pitch + x1]
            - workspace.prefix[static_cast<std::size_t>(y0) * pitch + x1]
            - workspace.prefix[static_cast<std::size_t>(y1) * pitch + x0]
            + workspace.prefix[static_cast<std::size_t>(y0) * pitch + x0];
    };
    for (int y = 0; y < height; ++y) {
        const int y0 = std::max(0, y - radius);
        const int y1 = std::min(height, y + radius + 1);
        for (int x = 0; x < width; ++x) {
            const int x0 = std::max(0, x - radius);
            const int x1 = std::min(width, x + radius + 1);
            const double smoothed = rectangle_sum(x0, y0, x1, y1) /
                denominator;
            if (smoothed >= density_threshold) {
                workspace.candidate[
                    static_cast<std::size_t>(y) * width + x] = 1U;
            }
        }
    }

    // Retain the largest 8-connected body. Equal-sized components are
    // resolved by row-major discovery order, keeping the result deterministic.
    for (std::size_t seed = 0; seed < sites; ++seed) {
        if (workspace.candidate[seed] != 1U) continue;
        workspace.queue.clear();
        workspace.queue.push_back(seed);
        workspace.candidate[seed] = 2U;
        for (std::size_t cursor = 0; cursor < workspace.queue.size(); ++cursor) {
            const std::size_t here = workspace.queue[cursor];
            const int x = static_cast<int>(here % static_cast<std::size_t>(width));
            const int y = static_cast<int>(here / static_cast<std::size_t>(width));
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dx = -1; dx <= 1; ++dx) {
                    if (dx == 0 && dy == 0) continue;
                    const int nx = x + dx;
                    const int ny = y + dy;
                    if (nx < 0 || nx >= width || ny < 0 || ny >= height) continue;
                    const std::size_t neighbour =
                        static_cast<std::size_t>(ny) * width + nx;
                    if (workspace.candidate[neighbour] != 1U) continue;
                    workspace.candidate[neighbour] = 2U;
                    workspace.queue.push_back(neighbour);
                }
            }
        }
        if (workspace.queue.size() > workspace.largest_component.size()) {
            workspace.largest_component.assign(
                workspace.queue.begin(), workspace.queue.end());
        }
    }
    if (workspace.largest_component.empty()) return {};

    int x0 = width;
    int y0 = height;
    int x1 = 0;
    int y1 = 0;
    for (const std::size_t location : workspace.largest_component) {
        tumour_mask[location] = 1U;
        const int x = static_cast<int>(
            location % static_cast<std::size_t>(width));
        const int y = static_cast<int>(
            location / static_cast<std::size_t>(width));
        x0 = std::min(x0, x);
        y0 = std::min(y0, y);
        x1 = std::max(x1, x + 1);
        y1 = std::max(y1, y + 1);
    }

    // Flood the non-tumour exterior from the component bounding box. Any
    // unvisited non-tumour site is an internal hole and is filled.
    workspace.queue.clear();
    const auto seed_exterior = [&](int x, int y) {
        const std::size_t location = static_cast<std::size_t>(y) * width + x;
        if (tumour_mask[location] != 0U ||
            workspace.exterior[location] != 0U) return;
        workspace.exterior[location] = 1U;
        workspace.queue.push_back(location);
    };
    for (int x = x0; x < x1; ++x) {
        seed_exterior(x, y0);
        seed_exterior(x, y1 - 1);
    }
    for (int y = y0; y < y1; ++y) {
        seed_exterior(x0, y);
        seed_exterior(x1 - 1, y);
    }
    constexpr int kDx[4]{-1, 1, 0, 0};
    constexpr int kDy[4]{0, 0, -1, 1};
    for (std::size_t cursor = 0; cursor < workspace.queue.size(); ++cursor) {
        const std::size_t here = workspace.queue[cursor];
        const int x = static_cast<int>(here % static_cast<std::size_t>(width));
        const int y = static_cast<int>(here / static_cast<std::size_t>(width));
        for (int direction = 0; direction < 4; ++direction) {
            const int nx = x + kDx[direction];
            const int ny = y + kDy[direction];
            if (nx < x0 || nx >= x1 || ny < y0 || ny >= y1) continue;
            const std::size_t neighbour =
                static_cast<std::size_t>(ny) * width + nx;
            if (tumour_mask[neighbour] != 0U ||
                workspace.exterior[neighbour] != 0U) continue;
            workspace.exterior[neighbour] = 1U;
            workspace.queue.push_back(neighbour);
        }
    }

    MovingTumorFrontSummary2D summary;
    for (int y = y0; y < y1; ++y) {
        for (int x = x0; x < x1; ++x) {
            const std::size_t location =
                static_cast<std::size_t>(y) * width + x;
            if (tumour_mask[location] == 0U &&
                workspace.exterior[location] == 0U) {
                tumour_mask[location] = 1U;
            }
            if (tumour_mask[location] == 0U) continue;
            ++summary.tumour_voxels;
            bool front = false;
            for (int direction = 0; direction < 4; ++direction) {
                const int nx = x + kDx[direction];
                const int ny = y + kDy[direction];
                if (nx < 0 || nx >= width || ny < 0 || ny >= height ||
                    tumour_mask[static_cast<std::size_t>(ny) * width + nx] ==
                        0U) {
                    front = true;
                    break;
                }
            }
            if (front) ++summary.front_voxels;
        }
    }
    return summary;
}

}  // namespace atcg3d::continuum
