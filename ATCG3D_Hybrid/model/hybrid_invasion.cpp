#include "model/hybrid_model.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>

namespace atcg3d::hybrid {
namespace {
constexpr std::array<Vec3i, 6> neighbours{{
    {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}}};
}

void HybridModel3D::rebuild_density_prefix() {
    const auto shape = config_.rules.continuum.grid.shape;
    const int nx = shape[0] + 1, ny = shape[1] + 1, nz = shape[2] + 1;
    density_prefix_.assign(std::size_t(nx) * ny * nz, {});
    const auto offset = [nx, ny](int x, int y, int z) {
        return (std::size_t(z) * ny + y) * nx + x;
    };
    for (int z = 1; z < nz; ++z)
        for (int y = 1; y < ny; ++y)
            for (int x = 1; x < nx; ++x) {
                const auto i = pde_->index(x - 1, y - 1, z - 1);
                density_prefix_[offset(x, y, z)] = {
                    own_r_[i], own_K_[i], own_occupied_[i] + pde_->external_occupied_[i]};
            }
    // Separable scans fix the floating-point order independently of threads.
    for (int z = 1; z < nz; ++z)
        for (int y = 1; y < ny; ++y)
            for (int x = 1; x < nx; ++x)
                for (int field = 0; field < 3; ++field)
                    density_prefix_[offset(x, y, z)][field] +=
                        density_prefix_[offset(x - 1, y, z)][field];
    for (int z = 1; z < nz; ++z)
        for (int y = 1; y < ny; ++y)
            for (int x = 1; x < nx; ++x)
                for (int field = 0; field < 3; ++field)
                    density_prefix_[offset(x, y, z)][field] +=
                        density_prefix_[offset(x, y - 1, z)][field];
    for (int z = 1; z < nz; ++z)
        for (int y = 1; y < ny; ++y)
            for (int x = 1; x < nx; ++x)
                for (int field = 0; field < 3; ++field)
                    density_prefix_[offset(x, y, z)][field] +=
                        density_prefix_[offset(x, y, z - 1)][field];
}

std::array<double, 3> HybridModel3D::box_counts(Vec3i point, int lower, int upper) const {
    if (density_prefix_.empty())
        return {};
    const auto& grid = config_.rules.continuum.grid;
    const bool thin = config_.rules.continuum.base.thin_layer;
    const std::array<int, 3> position{point.x, point.y, point.z};
    std::array<int, 3> lo{}, hi{};
    for (int axis = 0; axis < 3; ++axis) {
        if (thin && axis == 2) {
            hi[axis] = 1;
            continue;
        }
        const int center = int(std::floor(position[axis] - grid.origin[axis]));
        lo[axis] = std::clamp(center - lower, 0, grid.shape[axis]);
        hi[axis] = std::clamp(center + upper + 1, 0, grid.shape[axis]);
    }
    const int nx = grid.shape[0] + 1, ny = grid.shape[1] + 1;
    const auto offset = [nx, ny](int x, int y, int z) {
        return (std::size_t(z) * ny + y) * nx + x;
    };
    std::array<double, 3> result{};
    for (int field = 0; field < 3; ++field) {
        const auto value = [&](int x, int y, int z) {
            return density_prefix_[offset(x, y, z)][field];
        };
        result[field] = std::max(0.0,
            ((value(hi[0], hi[1], hi[2]) - value(lo[0], hi[1], hi[2])) -
             (value(hi[0], lo[1], hi[2]) - value(lo[0], lo[1], hi[2]))) -
            ((value(hi[0], hi[1], lo[2]) - value(lo[0], hi[1], lo[2])) -
             (value(hi[0], lo[1], lo[2]) - value(lo[0], lo[1], lo[2]))));
    }
    return result;
}

void HybridModel3D::classify_invasion_core() {
    const auto& grid = config_.rules.continuum.grid;
    const auto& nutrient = config_.rules.continuum.nutrient;
    const bool thin = config_.rules.continuum.base.thin_layer;
    const std::size_t directions = thin ? 4 : 6;
    std::vector<double> occupied(core_.size());
    for (std::size_t i = 0; i < core_.size(); ++i)
        occupied[i] = own_occupied_[i] + pde_->external_occupied_[i];
    if (thin) {
        continuum::build_moving_tumor_front_mask_2d(occupied, grid.shape[0], grid.shape[1],
            nutrient.tumor_front_smoothing_radius_voxels,
            nutrient.tumor_front_density_threshold, front_workspace_, front_mask_);
    } else {
        front_mask_.resize(core_.size());
        const int radius = nutrient.tumor_front_smoothing_radius_voxels;
        for (std::size_t i = 0; i < core_.size(); ++i) {
            const auto point = site(i);
            const auto density = box_counts(point, radius, radius)[2];
            int count = 1;
            const std::array<int, 3> coordinates{point.x, point.y, point.z};
            for (int axis = 0; axis < 3; ++axis) {
                const int center = int(std::floor(coordinates[axis] - grid.origin[axis]));
                count *= std::min(grid.shape[axis] - 1, center + radius) -
                         std::max(0, center - radius) + 1;
            }
            front_mask_[i] = density >= count * nutrient.tumor_front_density_threshold;
        }
    }
    front_distance_.assign(core_.size(), std::numeric_limits<int>::max());
    std::queue<std::size_t> queue;
    for (std::size_t i = 0; i < core_.size(); ++i) {
        if (!front_mask_[i]) {
            front_distance_[i] = 0;
            continue;
        }
        const auto point = site(i);
        for (std::size_t direction = 0; direction < directions; ++direction) {
            const auto j = location(point + neighbours[direction]);
            if (j >= core_.size() || !front_mask_[j]) {
                front_distance_[i] = 1;
                queue.push(i);
                break;
            }
        }
    }
    while (!queue.empty()) {
        const auto i = queue.front();
        queue.pop();
        const auto point = site(i);
        for (std::size_t direction = 0; direction < directions; ++direction) {
            const auto j = location(point + neighbours[direction]);
            if (j < core_.size() && front_mask_[j] &&
                front_distance_[j] > front_distance_[i] + 1) {
                front_distance_[j] = front_distance_[i] + 1;
                queue.push(j);
            }
        }
    }
    active_guard_.assign(core_.size(), 0);
    std::queue<std::pair<std::size_t, int>> active_queue;
    for (auto slot : abm_->cells_.alive_slots()) {
        if (!(abm_->cells_.flags(slot) & kMigrationActive))
            continue;
        for (const auto point : footprint(abm_->cells_.anchor(slot), abm_->cells_.stage(slot))) {
            const auto i = location(point);
            if (i < core_.size() && !active_guard_[i]) {
                active_guard_[i] = 1;
                active_queue.emplace(i, 0);
            }
        }
    }
    while (!active_queue.empty()) {
        const auto [i, distance] = active_queue.front();
        active_queue.pop();
        if (distance == config_.active_guard_voxels)
            continue;
        const auto point = site(i);
        for (std::size_t direction = 0; direction < directions; ++direction) {
            const auto j = location(point + neighbours[direction]);
            if (j < core_.size() && !active_guard_[j]) {
                active_guard_[j] = 1;
                active_queue.emplace(j, distance + 1);
            }
        }
    }
    for (std::size_t i = 0; i < core_.size(); ++i) {
        const auto point = site(i);
        double gradient = 0.0;
        for (std::size_t direction = 0; direction < directions; direction += 2) {
            const auto a = location(point + neighbours[direction]);
            const auto b = location(point + neighbours[direction + 1]);
            const double high = a < core_.size() ? pde_->nutrient_[a] : pde_->nutrient_[i];
            const double low = b < core_.size() ? pde_->nutrient_[b] : pde_->nutrient_[i];
            gradient = std::max(gradient, 0.5 * std::abs(high - low));
        }
        const int radius = config_.smoothing_radius;
        const double sum = box_counts(point, radius, radius)[2];
        int count = 1;
        const std::array<int, 3> coordinates{point.x, point.y, point.z};
        for (int axis = 0; axis < (thin ? 2 : 3); ++axis) {
            const int center = int(std::floor(coordinates[axis] - grid.origin[axis]));
            count *= std::min(grid.shape[axis] - 1, center + radius) -
                     std::max(0, center - radius) + 1;
        }
        const bool retain = core_[i] != 0;
        const int band = config_.front_band_voxels +
            (retain ? 0 : config_.front_hysteresis_voxels);
        core_[i] = !active_guard_[i] && front_distance_[i] > band &&
            sum >= count * (retain ? config_.core_off : config_.core_on) &&
            gradient < (retain ? config_.gradient_on : config_.gradient_off);
    }
}

bool HybridModel3D::is_core(Vec3i point) const {
    const auto i = location(point);
    return i < core_.size() && core_[i];
}

void HybridModel3D::update_representation_coverage() {
    const auto d = diagnostics();
    if (d.total_mass > 0.0) {
        minimum_abm_fraction_ = std::min(minimum_abm_fraction_, d.abm_mass / d.total_mass);
        maximum_pde_fraction_ = std::max(maximum_pde_fraction_, d.pde_mass / d.total_mass);
    } else {
        minimum_abm_fraction_ = 0.0;
    }
    if (d.r_mass > 0.0)
        minimum_active_fraction_ = std::min(minimum_active_fraction_, d.active_mass / d.r_mass);
    else
        minimum_active_fraction_ = 0.0;
}

void HybridModel3D::run_abm_until(double target_hours) {
    if (config_.mode != "all_abm" || target_hours < time_hours() ||
        target_hours > config_.rules.continuum.end_time_hours)
        throw std::invalid_argument("invalid hybrid native-ABM sample barrier");
    abm_->run([&](const Simulation3D& current) {
        if (current.clock().time_hours + 1e-9 >= target_hours)
            abm_->stop_requested_ = true;
    });
    abm_->stop_requested_ = false;
    if (config_.schema_version >= 4)
        update_representation_coverage();
}
} // namespace atcg3d::hybrid
