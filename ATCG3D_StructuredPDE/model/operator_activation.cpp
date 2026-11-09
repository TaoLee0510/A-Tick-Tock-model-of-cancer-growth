#include "model/structured_pde_model.hpp"
#include "model/beta_duration.hpp"
#include "model/feasible_jump.hpp"
#include "model/shared_resource.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <utility>

#include "common/density_growth_rule.hpp"
#include "engine/parallelism.hpp"
#include "engine/simulation.hpp"
#include "geometry/directions.hpp"
#include "geometry/footprint.hpp"
#include "rules/initial_rates.hpp"

#include "model/operator_math.hpp"

namespace atcg3d::structured_pde {

void StructuredPdeModel3D::build_activation_density() {
    const auto& continuum = config_.continuum;
    const int nx = continuum.grid.shape[0];
    const int ny = continuum.grid.shape[1];
    const int nz = continuum.grid.shape[2];
    const int edge = continuum.base.migration_activation_window_edge;
    const int query = continuum.base.migration_activation_block_edge;
    const int lower = (edge - 1) / 2;
    const int upper = edge - lower - 1;
    const int origin_x = static_cast<int>(std::floor(continuum.grid.origin[0]));
    const int origin_y = static_cast<int>(std::floor(continuum.grid.origin[1]));
    const int origin_z = static_cast<int>(std::floor(continuum.grid.origin[2]));

    if (continuum.base.thin_layer) {
        const auto clear_density = [&](const StructuredActiveBounds3D& bounds) {
            if (!bounds.valid) return;
            for (int y = bounds.y0; y < bounds.y1; ++y) {
                const std::size_t begin = index(bounds.x0, y, 0);
                const std::size_t end = index(bounds.x1 - 1, y, 0) + 1;
                for (auto& field : activation_density_) {
                    std::fill(field.begin() + begin, field.begin() + end, 0.0);
                }
            }
        };
        clear_density(activation_density_bounds_);
        activation_density_bounds_ = {};
        if (!population_bounds_.valid) return;

        // The ABM performs one query per 32x32 anchor block. Only blocks that
        // contain population can activate r mass, so use an exact local
        // integral image covering those blocks and their 70x70 query windows.
        const int first_qx = floor_div(
            origin_x + population_bounds_.x0, query);
        const int last_qx = floor_div(
            origin_x + population_bounds_.x1 - 1, query);
        const int first_qy = floor_div(
            origin_y + population_bounds_.y0, query);
        const int last_qy = floor_div(
            origin_y + population_bounds_.y1 - 1, query);
        const int source_x0 = std::clamp(
            first_qx * query + query / 2 - origin_x - lower, 0, nx);
        const int source_y0 = std::clamp(
            first_qy * query + query / 2 - origin_y - lower, 0, ny);
        const int source_x1 = std::clamp(
            last_qx * query + query / 2 - origin_x + upper + 1, 0, nx);
        const int source_y1 = std::clamp(
            last_qy * query + query / 2 - origin_y + upper + 1, 0, ny);
        const int width = source_x1 - source_x0;
        const int height = source_y1 - source_y0;
        const int pitch = width + 1;
        std::vector<double> prefix(
            static_cast<std::size_t>(pitch) * (height + 1), 0.0);
        for (int y = 1; y <= height; ++y) {
            double row = 0.0;
            for (int x = 1; x <= width; ++x) {
                const std::size_t source =
                    index(source_x0 + x - 1, source_y0 + y - 1, 0);
                row += (r_normal_[0][source] + r_normal_[1][source] +
                        r_active(StructuredStage3D::small, source) +
                        r_active(StructuredStage3D::large, source) +
                        K_[0][source] + K_[1][source]) * voxel_measure_;
                if(!external_r_.empty()) row+=(external_r_[source]+external_K_[source])*voxel_measure_;
                prefix[static_cast<std::size_t>(y) * pitch + x] = row;
            }
        }
        for (int x = 1; x <= width; ++x) {
            for (int y = 1; y <= height; ++y) {
                prefix[static_cast<std::size_t>(y) * pitch + x] +=
                    prefix[static_cast<std::size_t>(y - 1) * pitch + x];
            }
        }
        const auto sum = [&](int x0, int y0, int x1, int y1) {
            x0 -= source_x0;
            y0 -= source_y0;
            x1 -= source_x0;
            y1 -= source_y0;
            return prefix[static_cast<std::size_t>(y1) * pitch + x1]
                - prefix[static_cast<std::size_t>(y0) * pitch + x1]
                - prefix[static_cast<std::size_t>(y1) * pitch + x0]
                + prefix[static_cast<std::size_t>(y0) * pitch + x0];
        };
        for (int qy = first_qy; qy <= last_qy; ++qy) {
            for (int qx = first_qx; qx <= last_qx; ++qx) {
                const int center_x = qx * query + query / 2 - origin_x;
                const int center_y = qy * query + query / 2 - origin_y;
                const int x0 = std::clamp(center_x - lower, 0, nx);
                const int y0 = std::clamp(center_y - lower, 0, ny);
                const int x1 = std::clamp(center_x + upper + 1, 0, nx);
                const int y1 = std::clamp(center_y + upper + 1, 0, ny);
                const double count = x1 > x0 && y1 > y0
                    ? sum(x0, y0, x1, y1) : 0.0;
                const double small_density = count / static_cast<double>(edge * edge);
                const double large_density = count /
                    (static_cast<double>(edge * edge) / large_cell_volume_);
                const int bx0 = std::clamp(qx * query - origin_x, 0, nx);
                const int by0 = std::clamp(qy * query - origin_y, 0, ny);
                const int bx1 = std::clamp((qx + 1) * query - origin_x, 0, nx);
                const int by1 = std::clamp((qy + 1) * query - origin_y, 0, ny);
                if (bx1 <= bx0 || by1 <= by0) continue;
                if (!activation_density_bounds_.valid) {
                    activation_density_bounds_ =
                        {bx0, by0, 0, bx1, by1, 1, true};
                } else {
                    activation_density_bounds_.x0 = std::min(
                        activation_density_bounds_.x0, bx0);
                    activation_density_bounds_.y0 = std::min(
                        activation_density_bounds_.y0, by0);
                    activation_density_bounds_.x1 = std::max(
                        activation_density_bounds_.x1, bx1);
                    activation_density_bounds_.y1 = std::max(
                        activation_density_bounds_.y1, by1);
                }
                for (int y = by0; y < by1; ++y) {
                    for (int x = bx0; x < bx1; ++x) {
                        const std::size_t here = index(x, y, 0);
                        activation_density_[0][here] = small_density;
                        activation_density_[1][here] = large_density;
                    }
                }
            }
        }
        return;
    }

    activation_density_bounds_ = {0, 0, 0, nx, ny, nz, true};

    const int px = nx + 1;
    const int py = ny + 1;
    const int pz = nz + 1;
    const auto pindex = [px, py](int x, int y, int z) {
        return (static_cast<std::size_t>(z) * py + y) * px + x;
    };
    std::vector<double> prefix(static_cast<std::size_t>(px) * py * pz, 0.0);
    for (int z = 1; z <= nz; ++z) {
        for (int y = 1; y <= ny; ++y) {
            for (int x = 1; x <= nx; ++x) {
                const std::size_t source = index(x - 1, y - 1, z - 1);
                double value = (r_normal_[0][source] + r_normal_[1][source] +
                    r_active(StructuredStage3D::small, source) +
                    r_active(StructuredStage3D::large, source) +
                    K_[0][source] + K_[1][source]) * voxel_measure_;
                if(!external_r_.empty()) value+=(external_r_[source]+external_K_[source])*voxel_measure_;
                prefix[pindex(x, y, z)] = value
                    + prefix[pindex(x - 1, y, z)]
                    + prefix[pindex(x, y - 1, z)]
                    + prefix[pindex(x, y, z - 1)]
                    - prefix[pindex(x - 1, y - 1, z)]
                    - prefix[pindex(x - 1, y, z - 1)]
                    - prefix[pindex(x, y - 1, z - 1)]
                    + prefix[pindex(x - 1, y - 1, z - 1)];
            }
        }
    }
    const auto sum = [&](int x0, int y0, int z0, int x1, int y1, int z1) {
        return prefix[pindex(x1, y1, z1)]
            - prefix[pindex(x0, y1, z1)] - prefix[pindex(x1, y0, z1)]
            - prefix[pindex(x1, y1, z0)] + prefix[pindex(x0, y0, z1)]
            + prefix[pindex(x0, y1, z0)] + prefix[pindex(x1, y0, z0)]
            - prefix[pindex(x0, y0, z0)];
    };
    const int first_qx = floor_div(origin_x, query);
    const int last_qx = floor_div(origin_x + nx - 1, query);
    const int first_qy = floor_div(origin_y, query);
    const int last_qy = floor_div(origin_y + ny - 1, query);
    const int first_qz = floor_div(origin_z, query);
    const int last_qz = floor_div(origin_z + nz - 1, query);
    for (int qz = first_qz; qz <= last_qz; ++qz) {
        for (int qy = first_qy; qy <= last_qy; ++qy) {
            for (int qx = first_qx; qx <= last_qx; ++qx) {
                const int cx = qx * query + query / 2 - origin_x;
                const int cy = qy * query + query / 2 - origin_y;
                const int cz = qz * query + query / 2 - origin_z;
                const int x0 = std::clamp(cx - lower, 0, nx);
                const int y0 = std::clamp(cy - lower, 0, ny);
                const int z0 = std::clamp(cz - lower, 0, nz);
                const int x1 = std::clamp(cx + upper + 1, 0, nx);
                const int y1 = std::clamp(cy + upper + 1, 0, ny);
                const int z1 = std::clamp(cz + upper + 1, 0, nz);
                const double count = x1 > x0 && y1 > y0 && z1 > z0
                    ? sum(x0, y0, z0, x1, y1, z1) : 0.0;
                const double capacity = static_cast<double>(edge) * edge * edge;
                const double small_density = count / capacity;
                const double large_density = count / (capacity / large_cell_volume_);
                const int bx0 = std::clamp(qx * query - origin_x, 0, nx);
                const int by0 = std::clamp(qy * query - origin_y, 0, ny);
                const int bz0 = std::clamp(qz * query - origin_z, 0, nz);
                const int bx1 = std::clamp((qx + 1) * query - origin_x, 0, nx);
                const int by1 = std::clamp((qy + 1) * query - origin_y, 0, ny);
                const int bz1 = std::clamp((qz + 1) * query - origin_z, 0, nz);
                for (int z = bz0; z < bz1; ++z) {
                    for (int y = by0; y < by1; ++y) {
                        for (int x = bx0; x < bx1; ++x) {
                            const std::size_t here = index(x, y, z);
                            activation_density_[0][here] = small_density;
                            activation_density_[1][here] = large_density;
                        }
                    }
                }
            }
        }
    }
}

void StructuredPdeModel3D::refresh_activation(double dt) {
    build_activation_density();
    if (!population_bounds_.valid) return;
    const double threshold = config_.continuum.base.migration_activation_threshold;
    const double r_inherent = mean_growth_rate(CellType::r);
    const double full_cycle =
        config_.continuum.base.division_timing.base_cycle_hours /
        std::max(1.0e-12, r_inherent);
    const double mean_fraction =
        config_.continuum.base.migration_activation_duration_mean_fraction;
    const double mean_duration = mean_fraction * full_cycle;
    const auto bounds = population_bounds_;
    for (std::size_t stage = 0; stage < kStructuredStageCount3D; ++stage) {
        for (int z = bounds.z0; z < bounds.z1; ++z) {
            for (int y = bounds.y0; y < bounds.y1; ++y) {
                for (int x = bounds.x0; x < bounds.x1; ++x) {
                    const std::size_t location = index(x, y, z);
                    double mass = r_normal_[stage][location];
                    if (operators_.activation.transported_refractory) {
                        auto& refractory = r_refractory_[stage][location];
                        auto& clock = refractory_clock_[stage][location];
                        clock = std::max(0.0, clock - refractory * dt);
                        if (clock <= 1.0e-12 * std::max(1.0, refractory) &&
                            activation_density_[stage][location] <=
                                config_.migration.reactivation_density_threshold) {
                            refractory = 0.0;
                            clock = 0.0;
                        }
                        mass = std::max(0.0, mass - refractory);
                    } else if (operators_.nutrient.transient_resources) {
                        auto& cooldown = activation_cooldown_[stage][location];
                        auto& armed = activation_armed_[stage][location];
                        if (active_total_[stage][location] >=
                            config_.migration.minimum_density) {
                            armed = 0U;
                            cooldown = static_cast<float>(
                                config_.migration.reactivation_cooldown_hours);
                        } else {
                            cooldown = static_cast<float>(std::max(
                                0.0, static_cast<double>(cooldown) - dt));
                            if (cooldown <= 0.0F &&
                                activation_density_[stage][location] <=
                                    config_.migration
                                        .reactivation_density_threshold) {
                                armed = 1U;
                            }
                        }
                        if (armed == 0U) continue;
                    }
                    if (mass < config_.migration.minimum_density ||
                        activation_density_[stage][location] < threshold) continue;
                    r_normal_[stage][location] -= mass;
                    if (duration_) duration_->add_fresh(location, stage, mass);
                    if (velocity_) velocity_->add_fresh(location, stage, mass);
                    active_direction_[stage][0][location] +=
                        static_cast<float>(mass);
                    active_clock_[stage][0][location] +=
                        static_cast<float>(mass * mean_duration);
                    active_total_[stage][location] += static_cast<float>(mass);
                    if (operators_.nutrient.transient_resources) {
                        activation_armed_[stage][location] = 0U;
                        activation_cooldown_[stage][location] =
                            static_cast<float>(
                                config_.migration.reactivation_cooldown_hours);
                    }
                    include_active_location(stage, x, y, z);
                }
            }
        }
    }
}

void StructuredPdeModel3D::expire_active(std::size_t stage, double dt) {
    const double epsilon = config_.migration.minimum_density;
    const auto bounds = active_bounds_[stage];
    if (!bounds.valid) return;
    const auto expire_location = [&](std::size_t location) {
        if (duration_) {
            double before = 0.0;
            for (const auto& field : active_direction_[stage]) before += field[location];
            if (!(before > 0.0)) return;
            duration_->advance(location, stage, dt);
            const double remaining = duration_->mass(location, stage);
            const double factor = before > 0.0 ? std::clamp(remaining / before, 0.0, 1.0) : 0.0;
            const double mean = duration_->mean_work(location, stage);
            double after = 0.0;
            for (std::size_t bucket = 0; bucket < active_direction_[stage].size(); ++bucket) {
                auto& mass = active_direction_[stage][bucket][location];
                mass = static_cast<float>(mass * factor);
                active_clock_[stage][bucket][location] = static_cast<float>(mass * mean);
                after += mass;
            }
            const double expired = std::max(0.0, before - after);
            r_normal_[stage][location] += expired;
            r_refractory_[stage][location] += expired;
            refractory_clock_[stage][location] += expired * config_.migration.reactivation_cooldown_hours;
            active_total_[stage][location] = static_cast<float>(after);
            duration_->reconcile(location, stage, after);
            if (velocity_) velocity_->reconcile(location, stage, after);
            return;
        }
        double total = 0.0;
        for (std::size_t bucket = 0;
             bucket < active_direction_[stage].size(); ++bucket) {
            auto& mass = active_direction_[stage][bucket][location];
            auto& clock_value = active_clock_[stage][bucket][location];
            const double mass_value = mass;
            const double clock = std::max(
                0.0, static_cast<double>(clock_value));
            if (mass_value <= epsilon ||
                clock <= mass_value * dt +
                    1.0e-6 * std::max(1.0, clock)) {
                if (mass_value > 0.0) {
                    r_normal_[stage][location] += mass_value;
                    if (operators_.activation.transported_refractory) {
                        r_refractory_[stage][location] += mass_value;
                        refractory_clock_[stage][location] += mass_value *
                            config_.migration.reactivation_cooldown_hours;
                    } else if (operators_.nutrient.transient_resources) {
                        activation_armed_[stage][location] = 0U;
                        activation_cooldown_[stage][location] =
                            static_cast<float>(std::max(
                                static_cast<double>(
                                    activation_cooldown_[stage][location]),
                                config_.migration
                                    .reactivation_cooldown_hours));
                    }
                }
                mass = 0.0F;
                clock_value = 0.0F;
            } else {
                clock_value = static_cast<float>(
                    clock - mass_value * dt);
                total += mass_value;
            }
        }
        active_total_[stage][location] = static_cast<float>(total);
    };

    if (operators_.transport.cached_cohort_transport) {
        const std::size_t bx = static_cast<std::size_t>(
            bounds.x1 - bounds.x0);
        const std::size_t by = static_cast<std::size_t>(
            bounds.y1 - bounds.y0);
        const std::size_t bz = static_cast<std::size_t>(
            bounds.z1 - bounds.z0);
        const std::size_t box_size = bx * by * bz;
        const int workers = duration_ ? 1 : std::max(1, std::min(
            config_.continuum.base.threads, available_worker_threads()));
        deterministic_parallel_for(
            box_size, workers, [&](std::size_t offset) {
                const int x = bounds.x0 + static_cast<int>(offset % bx);
                const std::size_t yz = offset / bx;
                const int y = bounds.y0 + static_cast<int>(yz % by);
                const int z = bounds.z0 + static_cast<int>(yz / by);
                expire_location(index(x, y, z));
            });
        return;
    }

    for (int z = bounds.z0; z < bounds.z1; ++z) {
        for (int y = bounds.y0; y < bounds.y1; ++y) {
            for (int x = bounds.x0; x < bounds.x1; ++x) {
                expire_location(index(x, y, z));
            }
        }
    }
}

}  // namespace atcg3d::structured_pde
