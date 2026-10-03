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

double StructuredPdeModel3D::normal_diffusion(
    StructuredStage3D stage, CellType type) const noexcept {
    const auto& base = config_.continuum.base;
    const double rate = type == CellType::r
        ? beta_mean(base.normal_r_migration_beta)
        : (base.initial_K_migration_rate_model == "fixed"
               ? base.initial_K_migration_rate
               : beta_mean(base.initial_K_migration_beta));
    double result = (base.thin_layer &&
        config_.continuum.migration.mapping == "shared_fixed_lattice_means_v2"
        ? 3.0 / 8.0 : kFixed26DiffusionFactor) * rate *
        config_.continuum.migration.diffusion_scale;
    if (stage == StructuredStage3D::large) {
        result *= config_.continuum.migration.large_mobility_multiplier;
    }
    return result;
}

double StructuredPdeModel3D::active_rate(StructuredStage3D stage) const noexcept {
    double result = beta_mean(config_.continuum.base.normal_r_migration_beta) *
        config_.continuum.base.activated_r_normal_multiplier;
    if (stage == StructuredStage3D::large) {
        result *= config_.continuum.migration.large_mobility_multiplier;
    }
    return result;
}

void StructuredPdeModel3D::migrate_normal_and_K(double dt) {
    const auto& continuum = config_.continuum;
    if (!population_bounds_.valid) return;
    if (renewal_) renewal_->begin_transport();
    const int nx = continuum.grid.shape[0];
    const int ny = continuum.grid.shape[1];
    const int nz = continuum.grid.shape[2];
    const double inverse_h2 = 1.0 /
        (continuum.grid.spacing_voxels * continuum.grid.spacing_voxels);
    const int workers = renewal_ ? 1 : std::max(
        1, std::min(continuum.base.threads, available_worker_threads()));
    const StructuredActiveBounds3D old_bounds = population_bounds_;
    const StructuredActiveBounds3D target_bounds = expanded_bounds(old_bounds);
    const auto clear = [&](const StructuredActiveBounds3D& bounds) {
        if (!bounds.valid) return;
        for (int z = bounds.z0; z < bounds.z1; ++z) {
            for (int y = bounds.y0; y < bounds.y1; ++y) {
                const std::size_t begin = index(bounds.x0, y, z);
                const std::size_t end = index(bounds.x1 - 1, y, z) + 1;
                for (std::size_t stage = 0; stage < 2; ++stage) {
                    std::fill(r_normal_work_[stage].begin() + begin,
                              r_normal_work_[stage].begin() + end, 0.0);
                    std::fill(K_work_[stage].begin() + begin,
                              K_work_[stage].begin() + end, 0.0);
                    if (operators_.activation.transported_refractory) {
                        std::fill(refractory_work_[stage].begin() + begin,
                                  refractory_work_[stage].begin() + end, 0.0);
                        std::fill(refractory_clock_work_[stage].begin() + begin,
                                  refractory_clock_work_[stage].begin() + end, 0.0);
                    }
                }
            }
        }
    };
    clear(normal_work_dirty_bounds_);
    clear(target_bounds);
    const std::size_t bx = static_cast<std::size_t>(
        target_bounds.x1 - target_bounds.x0);
    const std::size_t by = static_cast<std::size_t>(
        target_bounds.y1 - target_bounds.y0);
    const std::size_t bz = static_cast<std::size_t>(
        target_bounds.z1 - target_bounds.z0);
    const std::size_t box_size = bx * by * bz;

    std::array<std::vector<std::array<double, 27>>, 2> feasible_probabilities;
    const bool feasible_jump = operators_.transport.feasible_normal_jumps;
    const auto box_offset = [&](int x, int y, int z) {
        return (static_cast<std::size_t>(z - target_bounds.z0) * by +
            static_cast<std::size_t>(y - target_bounds.y0)) * bx +
            static_cast<std::size_t>(x - target_bounds.x0);
    };
    if (feasible_jump) {
        for (std::size_t stage = 0; stage < kStructuredStageCount3D; ++stage) {
            auto& cache = feasible_probabilities[stage];
            cache.resize(box_size);
            for (int z = old_bounds.z0; z < old_bounds.z1; ++z) {
                for (int y = old_bounds.y0; y < old_bounds.y1; ++y) {
                    for (int x = old_bounds.x0; x < old_bounds.x1; ++x) {
                        const auto here = index(x, y, z);
                        if (vessel_blocks_cells(here)) continue;
                        if (r_normal_[stage][here] == 0.0 && K_[stage][here] == 0.0 &&
                            r_refractory_[stage][here] == 0.0 && refractory_clock_[stage][here] == 0.0) {
                            continue;
                        }
                        std::vector<double> availability(direction_ids_.size(), 0.0);
                        for (std::size_t d = 0; d < direction_ids_.size(); ++d) {
                            const auto jump = direction_vector(direction_ids_[d]);
                            const int ox = x + jump.x;
                            const int oy = y + jump.y;
                            const int oz = z + jump.z;
                            if (ox < 0 || ox >= nx || oy < 0 || oy >= ny || oz < 0 || oz >= nz) continue;
                            const auto other = index(ox, oy, oz);
                            if (transport_destination_blocks_cells(
                                    other, static_cast<StructuredStage3D>(stage))) continue;
                            const double vacancy = 1.0 - std::clamp(occupied_fraction(other) /
                                continuum.reaction.maximum_occupied_fraction, 0.0, 1.0);
                            const double overlapping_sites = (2 - std::abs(jump.x)) *
                                (2 - std::abs(jump.y)) *
                                (continuum.base.thin_layer ? 1 : 2 - std::abs(jump.z));
                            const double new_sites = stage == 0 ? 1.0 : large_cell_volume_ - overlapping_sites;
                            availability[d] = std::pow(vacancy, new_sites);
                        }
                        const auto probabilities = uniform_feasible_jump_probabilities(availability);
                        auto& target = cache[box_offset(x, y, z)];
                        for (std::size_t d = 0; d < direction_ids_.size(); ++d) {
                            target[direction_ids_[d]] = probabilities[d];
                        }
                    }
                }
            }
        }
    }

    const auto migrate_field = [&](const PagedField<double>& source,
                                   PagedField<double>& target,
                                   StructuredStage3D stage,
                                   CellType type, bool carry_cycle = false) {
        const bool large = stage == StructuredStage3D::large;
        const double vacancy_exponent = large ? large_cell_volume_ : 1.0;
        const double diffusion = normal_diffusion(stage, type);
        const bool lattice_jump = operators_.transport.fixed_normal_jumps;
        const double jump_coefficient = lattice_jump
            ? diffusion / (continuum.base.thin_layer ? 3.0 : 9.0) : diffusion;
        std::atomic<bool> invalid_negative{false};
        deterministic_parallel_for(box_size, workers, [&](std::size_t offset) {
            const int x = target_bounds.x0 + static_cast<int>(offset % bx);
            const std::size_t yz = offset / bx;
            const int y = target_bounds.y0 + static_cast<int>(yz % by);
            const int z = target_bounds.z0 + static_cast<int>(yz / by);
            const std::size_t here = index(x, y, z);
            const double local_occupied = std::clamp(
                occupied_fraction(here) /
                    continuum.reaction.maximum_occupied_fraction,
                0.0, 1.0);
            const double vacancy = 1.0 - local_occupied;
            const double availability = transport_destination_blocks_cells(here, stage)
                ? 0.0 : std::pow(vacancy, vacancy_exponent);
            double delta = 0.0;
            const auto exchange = [&](int ox, int oy, int oz, DirectionId direction = kStayDirection) {
                if (ox < 0 || ox >= nx || oy < 0 || oy >= ny ||
                    oz < 0 || oz >= nz) return;
                const std::size_t other = index(ox, oy, oz);
                const double other_occupied = std::clamp(
                    occupied_fraction(other) /
                        continuum.reaction.maximum_occupied_fraction,
                    0.0, 1.0);
                const double other_vacancy = 1.0 - other_occupied;
                const double other_availability = transport_destination_blocks_cells(other, stage)
                    ? 0.0 : std::pow(other_vacancy, vacancy_exponent);
                const double mobility = std::pow(
                    0.5 * (vacancy + other_vacancy),
                    continuum.migration.crowding_exponent);
                double incoming = jump_coefficient * inverse_h2 * mobility * availability;
                double outgoing = jump_coefficient * inverse_h2 * mobility * other_availability;
                if (feasible_jump) {
                    const double rate = diffusion / (continuum.base.thin_layer ? 3.0 / 8.0 : kFixed26DiffusionFactor);
                    const auto& cache = feasible_probabilities[std::size_t(stage)];
                    const bool in_box = ox >= target_bounds.x0 && ox < target_bounds.x1 &&
                        oy >= target_bounds.y0 && oy < target_bounds.y1 &&
                        oz >= target_bounds.z0 && oz < target_bounds.z1;
                    incoming = in_box ? rate * cache[box_offset(ox, oy, oz)][opposite_direction(direction)] : 0.0;
                    outgoing = rate * cache[box_offset(x, y, z)][direction];
                }
                if (feasible_jump) {
                    delta += source[other] * incoming - source[here] * outgoing;
                } else {
                    delta += jump_coefficient * inverse_h2 * mobility *
                        (source[other] * availability - source[here] * other_availability);
                }
                if (renewal_ && carry_cycle && here < other) {
                    const std::size_t channel = (type == CellType::r ? 0 : 2) + std::size_t(stage);
                    if (feasible_jump) {
                        renewal_->transfer(other, here, channel, dt * source[other] * incoming);
                        renewal_->transfer(here, other, channel, dt * source[here] * outgoing);
                    } else {
                        const double scale = dt * jump_coefficient * inverse_h2 * mobility;
                        renewal_->transfer(other, here, channel, scale * source[other] * availability);
                        renewal_->transfer(here, other, channel, scale * source[here] * other_availability);
                    }
                }
            };
            if (lattice_jump) {
                for (const auto direction_id : direction_ids_) {
                    const auto direction = direction_vector(direction_id);
                    exchange(x + direction.x, y + direction.y, z + direction.z, direction_id);
                }
            } else {
                exchange(x - 1, y, z);
                exchange(x + 1, y, z);
                exchange(x, y - 1, z);
                exchange(x, y + 1, z);
                if (!continuum.base.thin_layer) {
                    exchange(x, y, z - 1);
                    exchange(x, y, z + 1);
                }
            }
            const double value = source[here] + dt * delta;
            if (value < -1.0e-10) invalid_negative.store(true);
            target[here] = value >= config_.migration.minimum_density
                ? std::max(0.0, value) : 0.0;
        });
        if (invalid_negative.load()) {
            throw std::runtime_error(
                "structured normal migration produced negative density");
        }
    };

    for (std::size_t stage = 0; stage < kStructuredStageCount3D; ++stage) {
        const auto value = static_cast<StructuredStage3D>(stage);
        migrate_field(r_normal_[stage], r_normal_work_[stage], value, CellType::r, true);
        migrate_field(K_[stage], K_work_[stage], value, CellType::K, true);
        if (operators_.activation.transported_refractory) {
            migrate_field(r_refractory_[stage], refractory_work_[stage], value, CellType::r);
            migrate_field(refractory_clock_[stage], refractory_clock_work_[stage], value, CellType::r);
        }
    }
    r_normal_.swap(r_normal_work_);
    K_.swap(K_work_);
    if (operators_.activation.transported_refractory) {
        r_refractory_.swap(refractory_work_);
        refractory_clock_.swap(refractory_clock_work_);
        for (std::size_t stage = 0; stage < 2; ++stage) {
            for (int z = target_bounds.z0; z < target_bounds.z1; ++z) {
                for (int y = target_bounds.y0; y < target_bounds.y1; ++y) {
                    for (int x = target_bounds.x0; x < target_bounds.x1; ++x) {
                        const auto here = index(x, y, z);
                        const double mass = r_refractory_[stage][here];
                        if (mass > r_normal_[stage][here]) {
                            refractory_clock_[stage][here] *= r_normal_[stage][here] / mass;
                            r_refractory_[stage][here] = r_normal_[stage][here];
                        }
                        if (r_refractory_[stage][here] == 0.0) refractory_clock_[stage][here] = 0.0;
                    }
                }
            }
        }
    }
    normal_work_dirty_bounds_ = old_bounds;
    population_bounds_ = target_bounds;
    finish_division_transport();
}

std::vector<std::size_t> StructuredPdeModel3D::eligible_initial_directions(
    std::size_t location) const {
    const auto& continuum = config_.continuum;
    const int nx = continuum.grid.shape[0];
    const int ny = continuum.grid.shape[1];
    const int nz = continuum.grid.shape[2];
    const int x = static_cast<int>(location % static_cast<std::size_t>(nx));
    const std::size_t yz = location / static_cast<std::size_t>(nx);
    const int y = static_cast<int>(yz % static_cast<std::size_t>(ny));
    const int z = static_cast<int>(yz / static_cast<std::size_t>(ny));
    const int radius = continuum.base.direction_density_radius;
    const double minimum_cosine = std::cos(
        continuum.base.direction_density_half_angle_degrees *
        std::acos(-1.0) / 180.0);
    std::vector<std::size_t> result;
    for (std::size_t direction_index = 0;
         direction_index < direction_ids_.size(); ++direction_index) {
        const Vec3i forward = direction_vector(direction_ids_[direction_index]);
        const int target_x = x + forward.x;
        const int target_y = y + forward.y;
        const int target_z = z + forward.z;
        if (target_x < 0 || target_x >= nx || target_y < 0 || target_y >= ny ||
            target_z < 0 || target_z >= nz ||
            vessel_blocks_cells(index(target_x, target_y, target_z))) {
            continue;
        }
        const double forward_length =
            std::sqrt(static_cast<double>(squared_length(forward)));
        double count = 0.0;
        std::size_t sites = 0;
        for (int dz = continuum.base.thin_layer ? 0 : -radius;
             dz <= (continuum.base.thin_layer ? 0 : radius); ++dz) {
            for (int dy = -radius; dy <= radius; ++dy) {
                for (int dx = -radius; dx <= radius; ++dx) {
                    const int distance = std::max(
                        {std::abs(dx), std::abs(dy), std::abs(dz)});
                    if (distance == 0 || distance > radius) continue;
                    const Vec3i offset{dx, dy, dz};
                    const double offset_length =
                        std::sqrt(static_cast<double>(squared_length(offset)));
                    const double cosine = static_cast<double>(dot(offset, forward)) /
                        (offset_length * forward_length);
                    if (cosine + 1.0e-12 < minimum_cosine) continue;
                    ++sites;
                    const int ox = x + dx;
                    const int oy = y + dy;
                    const int oz = z + dz;
                    if (ox < 0 || ox >= nx || oy < 0 || oy >= ny ||
                        oz < 0 || oz >= nz) continue;
                    const std::size_t other = index(ox, oy, oz);
                    count += (r_normal_[0][other] + r_normal_[1][other] +
                              r_active(StructuredStage3D::small, other) +
                              r_active(StructuredStage3D::large, other) +
                              K_[0][other] + K_[1][other]) * voxel_measure_;
                }
            }
        }
        const double density = sites > 0 ? count / static_cast<double>(sites) : 1.0;
        if (density <= continuum.base.direction_density_threshold) {
            result.push_back(direction_index + 1);
        }
    }
    return result;
}

void StructuredPdeModel3D::build_guidance_prefix(double dt) {
    guidance_prefix_bounds_ = {};
    guidance_prefix_pitch_ = 0;
    if (sector_mean_) {
        std::vector<continuum::SectorQueryBox3D> boxes;
        const auto& continuum = config_.continuum;
        for (std::size_t stage = 0; stage < kStructuredStageCount3D; ++stage) {
            const auto& bounds = active_bounds_[stage];
            if (!bounds.valid) continue;
            const double rate = velocity_
                ? continuum.base.normal_r_migration_beta.scale *
                    continuum.migration.activated_r_mobility_multiplier *
                    (stage == 1 ? continuum.migration.large_mobility_multiplier : 1.0)
                : active_rate(static_cast<StructuredStage3D>(stage));
            const int margin = static_cast<int>(std::ceil(rate * dt /
                -std::log1p(-config_.migration.maximum_move_probability_per_substep))) + 1;
            continuum::SectorQueryBox3D box{
                {bounds.x0 - margin, bounds.y0 - margin, bounds.z0 - margin},
                {bounds.x1 + margin, bounds.y1 + margin, bounds.z1 + margin}};
            if (continuum.base.thin_layer) {
                box.lower[2] = 0;
                box.upper[2] = 1;
            }
            boxes.push_back(box);
        }
        sector_mean_->prepare(nutrient_, continuum.nutrient.vessel_value,
            boxes, step_count_ + 1);
        return;
    }
    if (!operators_.transport.directional_sectors ||
        !config_.continuum.base.thin_layer) return;

    StructuredActiveBounds3D bounds;
    for (const auto& stage_bounds : active_bounds_) {
        if (!stage_bounds.valid) continue;
        if (!bounds.valid) {
            bounds = stage_bounds;
        } else {
            bounds.x0 = std::min(bounds.x0, stage_bounds.x0);
            bounds.y0 = std::min(bounds.y0, stage_bounds.y0);
            bounds.x1 = std::max(bounds.x1, stage_bounds.x1);
            bounds.y1 = std::max(bounds.y1, stage_bounds.y1);
        }
    }
    if (!bounds.valid) return;

    int maximum_substeps = 1;
    for (std::size_t stage = 0; stage < kStructuredStageCount3D; ++stage) {
        const double rate = active_rate(static_cast<StructuredStage3D>(stage));
        if (!(rate > 0.0)) continue;
        maximum_substeps = std::max(maximum_substeps,
            static_cast<int>(std::ceil(
                rate * dt /
                -std::log1p(-config_.migration.maximum_move_probability_per_substep))));
    }
    const int edge = operators_.nutrient.transient_resources
        ? config_.migration.direction_nutrient_window_edge
        : config_.migration.direction_density_window_edge;
    const int lower = (edge - 1) / 2;
    const int upper = edge - lower - 1;
    const int margin = std::max(lower, upper) + maximum_substeps + 1;
    const int nx = config_.continuum.grid.shape[0];
    const int ny = config_.continuum.grid.shape[1];
    bounds.x0 = std::max(0, bounds.x0 - margin);
    bounds.y0 = std::max(0, bounds.y0 - margin);
    bounds.x1 = std::min(nx, bounds.x1 + margin);
    bounds.y1 = std::min(ny, bounds.y1 + margin);
    bounds.z0 = 0;
    bounds.z1 = 1;
    guidance_prefix_bounds_ = bounds;
    guidance_prefix_pitch_ = bounds.x1 - bounds.x0 + 1;
    const int height = bounds.y1 - bounds.y0;
    const std::size_t size = static_cast<std::size_t>(height) *
        static_cast<std::size_t>(guidance_prefix_pitch_);
    guidance_density_row_prefix_.assign(size, 0.0);
    guidance_resource_row_prefix_.assign(size, 0.0);
    const double vessel_value = config_.continuum.nutrient.vessel_value;
    for (int y = bounds.y0; y < bounds.y1; ++y) {
        const std::size_t row = static_cast<std::size_t>(y - bounds.y0) *
            static_cast<std::size_t>(guidance_prefix_pitch_);
        double density_sum = 0.0;
        double resource_sum = 0.0;
        for (int x = bounds.x0; x < bounds.x1; ++x) {
            const std::size_t location = index(x, y, 0);
            if (!operators_.nutrient.transient_resources) {
                density_sum += (r_normal_[0][location] +
                    r_normal_[1][location] +
                    r_active(StructuredStage3D::small, location) +
                    r_active(StructuredStage3D::large, location) +
                    K_[0][location] + K_[1][location]) * voxel_measure_;
            }
            resource_sum += std::clamp(
                nutrient_[location] / vessel_value, 0.0, 1.0);
            const std::size_t column =
                static_cast<std::size_t>(x - bounds.x0 + 1);
            guidance_density_row_prefix_[row + column] = density_sum;
            guidance_resource_row_prefix_[row + column] = resource_sum;
        }
    }
}

std::vector<double> StructuredPdeModel3D::guided_direction_weights(
    std::size_t location,
    bool allow_crowded_sectors) const {
    const auto& continuum = config_.continuum;
    const auto& base = continuum.base;
    const int nx = continuum.grid.shape[0];
    const int ny = continuum.grid.shape[1];
    const int nz = continuum.grid.shape[2];
    const int x = static_cast<int>(location % static_cast<std::size_t>(nx));
    const std::size_t yz = location / static_cast<std::size_t>(nx);
    const int y = static_cast<int>(yz % static_cast<std::size_t>(ny));
    const int z = static_cast<int>(yz / static_cast<std::size_t>(ny));
    const int edge = operators_.nutrient.transient_resources
        ? config_.migration.direction_nutrient_window_edge
        : operators_.transport.directional_sectors
        ? config_.migration.direction_density_window_edge
        : 2 * base.direction_density_radius + 1;
    const int lower = (edge - 1) / 2;
    const int upper = edge - lower - 1;
    const double minimum_cosine = std::cos(
        base.direction_density_half_angle_degrees * std::acos(-1.0) / 180.0);
    const double floor = base.direction_minimum_guidance_weight;
    std::vector<double> result(direction_ids_.size() + 1, 0.0);
    for (std::size_t direction_index = 0;
         direction_index < direction_ids_.size(); ++direction_index) {
        const Vec3i forward = direction_vector(direction_ids_[direction_index]);
        const int target_x = x + forward.x;
        const int target_y = y + forward.y;
        const int target_z = z + forward.z;
        if (target_x < 0 || target_x >= nx || target_y < 0 || target_y >= ny ||
            target_z < 0 || target_z >= nz ||
            vessel_blocks_cells(index(target_x, target_y, target_z))) {
            continue;
        }
        const double forward_length =
            std::sqrt(static_cast<double>(squared_length(forward)));
        if (sector_mean_) {
            if (sector_mean_->count(location, direction_ids_[direction_index]) == 0) continue;
            double gradient = (sector_mean_->mean(location, direction_ids_[direction_index]) -
                std::clamp(nutrient_[location] / continuum.nutrient.vessel_value, 0.0, 1.0)) *
                continuum.nutrient.vessel_value;
            if (std::abs(gradient) < config_.migration.zero_gradient_tolerance)
                gradient = 0.0;
            result[direction_index + 1] =
                std::pow(forward_length, -base.distance_weight_exponent) *
                std::exp(std::clamp(config_.migration.chemotaxis_strength * gradient, -40.0, 40.0));
            continue;
        }
        double count = 0.0;
        double resource = 0.0;
        std::size_t sites = 0;
        std::size_t resource_sites = 0;
        if (operators_.transport.directional_sectors && base.thin_layer &&
            guidance_prefix_bounds_.valid) {
            if (!operators_.activation.transported_refractory) {
                sites = direction_sector_site_counts_[direction_index];
            }
            for (const auto& span : direction_row_spans_[direction_index]) {
                const int oy = y + span.dy;
                if (oy < guidance_prefix_bounds_.y0 ||
                    oy >= guidance_prefix_bounds_.y1) continue;
                const int ox0 = std::max(
                    guidance_prefix_bounds_.x0, x + span.dx0);
                const int ox1 = std::min(
                    guidance_prefix_bounds_.x1, x + span.dx1 + 1);
                if (ox1 <= ox0) continue;
                const std::size_t row =
                    static_cast<std::size_t>(oy - guidance_prefix_bounds_.y0) *
                    static_cast<std::size_t>(guidance_prefix_pitch_);
                const std::size_t left =
                    static_cast<std::size_t>(ox0 - guidance_prefix_bounds_.x0);
                const std::size_t right =
                    static_cast<std::size_t>(ox1 - guidance_prefix_bounds_.x0);
                resource_sites += static_cast<std::size_t>(ox1 - ox0);
                if (operators_.activation.transported_refractory) sites += static_cast<std::size_t>(ox1 - ox0);
                count += guidance_density_row_prefix_[row + right] -
                    guidance_density_row_prefix_[row + left];
                resource += guidance_resource_row_prefix_[row + right] -
                    guidance_resource_row_prefix_[row + left];
            }
        } else {
            for (int dz = base.thin_layer ? 0 : -lower;
                 dz <= (base.thin_layer ? 0 : upper); ++dz) {
                for (int dy = -lower; dy <= upper; ++dy) {
                    for (int dx = -lower; dx <= upper; ++dx) {
                        if (dx == 0 && dy == 0 && dz == 0) continue;
                        const Vec3i offset{dx, dy, dz};
                        const double offset_length =
                            std::sqrt(static_cast<double>(squared_length(offset)));
                        const double cosine =
                            static_cast<double>(dot(offset, forward)) /
                            (offset_length * forward_length);
                        if (cosine + 1.0e-12 < minimum_cosine) continue;
                        if (!operators_.activation.transported_refractory) ++sites;
                        const int ox = x + dx;
                        const int oy = y + dy;
                        const int oz = z + dz;
                        if (ox < 0 || ox >= nx || oy < 0 || oy >= ny ||
                            oz < 0 || oz >= nz) continue;
                        if (operators_.activation.transported_refractory) ++sites;
                        const std::size_t other = index(ox, oy, oz);
                        ++resource_sites;
                        count += (r_normal_[0][other] + r_normal_[1][other] +
                                  r_active(StructuredStage3D::small, other) +
                                  r_active(StructuredStage3D::large, other) +
                                  K_[0][other] + K_[1][other]) * voxel_measure_;
                        resource += std::clamp(
                            nutrient_[other] / continuum.nutrient.vessel_value,
                            0.0, 1.0);
                    }
                }
            }
        }
        if (sites == 0) continue;
        if (operators_.nutrient.transient_resources) {
            if (resource_sites == 0) continue;
            const double local_resource = std::clamp(
                nutrient_[location] / continuum.nutrient.vessel_value,
                0.0, 1.0);
            const double directional_resource =
                resource / static_cast<double>(resource_sites);
            double gradient = directional_resource - local_resource;
            if (base.direction_guidance_model == "nutrient_gradient_shared_resource_v3" ||
                base.direction_guidance_model == "nutrient_gradient_shared_resource_v4") {
                gradient *= continuum.nutrient.vessel_value;
            }
            if (std::abs(gradient) <
                config_.migration.zero_gradient_tolerance) {
                gradient = 0.0;
            }
            const double exponent = std::clamp(
                config_.migration.chemotaxis_strength * gradient,
                -40.0, 40.0);
            result[direction_index + 1] =
                std::pow(forward_length, -base.distance_weight_exponent) *
                std::exp(exponent);
            continue;
        }
        const double density = count / static_cast<double>(sites);
        if (!allow_crowded_sectors &&
            density > base.direction_density_threshold) continue;
        const double density_fraction = std::clamp(
            density / base.direction_density_threshold, 0.0, 1.0);
        const double density_score =
            floor + (1.0 - floor) * (1.0 - density_fraction);
        const double resource_score = floor + (1.0 - floor) *
            resource / static_cast<double>(sites);
        result[direction_index + 1] =
            std::pow(forward_length, -base.distance_weight_exponent) *
            std::pow(density_score,
                     base.direction_density_guidance_exponent) *
            std::pow(resource_score,
                     base.direction_resource_guidance_exponent);
    }
    return result;
}

void StructuredPdeModel3D::migrate_active(double dt) {
    const auto& continuum = config_.continuum;
    const int nx = continuum.grid.shape[0];
    const int ny = continuum.grid.shape[1];
    const int nz = continuum.grid.shape[2];
    const bool conditional_vacancy = config_.migration.direction_transport ==
        "nutrient_gradient_feasible_direction_jump_exchange_v5";
    build_guidance_prefix(dt);
    if (operators_.transport.cached_cohort_transport) {
        ++guidance_weight_cache_generation_;
        if (guidance_weight_cache_generation_ == 0U) {
            fill_field(guidance_weight_cache_stamp_,0U);
            guidance_weight_cache_generation_ = 1U;
        }
    }
    const auto ensure_guidance_cache = [&](std::size_t location) {
        if (!operators_.transport.cached_cohort_transport ||
            guidance_weight_cache_stamp_[location] ==
                guidance_weight_cache_generation_) {
            return;
        }
        const std::vector<double> weights = guided_direction_weights(location);
        for (std::size_t bucket = 0; bucket < weights.size(); ++bucket) {
            guidance_weight_cache_[bucket][location] =
                static_cast<float>(weights[bucket]);
        }
        guidance_weight_cache_stamp_[location] =
            guidance_weight_cache_generation_;
    };

    for (std::size_t stage = 0; stage < kStructuredStageCount3D; ++stage) {
        if (duration_) { expire_active(stage, 0.0); shrink_active_bounds(stage); }
        if (!active_bounds_[stage].valid) continue;
        const auto stage_value = static_cast<StructuredStage3D>(stage);
        const double rate = active_rate(stage_value);
        if (!(rate > 0.0)) {
            if (duration_) { expire_active(stage, dt); shrink_active_bounds(stage); }
            continue;
        }
        const bool resource_guided = operators_.nutrient.transient_resources ||
            (operators_.transport.resource_guidance &&
             continuum.base.direction_guidance_model == "low_density_high_resource_v1");
        // Resolve the no-persistent-direction state before transport. V2/v3
        // rule uses the same low-density/high-resource directional cone score
        // as the nutrient-coupled ABM.
        const auto initial_bounds = active_bounds_[stage];
        for (int z = initial_bounds.z0; z < initial_bounds.z1; ++z) {
            for (int y = initial_bounds.y0; y < initial_bounds.y1; ++y) {
                for (int x = initial_bounds.x0; x < initial_bounds.x1; ++x) {
                    const std::size_t location = index(x, y, z);
                    const double mass = active_direction_[stage][0][location];
                    if (mass < config_.migration.minimum_density) continue;
                    std::vector<std::size_t> targets;
                    std::vector<double> weights(
                        active_direction_[stage].size(), 0.0);
                    double weight_sum = 0.0;
                    if (resource_guided) {
                        if (operators_.transport.cached_cohort_transport) {
                            ensure_guidance_cache(location);
                            for (std::size_t target = 0;
                                 target < weights.size(); ++target) {
                                weights[target] =
                                    guidance_weight_cache_[target][location];
                            }
                        } else {
                            weights = guided_direction_weights(location);
                        }
                        for (std::size_t target = 1;
                             target < weights.size(); ++target) {
                            if (conditional_vacancy) {
                                const auto step = direction_vector(direction_ids_[target - 1]);
                                if (x + step.x < 0 || x + step.x >= nx ||
                                    y + step.y < 0 || y + step.y >= ny ||
                                    z + step.z < 0 || z + step.z >= nz) weights[target] = 0.0;
                                else {
                                    const auto other = index(x + step.x, y + step.y, z + step.z);
                                    weights[target] *= vessel_blocks_cells(other) ? 0.0 :
                                        std::pow(std::clamp(1.0 - occupied_fraction(other) /
                                            continuum.reaction.maximum_occupied_fraction, 0.0, 1.0),
                                            stage == 1 ? large_cell_volume_ : 1.0);
                                }
                            }
                            if (weights[target] > 0.0) {
                                targets.push_back(target);
                                weight_sum += weights[target];
                            }
                        }
                    } else {
                        targets = eligible_initial_directions(location);
                        for (const std::size_t target : targets) {
                            const double length = std::sqrt(static_cast<double>(
                                squared_length(direction_vector(
                                    direction_ids_[target - 1]))));
                            weights[target] = std::pow(
                                length,
                                -continuum.base.distance_weight_exponent);
                            weight_sum += weights[target];
                        }
                    }
                    if (targets.empty() || !(weight_sum > 0.0)) continue;
                    const double clock = active_clock_[stage][0][location];
                    for (const std::size_t target : targets) {
                        const double fraction = weights[target] / weight_sum;
                        active_direction_[stage][target][location] +=
                            static_cast<float>(mass * fraction);
                        active_clock_[stage][target][location] +=
                            static_cast<float>(clock * fraction);
                    }
                    active_direction_[stage][0][location] = 0.0F;
                    active_clock_[stage][0][location] = 0.0F;
                }
            }
        }
        const double maximum_probability =
            config_.migration.maximum_move_probability_per_substep;
        const double rate_supremum = velocity_ ? continuum.base.normal_r_migration_beta.scale *
            continuum.migration.activated_r_mobility_multiplier *
            (stage == 1 ? continuum.migration.large_mobility_multiplier : 1.0) : rate;
        const int substeps = std::max(1, static_cast<int>(std::ceil(
            rate_supremum * dt / -std::log1p(-maximum_probability))));
        const double sub_dt = dt / substeps;
        const double move_probability = 1.0 - std::exp(-rate * sub_dt);
        const double vacancy_exponent = stage == 1 ? large_cell_volume_ : 1.0;
        // Active-r PDE tails occupy most sites inside their bounds, so a
        // contiguous dense traversal is faster than maintaining and sorting a
        // sparse index. Keep the sparse implementation available for future
        // genuinely sparse closures, but use the deterministic dense path for
        // this v4 model.
        constexpr bool sparse_transport = false;

        std::vector<std::size_t> sparse_locations;
        if (sparse_transport) {
            for (int z = initial_bounds.z0; z < initial_bounds.z1; ++z) {
                for (int y = initial_bounds.y0; y < initial_bounds.y1; ++y) {
                    for (int x = initial_bounds.x0; x < initial_bounds.x1; ++x) {
                        const std::size_t location = index(x, y, z);
                        if (active_total_[stage][location] >=
                            config_.migration.minimum_density) {
                            sparse_locations.push_back(location);
                        }
                    }
                }
            }
        }

        for (int substep = 0; substep < substeps; ++substep) {
            if (renewal_) renewal_->begin_transport();
            if (duration_) duration_->begin_transport();
            if (velocity_) {
                std::array<std::vector<double>, 4> probabilities;
                for (int c = 0; c < 2; ++c) {
                    probabilities[c].resize(config_.migration.activation_rate_bins + 1);
                    for (int b = 0; b <= config_.migration.activation_rate_bins; ++b)
                        probabilities[c][b] = -std::expm1(-rate_supremum * b /
                            config_.migration.activation_rate_bins * sub_dt);
                }
                velocity_->begin_weighted_transport(std::move(probabilities));
            }
            const StructuredActiveBounds3D old_bounds = active_bounds_[stage];
            const StructuredActiveBounds3D target_bounds =
                expanded_bounds(old_bounds);
            std::vector<double> proposed_incoming;
            bool collect_incoming = false;
            const auto incoming_offset = [&](int x, int y, int z) {
                return (static_cast<std::size_t>(z - target_bounds.z0) *
                    (target_bounds.y1 - target_bounds.y0) + y - target_bounds.y0) *
                    (target_bounds.x1 - target_bounds.x0) + x - target_bounds.x0;
            };
            if (conditional_vacancy) proposed_incoming.resize(bounds_size(target_bounds));
            StructuredActiveBounds3D clear_bounds = target_bounds;
            if (work_dirty_bounds_.valid) {
                clear_bounds.x0 = std::min(clear_bounds.x0, work_dirty_bounds_.x0);
                clear_bounds.y0 = std::min(clear_bounds.y0, work_dirty_bounds_.y0);
                clear_bounds.z0 = std::min(clear_bounds.z0, work_dirty_bounds_.z0);
                clear_bounds.x1 = std::max(clear_bounds.x1, work_dirty_bounds_.x1);
                clear_bounds.y1 = std::max(clear_bounds.y1, work_dirty_bounds_.y1);
                clear_bounds.z1 = std::max(clear_bounds.z1, work_dirty_bounds_.z1);
            }
            // V4 clears the complete prior dirty region once. After each
            // sparse swap below, only the actual source locations in the old
            // buffer need clearing. Legacy schemas retain their dense path.
            if (!sparse_transport || substep == 0) {
                clear_active_work(clear_bounds);
            }

            std::vector<std::size_t> sparse_targets;
            if (sparse_transport) {
                ++active_location_generation_;
                if (active_location_generation_ == 0U) {
                    fill_field(active_location_stamp_,0U);
                    active_location_generation_ = 1U;
                }
                sparse_targets.reserve(sparse_locations.size() * 2);
            }
            const auto mark_sparse_target = [&](std::size_t location) {
                if (!sparse_transport ||
                    active_location_stamp_[location] ==
                        active_location_generation_) {
                    return;
                }
                active_location_stamp_[location] =
                    active_location_generation_;
                sparse_targets.push_back(location);
            };
            const auto transport_location = [&](int x, int y, int z,
                                                std::size_t location) {
                if (active_total_[stage][location] <
                    config_.migration.minimum_density) {
                    return;
                }
                const double local_probability = velocity_ ? velocity_->weighted_transport_mass(location, stage) /
                    active_total_[stage][location] : move_probability;
                if (resource_guided && operators_.transport.cached_cohort_transport) {
                    ensure_guidance_cache(location);
                }
                const std::vector<double> guidance =
                    resource_guided && !operators_.transport.cached_cohort_transport
                    ? guided_direction_weights(location)
                    : std::vector<double>{};
                std::array<double, 27> vacancies{};
                if (conditional_vacancy) {
                    for (std::size_t b = 1; b < active_direction_[stage].size(); ++b) {
                        const auto step = direction_vector(direction_ids_[b - 1]);
                        if (x + step.x < 0 || x + step.x >= nx ||
                            y + step.y < 0 || y + step.y >= ny ||
                            z + step.z < 0 || z + step.z >= nz) continue;
                        const auto other = index(x + step.x, y + step.y, z + step.z);
                        if (!vessel_blocks_cells(other)) vacancies[b] = std::pow(
                            std::clamp(1.0 - occupied_fraction(other) /
                                continuum.reaction.maximum_occupied_fraction, 0.0, 1.0), vacancy_exponent);
                    }
                }
                const auto raw_guidance_value = [&](std::size_t bucket) {
                    const double weight = operators_.transport.cached_cohort_transport
                        ? static_cast<double>(
                              guidance_weight_cache_[bucket][location])
                        : guidance[bucket];
                    return weight;
                };
                const auto guidance_value = [&](std::size_t bucket) {
                    return raw_guidance_value(bucket) * (conditional_vacancy ? vacancies[bucket] : 1.0);
                };
                for (std::size_t bucket = 0;
                     bucket < active_direction_[stage].size(); ++bucket) {
                    const double mass =
                        active_direction_[stage][bucket][location];
                    if (!(mass > 0.0)) continue;
                    const double clock = std::max(0.0, static_cast<double>(
                        active_clock_[stage][bucket][location]));
                    if (mass < config_.migration.minimum_density) {
                        if (collect_incoming) continue;
                        active_work_[bucket][location] +=
                            static_cast<float>(mass);
                        clock_work_[bucket][location] +=
                            static_cast<float>(clock);
                        mark_sparse_target(location);
                        continue;
                    }
                    const double clock_per_mass = clock / mass;
                    double accepted_total = 0.0;
                    double reset_direction = 0.0;
                    const auto available = [&](std::size_t target_bucket) {
                        const Vec3i step = direction_vector(
                            direction_ids_[target_bucket - 1]);
                        if (!(x + step.x >= 0 && x + step.x < nx &&
                            y + step.y >= 0 && y + step.y < ny &&
                            z + step.z >= 0 && z + step.z < nz)) {
                            return false;
                        }
                        return conditional_vacancy ? vacancies[target_bucket] > 0.0 :
                            !vessel_blocks_cells(index(x + step.x, y + step.y, z + step.z));
                    };
                    const auto move = [&](std::size_t target_bucket,
                                          double direction_probability) {
                        const Vec3i step = direction_vector(
                            direction_ids_[target_bucket - 1]);
                        const int ox = x + step.x;
                        const int oy = y + step.y;
                        const int oz = z + step.z;
                        const std::size_t other = index(ox, oy, oz);
                        if (vessel_blocks_cells(other)) return;
                        const double occupied = std::clamp(
                            occupied_fraction(other) /
                                continuum.reaction.maximum_occupied_fraction,
                            0.0, 1.0);
                        const double availability = conditional_vacancy ? 1.0 : std::pow(
                            std::max(0.0, 1.0 - occupied),
                            vacancy_exponent);
                        double moved = mass * local_probability *
                            direction_probability * availability;
                        if (!(moved > 0.0)) return;
                        if (conditional_vacancy) {
                            auto& demand = proposed_incoming[incoming_offset(ox, oy, oz)];
                            if (collect_incoming) { demand += moved; return; }
                            const double capacity = std::max(0.0,
                                continuum.reaction.maximum_occupied_fraction - occupied_fraction(other)) /
                                (stage == 1 ? large_cell_volume_ : 1.0);
                            moved *= demand > 0.0 ? std::min(1.0, capacity / demand) : 0.0;
                        }
                        if (renewal_) renewal_->transfer(location, other, stage, moved);
                        if (duration_) duration_->transfer(location, other, stage, moved);
                        if (velocity_) velocity_->transfer(location, other, stage, moved);
                        accepted_total += moved;
                        active_work_[target_bucket][other] +=
                            static_cast<float>(moved);
                        clock_work_[target_bucket][other] +=
                            static_cast<float>(moved * clock_per_mass);
                        mark_sparse_target(other);
                    };
                    if (conditional_vacancy && bucket != 0) {
                        // Average the literal ABM choice over independent
                        // target-occupancy subsets. A blocked forward choice
                        // transfers its prior to feasible turns; normalizing
                        // mean vacancy weights would retain too much straight
                        // motion through a partially occupied neighbourhood.
                        std::vector<std::size_t> choices, uncertain;
                        if (available(bucket) && raw_guidance_value(bucket) > 0) choices.push_back(bucket);
                        for (auto turn : turn_buckets_[bucket])
                            if (available(turn) && raw_guidance_value(turn) > 0) choices.push_back(turn);
                        for (auto choice : choices) if (vacancies[choice] < 1.0) uncertain.push_back(choice);
                        std::array<double, 27> probabilities{};
                        for (std::size_t mask = 0; mask < (std::size_t(1) << uncertain.size()); ++mask) {
                            std::array<bool, 27> selected{};
                            for (auto choice : choices) selected[choice] = vacancies[choice] == 1.0;
                            double probability = 1.0;
                            for (std::size_t bit = 0; bit < uncertain.size(); ++bit) {
                                const auto choice = uncertain[bit];
                                selected[choice] = (mask >> bit) & 1U;
                                probability *= selected[choice] ? vacancies[choice] : 1 - vacancies[choice];
                            }
                            if (!(probability > 0.0)) continue;
                            std::size_t turns = 0;
                            for (auto choice : choices) if (choice != bucket && selected[choice]) ++turns;
                            std::array<double, 27> weights{};
                            double sum = 0.0;
                            for (auto choice : choices) {
                                if (!selected[choice]) continue;
                                const double prior = choice == bucket ? (turns ? continuum.base.continue_probability : 1.0) :
                                    (selected[bucket] ? 1 - continuum.base.continue_probability : 1.0) / turns;
                                weights[choice] = prior * raw_guidance_value(choice);
                                sum += weights[choice];
                            }
                            if (sum > 0.0) for (auto choice : choices)
                                probabilities[choice] += probability * weights[choice] / sum;
                        }
                        double choice_probability = 0.0;
                        for (auto choice : choices) {
                            choice_probability += probabilities[choice];
                            if (probabilities[choice] > 0.0) move(choice, probabilities[choice]);
                        }
                        // The ABM clears direction history after an attempted
                        // jump has no feasible forward/turn choice. Keep that
                        // mass at this site, ready to choose any direction on
                        // its next attempt, rather than trapping its history.
                        reset_direction = mass * local_probability *
                            std::clamp(1.0 - choice_probability, 0.0, 1.0);
                    } else if (bucket != 0) {
                        if (resource_guided) {
                            const bool forward_available = available(bucket) &&
                                guidance_value(bucket) > 0.0;
                            std::array<std::size_t, 26> valid_turns{};
                            std::size_t valid_turn_count = 0;
                            for (const std::size_t turn :
                                 turn_buckets_[bucket]) {
                                if (available(turn) &&
                                    guidance_value(turn) > 0.0) {
                                    valid_turns[valid_turn_count++] = turn;
                                }
                            }
                            double total_weight = 0.0;
                            double forward_weight = 0.0;
                            if (forward_available) {
                                forward_weight = valid_turn_count == 0
                                    ? guidance_value(bucket)
                                    : continuum.base.continue_probability *
                                        guidance_value(bucket);
                                total_weight += forward_weight;
                            }
                            const double turn_prior = valid_turn_count == 0
                                ? 0.0
                                : (forward_available
                                       ? 1.0 - continuum.base.continue_probability
                                       : 1.0) /
                                    static_cast<double>(valid_turn_count);
                            for (std::size_t turn_index = 0;
                                 turn_index < valid_turn_count; ++turn_index) {
                                const std::size_t turn = valid_turns[turn_index];
                                total_weight +=
                                    turn_prior * guidance_value(turn);
                            }
                            if (total_weight > 0.0) {
                                if (forward_weight > 0.0) {
                                    move(bucket, forward_weight / total_weight);
                                }
                                for (std::size_t turn_index = 0;
                                     turn_index < valid_turn_count; ++turn_index) {
                                    const std::size_t turn =
                                        valid_turns[turn_index];
                                    move(turn,
                                         turn_prior * guidance_value(turn) /
                                             total_weight);
                                }
                            }
                        } else {
                            const bool forward_available = available(bucket);
                            std::size_t valid_turns = 0;
                            for (const std::size_t turn :
                                 turn_buckets_[bucket]) {
                                if (available(turn)) ++valid_turns;
                            }
                            if (forward_available) {
                                move(bucket, valid_turns == 0 ? 1.0
                                    : continuum.base.continue_probability);
                            }
                            if (valid_turns > 0) {
                                const double total_turn_probability =
                                    forward_available
                                    ? 1.0 - continuum.base.continue_probability
                                    : 1.0;
                                const double turn_probability =
                                    total_turn_probability / valid_turns;
                                for (const std::size_t turn :
                                     turn_buckets_[bucket]) {
                                    if (available(turn)) {
                                        move(turn, turn_probability);
                                    }
                                }
                            }
                        }
                    } else if (conditional_vacancy) {
                        double weight_sum = 0.0;
                        for (std::size_t b = 1; b < active_direction_[stage].size(); ++b) {
                            weight_sum += guidance_value(b);
                        }
                        if (weight_sum > 0.0)
                            for (std::size_t b = 1; b < active_direction_[stage].size(); ++b)
                                if (guidance_value(b) > 0.0) move(b, guidance_value(b) / weight_sum);
                    }
                    const double retained =
                        std::max(0.0, mass - accepted_total);
                    if (collect_incoming) continue;
                    active_work_[bucket][location] +=
                        static_cast<float>(std::max(0.0, retained - reset_direction));
                    clock_work_[bucket][location] +=
                        static_cast<float>(std::max(0.0, retained - reset_direction) * clock_per_mass);
                    if (reset_direction > 0.0) {
                        active_work_[0][location] += static_cast<float>(reset_direction);
                        clock_work_[0][location] += static_cast<float>(reset_direction * clock_per_mass);
                    }
                    if (retained > 0.0) mark_sparse_target(location);
                }
            };

            if (conditional_vacancy) {
                // Reserve target vacancy against all simultaneous proposals.
                // Rejected incoming flux remains at its source, conserving
                // mass without relying on the reaction capacity clamp.
                collect_incoming = true;
                for (int z = old_bounds.z0; z < old_bounds.z1; ++z)
                    for (int y = old_bounds.y0; y < old_bounds.y1; ++y)
                        for (int x = old_bounds.x0; x < old_bounds.x1; ++x)
                            transport_location(x, y, z, index(x, y, z));
                collect_incoming = false;
            }

            if (sparse_transport) {
                for (const std::size_t location : sparse_locations) {
                    const int x = static_cast<int>(
                        location % static_cast<std::size_t>(nx));
                    const std::size_t yz =
                        location / static_cast<std::size_t>(nx);
                    const int y = static_cast<int>(
                        yz % static_cast<std::size_t>(ny));
                    const int z = static_cast<int>(
                        yz / static_cast<std::size_t>(ny));
                    transport_location(x, y, z, location);
                }
            } else if (operators_.transport.cached_cohort_transport &&
                       continuum.base.thin_layer) {
                // Sources on rows separated by three voxels have disjoint
                // Moore-neighbourhood write regions. Process each of the
                // three row colours in parallel without atomics; x remains
                // ordered within a row, preserving the local accumulation
                // order. This is algebraically the same transport operator.
                const int workers = (renewal_ || duration_) ? 1 : std::max(1, std::min(
                    continuum.base.threads, available_worker_threads()));
                for (int row_colour = 0; row_colour < 3; ++row_colour) {
                    const int first_y = old_bounds.y0 + row_colour;
                    if (first_y >= old_bounds.y1) continue;
                    const std::size_t row_count = static_cast<std::size_t>(
                        (old_bounds.y1 - first_y + 2) / 3);
                    const auto transport_row = [&](std::size_t row_offset) {
                        const int y = first_y +
                            3 * static_cast<int>(row_offset);
                        const int z = old_bounds.z0;
                        for (int x = old_bounds.x0;
                             x < old_bounds.x1; ++x) {
                            transport_location(x, y, z, index(x, y, z));
                        }
                    };
                    // Interleave rows across workers. Central rows contain
                    // more active mass and therefore more directional work;
                    // contiguous static chunks leave several workers idle at
                    // the barrier even though the row write sets are disjoint.
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1) num_threads(workers)
                    for (std::int64_t row_offset = 0;
                         row_offset < static_cast<std::int64_t>(row_count);
                         ++row_offset) {
                        transport_row(static_cast<std::size_t>(row_offset));
                    }
#else
                    for (std::size_t row_offset = 0;
                         row_offset < row_count; ++row_offset) {
                        transport_row(row_offset);
                    }
#endif
                }
            } else {
                for (int z = old_bounds.z0; z < old_bounds.z1; ++z) {
                    for (int y = old_bounds.y0; y < old_bounds.y1; ++y) {
                        for (int x = old_bounds.x0; x < old_bounds.x1; ++x) {
                            transport_location(x, y, z, index(x, y, z));
                        }
                    }
                }
            }
            active_direction_[stage].swap(active_work_);
            active_clock_[stage].swap(clock_work_);
            if (!population_bounds_.valid) {
                population_bounds_ = target_bounds;
            } else {
                population_bounds_.x0 = std::min(
                    population_bounds_.x0, target_bounds.x0);
                population_bounds_.y0 = std::min(
                    population_bounds_.y0, target_bounds.y0);
                population_bounds_.z0 = std::min(
                    population_bounds_.z0, target_bounds.z0);
                population_bounds_.x1 = std::max(
                    population_bounds_.x1, target_bounds.x1);
                population_bounds_.y1 = std::max(
                    population_bounds_.y1, target_bounds.y1);
                population_bounds_.z1 = std::max(
                    population_bounds_.z1, target_bounds.z1);
            }
            if (sparse_transport) {
                // After the swap, active_work_/clock_work_ contain the old
                // input. Clear exactly those source sites so the work fields
                // are globally clean for the next microstep.
                for (const std::size_t location : sparse_locations) {
                    active_total_[stage][location] = 0.0F;
                    for (std::size_t bucket = 0;
                         bucket < active_work_.size(); ++bucket) {
                        active_work_[bucket][location] = 0.0F;
                        clock_work_[bucket][location] = 0.0F;
                    }
                }
                work_dirty_bounds_ = {};
                std::sort(sparse_targets.begin(), sparse_targets.end());

                std::vector<std::size_t> remaining_locations;
                remaining_locations.reserve(sparse_targets.size());
                StructuredActiveBounds3D next_bounds;
                const double epsilon = config_.migration.minimum_density;
                for (const std::size_t location : sparse_targets) {
                    double total = 0.0;
                    for (std::size_t bucket = 0;
                         bucket < active_direction_[stage].size(); ++bucket) {
                        auto& mass =
                            active_direction_[stage][bucket][location];
                        auto& clock_value =
                            active_clock_[stage][bucket][location];
                        const double mass_value = mass;
                        const double clock = std::max(
                            0.0, static_cast<double>(clock_value));
                        if (mass_value <= epsilon ||
                            clock <= mass_value * sub_dt +
                                1.0e-6 * std::max(1.0, clock)) {
                            if (mass_value > 0.0) {
                                r_normal_[stage][location] += mass_value;
                            }
                            mass = 0.0F;
                            clock_value = 0.0F;
                        } else {
                            clock_value = static_cast<float>(
                                clock - mass_value * sub_dt);
                            total += mass_value;
                        }
                    }
                    active_total_[stage][location] =
                        static_cast<float>(total);
                    if (total < epsilon) continue;
                    remaining_locations.push_back(location);
                    const int x = static_cast<int>(
                        location % static_cast<std::size_t>(nx));
                    const std::size_t yz =
                        location / static_cast<std::size_t>(nx);
                    const int y = static_cast<int>(
                        yz % static_cast<std::size_t>(ny));
                    const int z = static_cast<int>(
                        yz / static_cast<std::size_t>(ny));
                    if (!next_bounds.valid) {
                        next_bounds = {x, y, z, x + 1, y + 1, z + 1, true};
                    } else {
                        next_bounds.x0 = std::min(next_bounds.x0, x);
                        next_bounds.y0 = std::min(next_bounds.y0, y);
                        next_bounds.z0 = std::min(next_bounds.z0, z);
                        next_bounds.x1 = std::max(next_bounds.x1, x + 1);
                        next_bounds.y1 = std::max(next_bounds.y1, y + 1);
                        next_bounds.z1 = std::max(next_bounds.z1, z + 1);
                    }
                }
                sparse_locations.swap(remaining_locations);
                active_bounds_[stage] = next_bounds;
            } else {
                work_dirty_bounds_ = old_bounds;
                active_bounds_[stage] = target_bounds;
                finish_duration_transport();
                finish_velocity_transport();
                expire_active(stage, sub_dt);
                shrink_active_bounds(stage);
            }
            finish_division_transport();
            if (!active_bounds_[stage].valid) break;
        }
    }
}

}  // namespace atcg3d::structured_pde
