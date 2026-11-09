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

double StructuredPdeModel3D::mean_growth_rate(CellType type) const {
    const auto& base = config_.continuum.base;
    if (operators_.reaction.true_normal_mean) {
        return expected_initial_growth_rate(base, type);
    }
    if (base.initial_growth_rate_model == "fixed") {
        return type == CellType::r ? base.initial_r_growth_rate
                                   : base.initial_K_growth_rate;
    }
    const auto& distribution = type == CellType::r
        ? base.initial_r_growth_truncated_normal
        : base.initial_K_growth_truncated_normal;
    return std::clamp(distribution.mean, distribution.minimum,
                      distribution.maximum);
}

void StructuredPdeModel3D::build_local_counts(
    std::vector<double>& r_counts,
    std::vector<double>& K_counts) const {
    const auto& continuum = config_.continuum;
    if (!population_bounds_.valid) {
        r_counts.clear();
        K_counts.clear();
        return;
    }
    const auto bounds = population_bounds_;
    const int width = bounds.x1 - bounds.x0;
    const int height = bounds.y1 - bounds.y0;
    const int depth = bounds.z1 - bounds.z0;
    const std::size_t box_size = static_cast<std::size_t>(width) * height * depth;
    const int radius = std::max(0, static_cast<int>(std::floor(
        0.5 * continuum.base.growth_density_window_edge /
        continuum.grid.spacing_voxels)));
    const int lower = operators_.reaction.exact_growth_window
        ? static_cast<int>(std::floor(((continuum.base.growth_density_window_edge - 1) / 2) / continuum.grid.spacing_voxels))
        : radius;
    const int upper = operators_.reaction.exact_growth_window
        ? static_cast<int>(std::floor((continuum.base.growth_density_window_edge -
            (continuum.base.growth_density_window_edge - 1) / 2 - 1) / continuum.grid.spacing_voxels))
        : radius;
    if (continuum.base.thin_layer) {
        const int pitch = width + 1;
        const std::size_t prefix_size =
            static_cast<std::size_t>(pitch) * (height + 1);
        std::vector<double> r_prefix(prefix_size, 0.0);
        std::vector<double> K_prefix(prefix_size, 0.0);
        for (int y = 1; y <= height; ++y) {
            double r_row = 0.0;
            double K_row = 0.0;
            for (int x = 1; x <= width; ++x) {
                const std::size_t source =
                    index(bounds.x0 + x - 1, bounds.y0 + y - 1, bounds.z0);
                r_row += (r_normal_[0][source] + r_normal_[1][source] +
                    r_active(StructuredStage3D::small, source) +
                    r_active(StructuredStage3D::large, source)) * voxel_measure_;
                K_row += (K_[0][source] + K_[1][source]) * voxel_measure_;
                if(!external_r_.empty()) { r_row+=external_r_[source]*voxel_measure_;K_row+=external_K_[source]*voxel_measure_; }
                const std::size_t here = static_cast<std::size_t>(y) * pitch + x;
                r_prefix[here] = r_row;
                K_prefix[here] = K_row;
            }
        }
        for (int x = 1; x <= width; ++x) {
            for (int y = 1; y <= height; ++y) {
                const std::size_t here = static_cast<std::size_t>(y) * pitch + x;
                r_prefix[here] += r_prefix[here - pitch];
                K_prefix[here] += K_prefix[here - pitch];
            }
        }
        const auto sum = [pitch](const std::vector<double>& prefix,
                                 int x0, int y0, int x1, int y1) {
            return prefix[static_cast<std::size_t>(y1) * pitch + x1]
                - prefix[static_cast<std::size_t>(y0) * pitch + x1]
                - prefix[static_cast<std::size_t>(y1) * pitch + x0]
                + prefix[static_cast<std::size_t>(y0) * pitch + x0];
        };
        r_counts.assign(box_size, 0.0);
        K_counts.assign(box_size, 0.0);
        const int workers = std::max(
            1, std::min(continuum.base.threads, available_worker_threads()));
        deterministic_parallel_for(box_size, workers, [&](std::size_t offset) {
            const int local_x = static_cast<int>(offset %
                static_cast<std::size_t>(width));
            const int local_y = static_cast<int>(offset /
                static_cast<std::size_t>(width));
            const std::size_t here = index(
                bounds.x0 + local_x, bounds.y0 + local_y, bounds.z0);
            if (!operators_.reaction.exact_growth_window && occupied_fraction(here) <= 0.0) return;
            const int x0 = std::max(0, local_x - lower);
            const int y0 = std::max(0, local_y - lower);
            const int x1 = std::min(width, local_x + upper + 1);
            const int y1 = std::min(height, local_y + upper + 1);
            r_counts[offset] = sum(r_prefix, x0, y0, x1, y1);
            K_counts[offset] = sum(K_prefix, x0, y0, x1, y1);
        });
        return;
    }

    const int px = width + 1;
    const int py = height + 1;
    const int pz = depth + 1;
    const auto pindex = [px, py](int x, int y, int z) {
        return (static_cast<std::size_t>(z) * py + y) * px + x;
    };
    std::vector<double> r_prefix(static_cast<std::size_t>(px) * py * pz, 0.0);
    std::vector<double> K_prefix(static_cast<std::size_t>(px) * py * pz, 0.0);
    for (int z = 1; z <= depth; ++z) {
        for (int y = 1; y <= height; ++y) {
            for (int x = 1; x <= width; ++x) {
                const std::size_t source = index(
                    bounds.x0 + x - 1, bounds.y0 + y - 1,
                    bounds.z0 + z - 1);
                double r_value = (r_normal_[0][source] + r_normal_[1][source] +
                    r_active(StructuredStage3D::small, source) +
                    r_active(StructuredStage3D::large, source)) * voxel_measure_;
                double K_value = (K_[0][source] + K_[1][source]) * voxel_measure_;
                if(!external_r_.empty()) { r_value+=external_r_[source]*voxel_measure_;K_value+=external_K_[source]*voxel_measure_; }
                const auto update = [&](std::vector<double>& prefix, double value) {
                    prefix[pindex(x, y, z)] = value
                        + prefix[pindex(x - 1, y, z)]
                        + prefix[pindex(x, y - 1, z)]
                        + prefix[pindex(x, y, z - 1)]
                        - prefix[pindex(x - 1, y - 1, z)]
                        - prefix[pindex(x - 1, y, z - 1)]
                        - prefix[pindex(x, y - 1, z - 1)]
                        + prefix[pindex(x - 1, y - 1, z - 1)];
                };
                update(r_prefix, r_value);
                update(K_prefix, K_value);
            }
        }
    }
    const auto sum = [&](const std::vector<double>& prefix,
                         int x0, int y0, int z0,
                         int x1, int y1, int z1) {
        return prefix[pindex(x1, y1, z1)]
            - prefix[pindex(x0, y1, z1)] - prefix[pindex(x1, y0, z1)]
            - prefix[pindex(x1, y1, z0)] + prefix[pindex(x0, y0, z1)]
            + prefix[pindex(x0, y1, z0)] + prefix[pindex(x1, y0, z0)]
            - prefix[pindex(x0, y0, z0)];
    };
    r_counts.assign(box_size, 0.0);
    K_counts.assign(box_size, 0.0);
    for (int z = 0; z < depth; ++z) {
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                const int x0 = std::max(0, x - lower);
                const int y0 = std::max(0, y - lower);
                const int z0 = std::max(0, z - lower);
                const int x1 = std::min(width, x + upper + 1);
                const int y1 = std::min(height, y + upper + 1);
                const int z1 = std::min(depth, z + upper + 1);
                const std::size_t offset =
                    (static_cast<std::size_t>(z) * height + y) * width + x;
                r_counts[offset] = sum(r_prefix, x0, y0, z0, x1, y1, z1);
                K_counts[offset] = sum(K_prefix, x0, y0, z0, x1, y1, z1);
            }
        }
    }
}

std::vector<StructuredPdeModel3D::SmallBirth3D> StructuredPdeModel3D::prepare_small_births(
    const StructuredActiveBounds3D& bounds) const {
    const auto& grid = config_.continuum.grid;
    const auto width = static_cast<std::size_t>(bounds.x1 - bounds.x0);
    const auto height = static_cast<std::size_t>(bounds.y1 - bounds.y0);
    const auto depth = static_cast<std::size_t>(bounds.z1 - bounds.z0);
    std::vector<SmallBirth3D> result(width * height * depth);
    for (int z = bounds.z0; z < bounds.z1; ++z) {
        for (int y = bounds.y0; y < bounds.y1; ++y) {
            for (int x = bounds.x0; x < bounds.x1; ++x) {
                const auto here = index(x, y, z);
                if (r_normal_[0][here] + r_active(StructuredStage3D::small, here) + K_[0][here] <= 0.0) continue;
                if (vessel_blocks_cells(here)) continue;
                std::vector<double> availability(direction_ids_.size(), 0.0);
                for (std::size_t d = 0; d < direction_ids_.size(); ++d) {
                    const auto direction = direction_vector(direction_ids_[d]);
                    const int ox = x + direction.x;
                    const int oy = y + direction.y;
                    const int oz = z + direction.z;
                    if (ox < 0 || ox >= grid.shape[0] || oy < 0 || oy >= grid.shape[1] ||
                        oz < 0 || oz >= grid.shape[2]) continue;
                    const auto other = index(ox, oy, oz);
                    if (transport_destination_blocks_cells(other, StructuredStage3D::small)) continue;
                    availability[d] = 1.0 - std::clamp(occupied_fraction(other) /
                        config_.continuum.reaction.maximum_occupied_fraction, 0.0, 1.0);
                }
                const auto probabilities = uniform_feasible_jump_probabilities(availability);
                const auto offset = (static_cast<std::size_t>(z - bounds.z0) * height +
                    static_cast<std::size_t>(y - bounds.y0)) * width + static_cast<std::size_t>(x - bounds.x0);
                for (std::size_t d = 0; d < direction_ids_.size(); ++d) {
                    result[offset].probabilities[direction_ids_[d]] = probabilities[d];
                    result[offset].success += probabilities[d];
                }
                result[offset].success = std::clamp(result[offset].success, 0.0, 1.0);
            }
        }
    }
    return result;
}

void StructuredPdeModel3D::place_small_births(const StructuredActiveBounds3D& bounds,
                                            const std::vector<SmallBirth3D>& births) {
    const auto& grid = config_.continuum.grid;
    const auto width = static_cast<std::size_t>(bounds.x1 - bounds.x0);
    const auto height = static_cast<std::size_t>(bounds.y1 - bounds.y0);
    std::map<std::size_t, std::array<double, 2>> incoming;
    for (std::size_t offset = 0; offset < births.size(); ++offset) {
        const auto& birth = births[offset];
        if (birth.mass[0] + birth.mass[1] <= 0.0) continue;
        const int x = bounds.x0 + static_cast<int>(offset % width);
        const int y = bounds.y0 + static_cast<int>((offset / width) % height);
        const int z = bounds.z0 + static_cast<int>(offset / (width * height));
        for (const auto id : direction_ids_) {
            const double probability = birth.probabilities[id];
            if (probability <= 0.0) continue;
            const auto direction = direction_vector(id);
            auto& target = incoming[index(x + direction.x, y + direction.y, z + direction.z)];
            for (std::size_t type = 0; type < 2; ++type) target[type] += birth.mass[type] * probability;
        }
    }
    for (const auto& [here, mass] : incoming) {
        const double available = std::max(0.0, config_.continuum.reaction.maximum_occupied_fraction -
            occupied_fraction(here));
        const double total = mass[0] + mass[1];
        const double scale = total > available ? available / total : 1.0;
        const double r_birth = scale * mass[0];
        const double K_birth = scale * mass[1];
        r_normal_[0][here] += r_birth;
        K_[0][here] += K_birth;
        if (renewal_) {
            if (r_birth > 0.0) renewal_->add_fresh(here, 0, r_birth);
            if (K_birth > 0.0) renewal_->add_fresh(here, 2, K_birth);
        }
        const int x = static_cast<int>(here % grid.shape[0]);
        const int y = static_cast<int>((here / grid.shape[0]) % grid.shape[1]);
        const int z = static_cast<int>(here / (static_cast<std::size_t>(grid.shape[0]) * grid.shape[1]));
        if (r_birth + K_birth > 0.0) include_population_location(x, y, z);
    }
}

void StructuredPdeModel3D::react(double dt) {
    if (!population_bounds_.valid) return;
    std::vector<double> r_counts;
    std::vector<double> K_counts;
    build_local_counts(r_counts, K_counts);
    const auto& continuum = config_.continuum;
    const double r_inherent = mean_growth_rate(CellType::r);
    const double K_inherent = mean_growth_rate(CellType::K);
    const double legacy_r_limit = continuum.base.thin_layer
        ? continuum.base.legacy_mapping.source_r_limit : continuum.base.r_limit;
    const double legacy_K_limit = continuum.base.thin_layer
        ? continuum.base.legacy_mapping.source_K_limit : continuum.base.K_limit;
    const double legacy_r_capacity = continuum.base.thin_layer
        ? continuum.base.legacy_mapping.source_carrying_capacity_r
        : continuum.base.carrying_capacity_r;
    const double legacy_K_capacity = continuum.base.thin_layer
        ? continuum.base.legacy_mapping.source_carrying_capacity_K
        : continuum.base.carrying_capacity_K;
    const double r_limit = operators_.nutrient.transient_resources
        ? continuum.nutrient.common_density_limit : legacy_r_limit;
    const double K_limit = operators_.nutrient.transient_resources
        ? continuum.nutrient.common_density_limit : legacy_K_limit;
    const double r_capacity = operators_.nutrient.transient_resources
        ? continuum.nutrient.common_carrying_capacity : legacy_r_capacity;
    const double K_capacity = operators_.nutrient.transient_resources
        ? continuum.nutrient.common_carrying_capacity : legacy_K_capacity;
    const int workers = (renewal_ || duration_) ? 1 : std::max(
        1, std::min(continuum.base.threads, available_worker_threads()));
    const auto bounds = population_bounds_;
    const bool spatial_birth = operators_.reaction.neighbour_births &&
        config_.division_operator_enabled;
    auto births = spatial_birth ? prepare_small_births(bounds) : std::vector<SmallBirth3D>{};
    const std::size_t bx = static_cast<std::size_t>(bounds.x1 - bounds.x0);
    const std::size_t by = static_cast<std::size_t>(bounds.y1 - bounds.y0);
    const std::size_t bz = static_cast<std::size_t>(bounds.z1 - bounds.z0);
    deterministic_parallel_for(bx * by * bz, workers, [&](std::size_t offset) {
        const int x = bounds.x0 + static_cast<int>(offset % bx);
        const std::size_t yz = offset / bx;
        const int y = bounds.y0 + static_cast<int>(yz % by);
        const int z = bounds.z0 + static_cast<int>(yz / by);
        const std::size_t location = index(x, y, z);
        if (vessel_blocks_cells(location)) return;
        const double occupied = occupied_fraction(location);
        if (occupied <= 0.0) return;
        const double multiplier = capacity_multiplier(nutrient_[location]);
        const double r_count = r_counts[offset] / multiplier;
        const double K_count = K_counts[offset] / multiplier;
        const double total_count = r_count + K_count;
        const double raw_r_growth = calculate_density_growth_rate_continuous(
            static_cast<int>(CellType::r), r_inherent,
            r_count, K_count, total_count, r_limit, K_limit,
            continuum.base.alpha, continuum.base.beta, r_capacity, K_capacity);
        const double raw_K_growth = calculate_density_growth_rate_continuous(
            static_cast<int>(CellType::K), K_inherent,
            r_count, K_count, total_count, r_limit, K_limit,
            continuum.base.alpha, continuum.base.beta, r_capacity, K_capacity);
        const double nutrient_factor = operators_.nutrient.transient_resources
            ? nutrient_[location] /
                (continuum.nutrient.growth_half_saturation +
                 nutrient_[location])
            : 1.0;
        const auto resource_limited = [&](double growth) {
            return growth > 0.0 ? growth * nutrient_factor : growth;
        };
        const double r_growth = resource_limited(raw_r_growth);
        const double K_growth = resource_limited(raw_K_growth);
        const bool abm_work_clock = continuum.reaction.model ==
            "abm_work_clock_neighbor_availability_v2";
        const double r_division_rate = (config_.division_operator_enabled ? positive_part(r_growth) : 0.0) /
            (continuum.base.division_timing.base_cycle_hours *
             (abm_work_clock ? 1.0 : std::max(1.0e-12, r_inherent)));
        const double K_division_rate = (config_.division_operator_enabled ? positive_part(K_growth) : 0.0) /
            (continuum.base.division_timing.base_cycle_hours *
             (abm_work_clock ? 1.0 : std::max(1.0e-12, K_inherent)));
        const double r_death_rate =
            r_growth <= continuum.base.death_growth_rate_threshold
            ? 1.0 / continuum.base.r_death_delay_hours : 0.0;
        const double K_death_rate =
            K_growth <= continuum.base.death_growth_rate_threshold
            ? 1.0 / continuum.base.K_death_delay_hours : 0.0;
        const double vacancy = std::clamp(
            1.0 - occupied / continuum.reaction.maximum_occupied_fraction,
            0.0, 1.0);
        const double large_success = external_transport_destination_ &&
            transport_destination_blocks_cells(location, StructuredStage3D::large)
            ? 0.0 : std::pow(vacancy, continuum.reaction.large_daughter_vacancy_exponent);
        const double small_success = external_transport_destination_ &&
            !spatial_birth && !external_transport_destination_(location)
            ? 0.0 : spatial_birth ? births[offset].success : abm_work_clock
            ? 1.0 - std::pow(
                  1.0 - vacancy,
                  continuum.reaction.small_daughter_vacancy_exponent)
            : std::pow(
                  vacancy,
                  continuum.reaction.small_daughter_vacancy_exponent);
        const auto conversion_probability = [&](std::size_t stage) {
            if (!continuum.base.r_to_K_conversion.enabled) return 0.0;
            const double density = operators_.reaction.local_density_conversion
                ? activation_density_[stage][location] : occupied;
            return density >= continuum.base.r_to_K_conversion.density_threshold
                ? continuum.base.r_to_K_conversion.probability_per_division
                : 0.0;
        };

        const std::array<double, 2> old_r{
            r_normal_[0][location] + r_active(StructuredStage3D::small, location),
            r_normal_[1][location] + r_active(StructuredStage3D::large, location)};
        const std::array<double, 2> old_K{K_[0][location], K_[1][location]};
        std::array<double, 4> completed{};
        if (renewal_ && config_.division_operator_enabled) {
            for (std::size_t stage = 0; stage < 2; ++stage) {
                completed[stage] = renewal_->advance(location, stage, dt * positive_part(r_growth));
                completed[stage + 2] = renewal_->advance(location, stage + 2, dt * positive_part(K_growth));
            }
        }
        std::array<double, 2> delta_r{
            -r_death_rate * old_r[0], -r_death_rate * old_r[1]};
        std::array<double, 2> delta_K{
            -K_death_rate * old_K[0], -K_death_rate * old_K[1]};

        const double r_large_events = renewal_ ? completed[1] / dt : r_division_rate * old_r[1];
        const double r_large_daughters = large_success * r_large_events;
        const double r_shape_reductions =
            (1.0 - large_success) * r_large_events;
        const double large_conversion = conversion_probability(1);
        delta_r[1] += (1.0 - large_conversion) * r_large_daughters
            - r_shape_reductions;
        delta_K[1] += large_conversion * r_large_daughters;
        delta_r[0] += (2.0 - large_conversion) * r_shape_reductions;
        delta_K[0] += large_conversion * r_shape_reductions;

        const double K_large_events = renewal_ ? completed[3] / dt : K_division_rate * old_K[1];
        const double K_large_daughters = large_success * K_large_events;
        const double K_shape_reductions =
            (1.0 - large_success) * K_large_events;
        delta_K[1] += K_large_daughters - K_shape_reductions;
        delta_K[0] += 2.0 * K_shape_reductions;

        const double r_small_events = renewal_ ? completed[0] / dt : r_division_rate * old_r[0];
        const double r_small_births = small_success * r_small_events;
        const double small_conversion = conversion_probability(0);
        delta_r[0] += (1.0 - small_conversion) * r_small_births
            - (1.0 - small_success) * r_small_events *
                continuum.reaction.failed_r_division_death_fraction;
        delta_K[0] += small_conversion * r_small_births;
        delta_K[0] += small_success * (renewal_ ? completed[2] / dt : K_division_rate * old_K[0]);
        if (spatial_birth) {
            const double K_small_events = renewal_ ? completed[2] / dt : K_division_rate * old_K[0];
            births[offset].mass[0] = dt * (1.0 - small_conversion) * r_small_events;
            births[offset].mass[1] = dt * (small_conversion * r_small_events + K_small_events);
            delta_r[0] -= (1.0 - small_conversion) * r_small_births;
            delta_K[0] -= small_conversion * r_small_births + small_success * K_small_events;
        }

        std::array<double, 2> positive_r{};
        std::array<double, 2> positive_K{};
        double base_occupied = 0.0;
        double positive_occupied = 0.0;
        for (std::size_t stage = 0; stage < 2; ++stage) {
            const double r_negative = dt * std::min(0.0, delta_r[stage]);
            const double r_factor = old_r[stage] > 0.0
                ? std::clamp((old_r[stage] + r_negative) / old_r[stage], 0.0, 1.0)
                : 0.0;
            r_normal_[stage][location] *= r_factor;
            if (operators_.reaction.exact_growth_window) {
                r_refractory_[stage][location] *= r_factor;
                refractory_clock_[stage][location] *= r_factor;
            }
            for (std::size_t bucket = 0;
                 bucket < active_direction_[stage].size(); ++bucket) {
                active_direction_[stage][bucket][location] = static_cast<float>(
                    active_direction_[stage][bucket][location] * r_factor);
                active_clock_[stage][bucket][location] = static_cast<float>(
                    active_clock_[stage][bucket][location] * r_factor);
            }
            if (operators_.reaction.exact_growth_window) {
                double active_sum = 0.0;
                for (const auto& bucket : active_direction_[stage]) {
                    active_sum += bucket[location];
                }
                active_total_[stage][location] = static_cast<float>(active_sum);
            } else {
                active_total_[stage][location] = static_cast<float>(
                    active_total_[stage][location] * r_factor);
            }
            K_[stage][location] = std::max(
                0.0, old_K[stage] + dt * std::min(0.0, delta_K[stage]));
            positive_r[stage] = dt * std::max(0.0, delta_r[stage]);
            positive_K[stage] = dt * std::max(0.0, delta_K[stage]);
            const double volume = stage == 0 ? 1.0 : large_cell_volume_;
            base_occupied += volume *
                (r_normal_[stage][location] +
                 r_active(static_cast<StructuredStage3D>(stage), location) +
                 K_[stage][location]);
            positive_occupied += volume *
                (positive_r[stage] + positive_K[stage]);
        }
        const double available = std::max(
            0.0, continuum.reaction.maximum_occupied_fraction - base_occupied -
                (external_occupied_.empty()?0.0:external_occupied_[location]));
        const double scale = positive_occupied > available && positive_occupied > 0.0
            ? available / positive_occupied : 1.0;
        for (std::size_t stage = 0; stage < 2; ++stage) {
            // All newly created r mass begins in the ordinary state, matching
            // the ABM division-cycle reset before a later density refresh.
            r_normal_[stage][location] += scale * positive_r[stage];
            K_[stage][location] += scale * positive_K[stage];
            if (renewal_ || duration_) {
                // Successful mothers reset activity as well as their work.
                const double reset_fraction = old_r[stage] > 0.0
                    ? std::clamp(renewal_ ? completed[stage] / old_r[stage] : dt * r_division_rate, 0.0, 1.0) : 0.0;
                double active_sum = 0.0;
                for (std::size_t bucket = 0; bucket < active_direction_[stage].size(); ++bucket) {
                    const double before = active_direction_[stage][bucket][location];
                    active_direction_[stage][bucket][location] = static_cast<float>(before * (1.0 - reset_fraction));
                    active_clock_[stage][bucket][location] = static_cast<float>(
                        active_clock_[stage][bucket][location] * (1.0 - reset_fraction));
                    r_normal_[stage][location] += before - active_direction_[stage][bucket][location];
                    active_sum += active_direction_[stage][bucket][location];
                }
                active_total_[stage][location] = static_cast<float>(active_sum);
                if (duration_ && (active_sum > 0.0 || duration_->mass(location, stage) > 0.0))
                    duration_->reconcile(location, stage, active_sum);
                if (velocity_ && (active_sum > 0.0 || velocity_->mass(location, stage) > 0.0))
                    velocity_->reconcile(location, stage, active_sum);
                if (renewal_) renewal_->reconcile(location, stage, division_channel_mass(location, stage));
                if (renewal_ && stage == 0) {
                    // A K cell with no free daughter site retries; it does not
                    // draw a fresh biological cycle after a failed attempt.
                    const double retries = completed[2] * (1.0 - small_success * scale);
                    if (retries > 0.0) renewal_->add(location, 2, retries,
                        positive_part(K_growth) * continuum.base.division_timing.retry_delay_hours);
                }
                if (renewal_) renewal_->reconcile(location, stage + 2, K_[stage][location]);
            }
        }
    });
    if (spatial_birth) place_small_births(bounds, births);
}

}  // namespace atcg3d::structured_pde
