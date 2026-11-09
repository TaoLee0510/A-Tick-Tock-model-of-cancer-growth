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

void StructuredPdeModel3D::exchange_active_r_with_K(double dt) {
    if (!operators_.exchange.direction_flux ||
        config_.migration.crowding_exchange !=
            "active_r_K_stage1_conservative_v1" ||
        !active_bounds_[0].valid) {
        return;
    }

    // This is the structured-PDE closure of the ABM singleton swap. Only
    // small active-r and small K participate. A swap moves equal mass in
    // opposite directions, so both global species mass and the occupied
    // fraction at each endpoint are conserved by this operator.
    struct SwapProposal {
        std::size_t source{};
        std::size_t target{};
        std::size_t direction_bucket{};
        double requested{};
        double clock_per_mass{};
    };

    const auto& continuum = config_.continuum;
    const int nx = continuum.grid.shape[0];
    const int ny = continuum.grid.shape[1];
    const int nz = continuum.grid.shape[2];
    const auto bounds = active_bounds_[0];
    const double rate = active_rate(StructuredStage3D::small);
    const double move_probability = 1.0 - std::exp(-rate * dt);
    const double maximum = continuum.reaction.maximum_occupied_fraction;
    const double epsilon = config_.migration.minimum_density;
    std::vector<SwapProposal> proposals;
    std::unordered_map<std::size_t, double> demand_by_target;

    for (int z = bounds.z0; z < bounds.z1; ++z) {
        for (int y = bounds.y0; y < bounds.y1; ++y) {
            for (int x = bounds.x0; x < bounds.x1; ++x) {
                const std::size_t source = index(x, y, z);
                const double active = active_total_[0][source];
                if (active < epsilon || vessel_blocks_cells(source)) continue;

                // The ABM only enters its swap path when no empty direction is
                // available. In a density closure, the minimum neighbour
                // occupancy is a continuous approximation of that condition:
                // one empty neighbour makes blockedness zero; a fully packed
                // neighbourhood makes it one.
                double maximum_vacancy = 0.0;
                bool has_cell_site_neighbor = false;
                for (const DirectionId direction : direction_ids_) {
                    const Vec3i step = direction_vector(direction);
                    const int tx = x + step.x;
                    const int ty = y + step.y;
                    const int tz = z + step.z;
                    if (tx < 0 || tx >= nx || ty < 0 || ty >= ny ||
                        tz < 0 || tz >= nz) continue;
                    const std::size_t target = index(tx, ty, tz);
                    if (vessel_blocks_cells(target)) continue;
                    has_cell_site_neighbor = true;
                    const double occupied = std::clamp(
                        occupied_fraction(target) / maximum, 0.0, 1.0);
                    maximum_vacancy = std::max(maximum_vacancy, 1.0 - occupied);
                }
                if (!has_cell_site_neighbor) continue;
                const double blockedness = 1.0 - maximum_vacancy;
                if (blockedness <= epsilon) continue;

                const std::vector<double> weights =
                    guided_direction_weights(source, true);
                double weight_sum = 0.0;
                for (std::size_t bucket = 1; bucket < weights.size(); ++bucket) {
                    if (!(weights[bucket] > 0.0)) continue;
                    const Vec3i step = direction_vector(direction_ids_[bucket - 1]);
                    const std::size_t target = index(
                        x + step.x, y + step.y, z + step.z);
                    if (K_[0][target] >= epsilon) weight_sum += weights[bucket];
                }
                if (!(weight_sum > 0.0)) continue;

                double clock = 0.0;
                for (const auto& field : active_clock_[0]) clock += field[source];
                const double clock_per_mass = clock / active;
                for (std::size_t bucket = 1; bucket < weights.size(); ++bucket) {
                    if (!(weights[bucket] > 0.0)) continue;
                    const Vec3i step = direction_vector(direction_ids_[bucket - 1]);
                    const std::size_t target = index(
                        x + step.x, y + step.y, z + step.z);
                    const double partner = std::clamp(K_[0][target], 0.0, 1.0);
                    if (partner < epsilon) continue;
                    const double requested = active * move_probability *
                        blockedness * weights[bucket] / weight_sum * partner;
                    if (requested < epsilon) continue;
                    proposals.push_back(
                        {source, target, bucket, requested, clock_per_mass});
                    demand_by_target[target] += requested;
                }
            }
        }
    }
    if (proposals.empty()) return;
    if (renewal_) renewal_->begin_transport();
    if (duration_) duration_->begin_transport();
    if (velocity_) velocity_->begin_transport();

    std::unordered_map<std::size_t, double> outgoing_by_source;
    std::unordered_map<std::size_t, double> K_delta;
    StructuredActiveBounds3D clear_bounds = expanded_bounds(bounds);
    if (work_dirty_bounds_.valid) {
        clear_bounds.x0 = std::min(clear_bounds.x0, work_dirty_bounds_.x0);
        clear_bounds.y0 = std::min(clear_bounds.y0, work_dirty_bounds_.y0);
        clear_bounds.z0 = std::min(clear_bounds.z0, work_dirty_bounds_.z0);
        clear_bounds.x1 = std::max(clear_bounds.x1, work_dirty_bounds_.x1);
        clear_bounds.y1 = std::max(clear_bounds.y1, work_dirty_bounds_.y1);
        clear_bounds.z1 = std::max(clear_bounds.z1, work_dirty_bounds_.z1);
    }
    clear_active_work(clear_bounds);
    for (const SwapProposal& proposal : proposals) {
        const double demand = demand_by_target.at(proposal.target);
        const double target_scale = demand > K_[0][proposal.target]
            ? K_[0][proposal.target] / demand : 1.0;
        const double accepted = proposal.requested * target_scale;
        if (!(accepted > 0.0)) continue;
        if (duration_) duration_->transfer(proposal.source, proposal.target, 0, accepted);
        if (velocity_) velocity_->transfer(proposal.source, proposal.target, 0, accepted);
        if (renewal_) {
            renewal_->transfer(proposal.source, proposal.target, 0, accepted);
            renewal_->transfer(proposal.target, proposal.source, 2, accepted);
        }
        outgoing_by_source[proposal.source] += accepted;
        K_delta[proposal.source] += accepted;
        K_delta[proposal.target] -= accepted;
        active_work_[proposal.direction_bucket][proposal.target] +=
            static_cast<float>(accepted);
        clock_work_[proposal.direction_bucket][proposal.target] +=
            static_cast<float>(accepted * proposal.clock_per_mass);
    }

    const auto target_bounds = expanded_bounds(bounds);
    for (int z = target_bounds.z0; z < target_bounds.z1; ++z) {
        for (int y = target_bounds.y0; y < target_bounds.y1; ++y) {
            for (int x = target_bounds.x0; x < target_bounds.x1; ++x) {
                const std::size_t location = index(x, y, z);
                const auto found = outgoing_by_source.find(location);
                const double active = active_total_[0][location];
                const double fraction = found == outgoing_by_source.end() ||
                    !(active > 0.0) ? 0.0 :
                    std::clamp(found->second / active, 0.0, 1.0);
                double total = 0.0;
                for (std::size_t bucket = 0;
                     bucket < active_direction_[0].size(); ++bucket) {
                    const double mass =
                        active_direction_[0][bucket][location] * (1.0 - fraction) +
                        active_work_[bucket][location];
                    const double clock =
                        active_clock_[0][bucket][location] * (1.0 - fraction) +
                        clock_work_[bucket][location];
                    active_direction_[0][bucket][location] =
                        static_cast<float>(std::max(0.0, mass));
                    active_clock_[0][bucket][location] =
                        static_cast<float>(std::max(0.0, clock));
                    total += std::max(0.0, mass);
                }
                active_total_[0][location] = static_cast<float>(total);
            }
        }
    }
    for (const auto& [location, delta] : K_delta) {
        K_[0][location] = std::max(0.0, K_[0][location] + delta);
    }
    work_dirty_bounds_ = clear_bounds;
    active_bounds_[0] = target_bounds;
    shrink_active_bounds(0);
    finish_division_transport();
    finish_duration_transport();
    finish_velocity_transport();
}

}  // namespace atcg3d::structured_pde
