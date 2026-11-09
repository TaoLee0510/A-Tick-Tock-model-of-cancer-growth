#include "model/operator_strategies.hpp"
#include "model/structured_pde_model.hpp"

#include <algorithm>
#include <cmath>

namespace {
bool same_time(double lhs, double rhs) noexcept {
    return std::abs(lhs - rhs) <=
        1.0e-10 * std::max({1.0, std::abs(lhs), std::abs(rhs)});
}
}  // namespace

namespace atcg3d::structured_pde {
namespace {
std::uint32_t select_checkpoint_version(int schema) noexcept {
    if (schema >= 17) {
        return 10U;
    }
    if (schema >= 13) {
        return 9U;
    }
    if (schema >= 12) {
        return 8U;
    }
    if (schema >= 11) {
        return 7U;
    }
    if (schema >= 10) {
        return 6U;
    }
    if (schema >= 8) {
        return 5U;
    }
    if (schema >= 7) {
        return 4U;
    }
    if (schema >= 5) {
        return 3U;
    }
    return 2U;
}
}  // namespace

StructuredOperatorStrategies3D::StructuredOperatorStrategies3D(
    const StructuredPdeConfig3D& config) {
    const int schema = config.schema_version;
    activation.transported_refractory = schema >= 7;
    activation.grid_local_refractory = schema >= 5;
    transport.resource_guidance = schema >= 2;
    transport.directional_sectors = schema >= 3;
    transport.cached_cohort_transport = schema >= 4;
    transport.nutrient_sectors = schema >= 5;
    transport.feasible_normal_jumps = config.normal_transport == "feasible_fixed_lattice_jump_v3";
    transport.fixed_normal_jumps = config.normal_transport != "axial_diffusion_v1";
    exchange.direction_flux = schema >= 3;
    reaction.local_density_conversion = schema >= 3;
    reaction.post_transport_conversion = schema >= 4;
    reaction.exact_growth_window = schema >= 7;
    reaction.neighbour_births = config.small_daughter_placement == "feasible_neighbor_birth_v2";
    reaction.true_normal_mean = config.growth_rate_closure == "truncated_normal_expectation_v2";
    nutrient.transient_resources = schema >= 5;
    nutrient.moving_front = schema >= 6;
    nutrient.footprint_consumers = schema >= 7;
    nutrient.persistent_workspaces = schema >= 17;
    vascular.density_exclusion = schema >= 8;
    vascular.shared_geometry = schema >= 7;
    vascular.skip_empty_cells = schema >= 9;
    vascular.removal_diagnostics = schema >= 13;
    checkpoint_version = select_checkpoint_version(schema);
}

void ActivationOperator3D::advance(StructuredPdeModel3D& state, double dt) const {
    state.refresh_activation(dt);
}

void TransportOperator3D::advance(StructuredPdeModel3D& state, double dt) const {
    state.migrate_normal_and_K(dt);
    state.migrate_active(dt);
}

void TransportOperator3D::expire(StructuredPdeModel3D& state, double dt) const {
    for (std::size_t stage = 0; stage < kStructuredStageCount3D; ++stage) {
        state.expire_active(stage, dt);
    }
}

void ExchangeOperator3D::advance(StructuredPdeModel3D& state, double dt) const {
    state.exchange_active_r_with_K(dt);
}

void ReactionOperator3D::advance(StructuredPdeModel3D& state, double dt) const {
    state.react(dt);
}

void NutrientOperator3D::advance(StructuredPdeModel3D& state, double dt) const {
    if (transient_resources) {
        state.rebuild_moving_tumour_front();
        state.advance_transient_nutrient(dt);
        state.next_nutrient_refresh_hours_ =
            state.time_hours_ + state.config_.continuum.time_step_hours;
    } else if (state.time_hours_ > state.next_nutrient_refresh_hours_ ||
               same_time(state.time_hours_, state.next_nutrient_refresh_hours_)) {
        state.solve_nutrient();
        do {
            state.next_nutrient_refresh_hours_ +=
                state.config_.continuum.nutrient.refresh_every_hours;
        } while (state.time_hours_ > state.next_nutrient_refresh_hours_ ||
                 same_time(state.time_hours_, state.next_nutrient_refresh_hours_));
    }
}

void VascularOperator3D::advance(StructuredPdeModel3D& state, double dt) const {
    state.advance_angiogenesis(dt);
}
}  // namespace atcg3d::structured_pde
