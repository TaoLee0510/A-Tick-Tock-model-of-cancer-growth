#pragma once

#include <cstdint>

#include "config/structured_config.hpp"

namespace atcg3d::structured_pde {
class StructuredPdeModel3D;

// These values describe arithmetic contracts, rather than checkpoint formats.
// A newer schema composes existing contracts until it explicitly selects a
// new rule. Selection happens once, before any state is initialized.
struct ActivationOperator3D {
    bool transported_refractory{};
    bool grid_local_refractory{};
    void advance(StructuredPdeModel3D& state, double dt) const;
};

struct TransportOperator3D {
    bool resource_guidance{};
    bool directional_sectors{};
    bool cached_cohort_transport{};
    bool nutrient_sectors{};
    bool feasible_normal_jumps{};
    bool fixed_normal_jumps{};
    void advance(StructuredPdeModel3D& state, double dt) const;
    void expire(StructuredPdeModel3D& state, double dt) const;
};

struct ExchangeOperator3D {
    bool direction_flux{};
    void advance(StructuredPdeModel3D& state, double dt) const;
};

struct ReactionOperator3D {
    bool local_density_conversion{};
    bool post_transport_conversion{};
    bool exact_growth_window{};
    bool neighbour_births{};
    bool true_normal_mean{};
    void advance(StructuredPdeModel3D& state, double dt) const;
};

struct NutrientOperator3D {
    bool transient_resources{};
    bool moving_front{};
    bool footprint_consumers{};
    bool persistent_workspaces{};
    void advance(StructuredPdeModel3D& state, double dt) const;
};

struct VascularOperator3D {
    bool density_exclusion{};
    bool shared_geometry{};
    bool skip_empty_cells{};
    bool removal_diagnostics{};
    void advance(StructuredPdeModel3D& state, double dt) const;
};

struct StructuredOperatorStrategies3D {
    explicit StructuredOperatorStrategies3D(const StructuredPdeConfig3D& config);
    ActivationOperator3D activation;
    TransportOperator3D transport;
    ExchangeOperator3D exchange;
    ReactionOperator3D reaction;
    NutrientOperator3D nutrient;
    VascularOperator3D vascular;
    std::uint32_t checkpoint_version{};
};
}  // namespace atcg3d::structured_pde
