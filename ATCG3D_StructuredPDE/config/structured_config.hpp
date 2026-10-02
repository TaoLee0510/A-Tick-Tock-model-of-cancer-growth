#pragma once

#include <cstdint>
#include <filesystem>
#include <string>

#include "config/continuum_config.hpp"

namespace atcg3d::structured_pde {

struct StructuredMigrationConfig3D {
    std::string model{"abm_activation_clock_discrete_velocity_v1"};
    std::string activation_density{"abm_anchor_box_v1"};
    std::string activation_clock{"beta_mean_remaining_cycle_v1"};
    std::string activation_stop{"clock_expiry_v1"};
    std::string direction_transport{"fixed_direction_jump_v1"};
    int direction_density_window_edge{};
    int direction_nutrient_window_edge{};
    double chemotaxis_strength{};
    double zero_gradient_tolerance{1.0e-9};
    double reactivation_cooldown_hours{};
    double reactivation_density_threshold{};
    std::string crowding_exchange{"none"};
    bool vessel_exclusion{};
    double maximum_move_probability_per_substep{0.15};
    double minimum_density{1.0e-12};
};

struct StructuredPdeConfig3D {
    int schema_version{1};
    std::string profile;
    std::string storage_model{"dense_v1"};
    std::uint64_t maximum_active_voxels{1000000};
    std::filesystem::path source_path;
    std::filesystem::path continuum_config_path;
    continuum::ContinuumModelConfig3D continuum;
    StructuredMigrationConfig3D migration;

    static StructuredPdeConfig3D load(const std::filesystem::path& path);
    void validate() const;
    std::uint64_t dynamics_fingerprint() const;
    std::string to_json() const;
};

}  // namespace atcg3d::structured_pde
