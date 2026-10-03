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
    double activation_time_bin_width_hours{0.5};
    double activation_maximum_hours{32.0};
    std::string activation_rate_model{"phenotype_mean_v1"};
    int activation_rate_bins{8};
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
    std::string division_clock_model{"mean_rate_v1"};
    double division_work_bin_width{0.5};
    double division_maximum_work{128.0};
    std::string operator_model{"published_operators_v1"};
    bool migration_operator_enabled{true};
    bool activation_operator_enabled{true};
    bool division_operator_enabled{true};
    bool exchange_operator_enabled{true};
    std::string growth_rate_closure{"clipped_location_mean_v1"};
    std::string normal_transport{"axial_diffusion_v1"};
    std::string sector_mean_model{"published_sector_sums_v1"};
    std::string small_daughter_placement{"local_growth_v1"};
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
