#pragma once

#include <array>
#include <cstdint>
#include <filesystem>
#include <string>

#include "config/model_config.hpp"

namespace atcg3d::continuum {

struct ContinuumGridConfig3D {
    std::array<int, 3> shape{12, 12, 12};
    std::array<double, 3> origin{-12.0, -12.0, -12.0};
    double spacing_voxels{2.0};
};

struct ContinuumMigrationConfig3D {
    std::string mapping{"fixed_26_from_base_means_v1"};
    double diffusion_scale{1.0};
    double large_mobility_multiplier{0.5};
    double activated_r_mobility_multiplier{4.0};
    double crowding_exponent{1.0};
};

struct ContinuumReactionConfig3D {
    std::string model{"abm_mean_clock_volume_filling_v1"};
    double large_daughter_vacancy_exponent{8.0};
    double small_daughter_vacancy_exponent{1.0};
    double failed_r_division_death_fraction{1.0};
    double maximum_occupied_fraction{1.0};
};

struct ContinuumNutrientConfig3D {
    std::string model{"effective_resource_surplus_v1"};
    std::string solver{"deterministic_quasi_steady_jacobi_v1"};
    double diffusion_voxels2_per_hour{1.0};
    double decay_per_hour{1.0 / 144.0};
    double vessel_exchange_per_hour{10.0};
    double vessel_value{1.0};
    double r_consumption_per_occupied_voxel_hour{0.01};
    double K_consumption_per_occupied_voxel_hour{0.01};
    double r_consumption_half_saturation{0.25};
    double K_consumption_half_saturation{0.25};
    double capacity_half_saturation{0.25};
    double maximum_capacity_multiplier{2.0};
    double refresh_every_hours{0.5};
    int solver_iterations{32};
    double relaxation{0.8};
};

struct ContinuumVascularConfig3D {
    std::string source_mode{"abm_perfusion"};
    std::string synthetic_axis{"z"};
    std::array<double, 3> synthetic_center{0.0, 0.0, 0.0};
    double synthetic_radius_voxels{1.5};
};

struct ContinuumOutputConfig3D {
    bool enabled{true};
    std::filesystem::path directory{"atcg3d_continuum_run"};
    double metrics_every_hours{1.0};
    double field_every_hours{4.0};
    double radial_profile_every_hours{4.0};
    double checkpoint_every_hours{4.0};
};

struct ContinuumModelConfig3D {
    std::string schema_name{"atcg3d.continuum_model_config"};
    int schema_version{1};
    std::string profile{"continuum_v1"};
    std::filesystem::path source_path;
    std::filesystem::path base_config_path;
    Model3DConfig base;

    std::string run_mode{"new"};
    std::filesystem::path resume_checkpoint;
    std::string initialization_mode{"base_model"};
    std::filesystem::path abm_checkpoint;
    double start_time_hours{};
    double end_time_hours{24.0};
    double time_step_hours{0.05};

    ContinuumGridConfig3D grid;
    ContinuumMigrationConfig3D migration;
    ContinuumReactionConfig3D reaction;
    ContinuumNutrientConfig3D nutrient;
    ContinuumVascularConfig3D vascular;
    ContinuumOutputConfig3D output;

    static ContinuumModelConfig3D load(const std::filesystem::path& path);
    void validate() const;
    std::uint64_t dynamics_fingerprint() const;
    std::string to_json() const;
};

}  // namespace atcg3d::continuum
