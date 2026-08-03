#pragma once

#include <cstdint>
#include <filesystem>
#include <string>

#include "config/model_config.hpp"

namespace atcg3d::nutrient {

struct NutrientFieldConfig3D {
    std::string model{"effective_resource_surplus_v1"};
    std::string solver{"deterministic_quasi_steady_jacobi_v1"};
    int block_edge{4};
    int halo_voxels{16};
    double diffusion_voxels2_per_hour{1.0};
    double decay_per_hour{1.0 / 144.0};
    double vessel_exchange_per_hour{10.0};
    double vessel_value{1.0};
    double r_consumption_per_voxel_hour{0.01};
    double K_consumption_per_voxel_hour{0.01};
    double r_consumption_half_saturation{0.25};
    double K_consumption_half_saturation{0.25};
    double capacity_half_saturation{0.25};
    double maximum_capacity_multiplier{2.0};
    double refresh_every_hours{0.25};
    int solver_iterations{128};
    double relaxation{0.8};
    double metrics_every_hours{1.0};
    double field_snapshot_every_hours{0.0};

    void validate() const;
    std::uint64_t fingerprint() const noexcept;
    std::string to_json() const;
};

struct NutrientModelConfig3D {
    std::string schema_name{"atcg3d.nutrient_model_config"};
    int schema_version{1};
    std::string profile{"nutrient_v1"};
    std::filesystem::path source_path;
    std::filesystem::path base_config_path;
    Model3DConfig base;
    NutrientFieldConfig3D nutrient;

    static NutrientModelConfig3D load(const std::filesystem::path& path);
    std::string to_json() const;
};

}  // namespace atcg3d::nutrient
