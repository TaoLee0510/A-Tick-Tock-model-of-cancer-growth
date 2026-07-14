#pragma once

#include <cstdint>
#include <filesystem>
#include <string>

#include "core/types.hpp"

namespace atcg3d {

struct Model3DConfig {
    std::uint32_t schema_version{1};
    std::string profile{"legacy_like_v1"};
    std::uint64_t seed{1};

    std::string domain_policy{"expandable_sparse"};
    int chunk_edge{32};
    bool bounded_domain{false};
    Vec3i domain_min{-1000000, -1000000, -1000000};
    Vec3i domain_max{1000000, 1000000, 1000000};
    bool thin_layer{false};

    std::string direction_set{"fixed_26_v1"};
    double continue_probability{0.90};
    double turn_half_angle_degrees{45.0};
    int direction_density_radius{5};
    double direction_density_half_angle_degrees{45.0};
    double direction_density_threshold{0.60};
    bool persistence_uses_density{false};
    double distance_weight_exponent{0.0};
    bool migration_activation_enabled{true};
    int migration_activation_window_edge{70};
    int migration_activation_block_edge{32};
    double migration_activation_threshold{0.90};

    std::string density_backend{"block_anchor_v1"};
    int density_block_edge{4};
    int growth_density_window_edge{6};
    double carrying_capacity_scale_2d_to_3d{6.0};
    double r_limit{24.0};
    double K_limit{24.0};
    double carrying_capacity_r{186.0};
    double carrying_capacity_K{216.0};
    double alpha{0.0};
    double beta{0.0};

    int large_footprint_edge{2};
    int small_footprint_voxels{1};
    int division_shell_radius{2};
    bool allow_shape_reduction{true};
    bool ultrasmall_enabled{true};

    std::uint64_t initial_r_cells{32};
    std::uint64_t initial_K_cells{32};
    int initial_radius{8};
    int initial_shell_thickness{2};
    double initial_r_growth_rate{1.0};
    double initial_K_growth_rate{1.0};
    double initial_r_migration_rate{0.25};
    double initial_K_migration_rate{0.25};
    double r_death_delay_hours{72.0};
    double K_death_delay_hours{72.0};

    double end_time_hours{24.0};
    std::uint64_t max_events{1000000};
    int threads{1};
    std::string scheduler_backend{"event_queue_v1"};
    double conflict_bucket_hours{0.0};

    bool output_enabled{false};
    std::string preview_mode{"stable_uid_hash_v1"};
    std::string full_format{"vtkhdf_points_v1"};
    std::string checkpoint_format{"hdf5_v1"};
    std::filesystem::path output_directory{"atcg3d_run"};
    double preview_every_hours{1.0};
    double full_every_hours{6.0};
    double checkpoint_every_hours{24.0};
    std::uint64_t preview_max_cells{1000000};
    std::uint64_t preview_seed{0};

    static Model3DConfig load(const std::filesystem::path& path);
    void apply_override(const std::string& assignment);
    void validate() const;
    std::string to_json() const;
    std::string dynamics_json() const;
};

}  // namespace atcg3d
