#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "core/types.hpp"

namespace atcg3d {

struct Legacy2DMappingConfig {
    std::string source_profile{"legacy_2d_default_no_ultrasmall"};
    std::string conversion{"density_window_volume_ratio_v1"};
    double density_count_scale_2d_to_3d{6.0};
    double source_r_limit{15.5};
    double source_K_limit{18.0};
    double source_carrying_capacity_r{31.0};
    double source_carrying_capacity_K{36.0};
    int source_outer_radius{60};
    int source_inner_radius{50};
    bool scaled_down_initialization{true};
};

struct DivisionTimingConfig {
    double base_cycle_hours{24.0};
    double minimum_fraction{0.90};
    double stochastic_tail_fraction{0.10};
    double stochastic_time_quantum_hours{1.0};
    double retry_delay_hours{1.0};
    double inherited_growth_multiplier_min{0.95};
    double inherited_growth_multiplier_max{1.05};
    double r_max_inherent_growth_rate{1.3171805};
    double K_max_inherent_growth_rate{0.99505180};
};

struct RToKConversionConfig {
    bool enabled{false};
    int density_window_edge{70};
    int query_block_edge{32};
    double density_threshold{0.50};
    double probability_per_division{0.05};
};

struct DisplayRadiusConfig {
    float large{1.0F};
    float small{0.5F};
    float ultrasmall{0.25F};
};

struct ParallelThreadThresholdConfig {
    std::uint64_t minimum_cells{};
    double max_thread_fraction{1.0};
};

struct TruncatedNormalRateConfig {
    double mean{};
    double standard_deviation{1.0};
    double minimum{};
    double maximum{1.0};
};

struct BetaRateConfig {
    double alpha{1.0};
    double beta{1.0};
    double scale{1.0};
    bool lower_clamp_enabled{false};
    double lower_clamp_threshold{};
    double lower_clamp_value{};
};

struct AngiogenesisConfig {
    bool enabled{false};

    std::string lesion_detection_backend{"sparse_coarse_blocks_v1"};
    int lesion_block_edge{8};
    int lesion_connectivity{26};
    double lesion_core_activation_occupied_fraction{0.15};
    double lesion_core_deactivation_occupied_fraction{0.10};
    std::uint64_t lesion_minimum_cells_per_core_block{8};
    double lesion_minimum_biological_volume_per_core_block{0.0};
    int lesion_halo_blocks{1};
    double lesion_refresh_interval_hours{1.0};

    std::string trigger_metric{"per_lesion_biological_cell_volume_voxels3"};
    double trigger_activation_volume_voxels3{100000.0};
    double trigger_deactivation_volume_voxels3{80000.0};
    double trigger_delay_hours{0.0};
    std::uint64_t trigger_minimum_core_blocks{4};
    double stage0_biological_volume_voxels3{8.0};
    double stage1_biological_volume_voxels3{1.0};
    double stage2_biological_volume_voxels3{0.5};

    std::string seed_process_model{"homogeneous_poisson"};
    std::string seed_process_scope{"per_eligible_lesion"};
    double seed_rate_sites_per_30_days{10.0};
    double seed_rate_sites_per_hour{10.0 / 720.0};
    // Schema-v3 invariant: one Poisson site arrival attempts at most one root.
    // Multiple roots arise from multiple arrivals, so this must remain 1.
    std::uint32_t roots_per_event{1};
    int surface_min_separation_voxels{8};
    std::uint32_t surface_max_sampling_attempts{4096};
    std::uint64_t max_total_roots{64};
    std::uint64_t max_active_tips{128};
    std::uint64_t max_roots_per_lesion{64};
    std::uint64_t max_active_tips_per_lesion{128};
    std::string root_position_policy{"inside_surface_voxel"};
    std::uint64_t surface_min_local_cells{1};

    double diameter_voxels{3.0};
    double inward_speed_voxels_per_hour{0.50};
    // Must exceed the maximum possible activated-r 3D path speed.  Cell
    // migration rates are moves/hour and one fixed-26 step can span sqrt(3)
    // voxels, whereas vessel speed is already expressed in voxels/hour.
    double outward_speed_voxels_per_hour{2.0};
    int inward_max_length_voxels{128};
    int outward_max_length_voxels{128};
    double inward_target_tolerance_voxels{2.0};
    double outward_external_connection_distance_voxels{64.0};

    std::string direction_model{"forward_biased_26_v1"};
    double direction_forward_bias{1.0};
    double direction_half_angle_degrees{45.0};
    double direction_turn_half_angle_degrees{45.0};
    double direction_persistence_probability{0.90};
    double direction_distance_weight_exponent{0.0};

    std::string inward_replacement_policy{"remove_whole_cells"};
    std::string outward_occupancy_policy{"permanent_reserved_placeholder"};
    std::string vessel_collision_policy{"anastomose_and_stop"};
    std::string boundary_policy{"stop"};

    std::string influence_profile{"linear_cutoff"};
    double influence_max_relief_fraction{0.50};
    double influence_decay_length_voxels{4.0};
    double influence_cutoff_radius_voxels{12.0};
    std::string influence_scope{"growth_density"};
    // Every generated vessel is perfused immediately. The field remains in
    // the versioned schema so older run metadata remains self-describing.
    std::string influence_activation{"immediate"};
};

struct Model3DConfig {
    std::string schema_name{"atcg3d.model_config"};
    std::uint32_t schema_version{3};
    std::string profile{"legacy_2d_mapped_v3"};
    Legacy2DMappingConfig legacy_mapping{};
    std::uint64_t seed{1};

    std::string run_mode{"new"};
    std::filesystem::path resume_checkpoint{};

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
    BetaRateConfig normal_r_migration_beta{
        5.0, 5.0, 0.5, false, 0.0, 0.0};
    // Inherent r-cell migration rate used while density activation is active.
    // This is deliberately independent of initialization: initial r cells and
    // each newly committed r division cycle draw from the same configured law.
    std::string activated_r_migration_rate_model{"beta"};
    BetaRateConfig activated_r_migration_beta{
        0.01, 0.0566666667, 1.0, true, 0.5, 0.25};
    double migration_activation_duration_alpha{0.005};
    double migration_activation_duration_mean_fraction{0.30};
    double migration_activation_duration_beta{0.011666666666666667};

    std::string density_backend{"block_anchor_v1"};
    int density_block_edge{4};
    int growth_density_window_edge{6};
    // The YAML loader performs the 2D-to-3D conversion once. This multiplier
    // remains one so hot-loop biology code never applies the mapping twice.
    double carrying_capacity_scale_2d_to_3d{1.0};
    double r_limit{93.0};
    double K_limit{108.0};
    double carrying_capacity_r{186.0};
    double carrying_capacity_K{216.0};
    double alpha{2.2};
    double beta{0.0};

    int large_footprint_edge{2};
    int small_footprint_voxels{1};
    int division_shell_radius{2};
    bool allow_shape_reduction{true};
    bool ultrasmall_enabled{true};

    std::string initialization_mode{"explicit_counts"};
    std::uint64_t initial_r_cells{32};
    std::uint64_t initial_K_cells{32};
    int initial_radius{8};
    int initial_shell_thickness{2};
    int initial_shell_inner_radius{6};
    int initial_inner_small_radius{6};
    double initial_r_fraction{0.50};
    std::string initial_growth_rate_model{"fixed"};
    double initial_r_growth_rate{1.0};
    double initial_K_growth_rate{1.0};
    TruncatedNormalRateConfig initial_r_growth_truncated_normal{
        1.1832, 0.2441, 1.0722619, 1.3171805};
    TruncatedNormalRateConfig initial_K_growth_truncated_normal{
        0.6832, 0.3764, 0.33963482, 0.99505180};
    // Initial/cycle K-cell migration is configured separately because K cells
    // do not use the r-cell density-activation state machine.
    std::string initial_K_migration_rate_model{"fixed"};
    double initial_K_migration_rate{0.25};
    BetaRateConfig initial_K_migration_beta{
        5.0, 5.0, 0.25, false, 0.0, 0.0};
    std::string death_delay_model{"legacy_geometric_mean_v1"};
    double death_growth_rate_threshold{0.005};
    double r_death_delay_hours{48.0};
    double K_death_delay_hours{120.0};
    RToKConversionConfig r_to_K_conversion{};
    double initial_large_fraction{0.50};
    DivisionTimingConfig division_timing{};

    double end_time_hours{24.0};
    std::uint64_t max_events{1000000};
    // Maximum worker count. The adaptive policy may select fewer workers for
    // small populations or small same-time event batches.
    int threads{1};
    std::string parallel_mode{"adaptive_cells_and_events_v1"};
    int parallel_min_threads{1};
    std::uint64_t parallel_min_events_per_thread{1024};
    std::vector<ParallelThreadThresholdConfig> parallel_thread_thresholds{
        {0, 0.50}, {5000, 0.60}, {10000, 0.70},
        {15000, 0.80}, {20000, 0.90}, {25000, 1.00}};
    std::string scheduler_backend{"event_queue_v1"};
    double conflict_bucket_hours{0.0};

    bool output_enabled{false};
    std::string preview_mode{"stable_uid_hash_v1"};
    std::string full_format{"vtkhdf_points_v1"};
    std::string checkpoint_format{"hdf5_v3"};
    std::filesystem::path output_directory{"atcg3d_run"};
    double preview_every_hours{1.0};
    double full_every_hours{6.0};
    double checkpoint_every_hours{24.0};
    std::uint64_t preview_max_cells{1000000};
    std::uint64_t preview_seed{0};
    DisplayRadiusConfig display_radius{};

    AngiogenesisConfig angiogenesis{};

    static Model3DConfig load(const std::filesystem::path& path);
    void validate() const;
    std::string to_json() const;
    std::string dynamics_json() const;
};

}  // namespace atcg3d
