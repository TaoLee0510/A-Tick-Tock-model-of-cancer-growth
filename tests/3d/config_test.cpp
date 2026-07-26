#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include "config/model_config.hpp"

namespace {

std::string read_text(const std::filesystem::path& path) {
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("unable to read test config: " + path.string());
    return {std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};
}

std::string replace_once(std::string text,
                         const std::string& needle,
                         const std::string& replacement) {
    const std::size_t offset = text.find(needle);
    if (offset == std::string::npos) throw std::runtime_error("test fixture text was not found: " + needle);
    text.replace(offset, needle.size(), replacement);
    return text;
}

std::filesystem::path write_case(const std::string& name, const std::string& text) {
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / ("atcg3d_config_" + name + ".yaml");
    std::ofstream stream(path, std::ios::binary | std::ios::trunc);
    if (!stream) throw std::runtime_error("unable to create test config: " + path.string());
    stream << text;
    stream.close();
    return path;
}

void expect_rejected(const std::string& name, const std::string& text) {
    const std::filesystem::path path = write_case(name, text);
    bool rejected = false;
    try {
        (void)atcg3d::Model3DConfig::load(path);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    std::filesystem::remove(path);
    assert(rejected);
}

}  // namespace

int main() {
    using atcg3d::Model3DConfig;

    const std::filesystem::path smoke_path =
        std::filesystem::path("configs") / "atcg3d_smoke_test_v3.yaml";
    const std::filesystem::path production_path =
        std::filesystem::path("configs") / "atcg3d_legacy_2d_mapped_v3.yaml";
    const std::filesystem::path single_cell_v5_path =
        std::filesystem::path("configs") /
        "single_r_stage0_2160h_density_vascular_v5.yaml";
    const std::filesystem::path requested_single_cell_path =
        std::filesystem::path("configs") /
        "single_r_stage0_2160h_seed1.yaml";
    const std::string smoke_yaml = read_text(smoke_path);
    const std::string production_yaml = read_text(production_path);

    const Model3DConfig smoke = Model3DConfig::load(smoke_path);
    assert(smoke.schema_name == "atcg3d.model_config");
    assert(smoke.schema_version == 3);
    assert(smoke.profile == "smoke_test_v3");
    assert(smoke.initialization_mode == "explicit_counts");
    assert(smoke.initial_r_cells == 32 && smoke.initial_K_cells == 32);
    assert(smoke.initial_growth_rate_model == "fixed");
    assert(smoke.initial_r_growth_rate == 1.0 && smoke.initial_K_growth_rate == 1.0);
    assert(smoke.activated_r_migration_rate_model == "beta");
    assert(smoke.activated_r_migration_beta.alpha == 0.01);
    assert(smoke.activated_r_migration_beta.beta == 0.0566666667);
    assert(smoke.activated_r_migration_beta.scale == 1.0);
    assert(smoke.activated_r_migration_beta.lower_clamp_enabled);
    assert(smoke.activated_r_migration_beta.lower_clamp_threshold == 0.5);
    assert(smoke.activated_r_migration_beta.lower_clamp_value == 0.25);
    assert(smoke.initial_K_migration_rate_model == "fixed");
    assert(smoke.initial_K_migration_rate == 0.25);
    assert(smoke.normal_r_migration_beta.alpha == 5.0);
    assert(smoke.normal_r_migration_beta.beta == 5.0);
    assert(smoke.normal_r_migration_beta.scale == 0.5);
    assert(smoke.migration_activation_duration_alpha == 0.005);
    assert(smoke.migration_activation_duration_mean_fraction == 0.30);
    assert(std::abs(smoke.migration_activation_duration_beta -
                    0.011666666666666667) < 1e-16);
    assert(smoke.migration_swap_enabled);
    assert(smoke.migration_swap_stage_policy == "stage1_singleton_v1");
    assert(smoke.migration_swap_wait_fraction == 0.20);
    assert(smoke.migration_swap_post_cooldown_fraction == 0.20);
    assert(smoke.scheduler_backend == "event_queue_v1");
    assert(smoke.r_limit == 93.0);
    assert(smoke.K_limit == 108.0);
    assert(smoke.carrying_capacity_r == 186.0);
    assert(smoke.carrying_capacity_K == 216.0);
    assert(smoke.carrying_capacity_scale_2d_to_3d == 1.0);
    assert(smoke.alpha == 2.2 && smoke.beta == 0.0);
    assert(smoke.death_delay_model == "legacy_geometric_mean_v1");
    assert(smoke.death_growth_rate_threshold == 0.005);
    assert(smoke.r_death_delay_hours == 48.0);
    assert(smoke.K_death_delay_hours == 120.0);
    assert(!smoke.r_to_K_conversion.enabled);
    assert(smoke.r_to_K_conversion.density_window_edge == 70);
    assert(smoke.r_to_K_conversion.query_block_edge == 32);
    assert(smoke.r_to_K_conversion.density_threshold == 0.50);
    assert(smoke.r_to_K_conversion.probability_per_division == 0.05);
    assert(smoke.division_timing.r_max_inherent_growth_rate == 1.3171805);
    assert(smoke.division_timing.K_max_inherent_growth_rate == 0.99505180);
    assert(!smoke.ultrasmall_enabled);
    assert(std::abs(smoke.angiogenesis.seed_rate_sites_per_hour - 10.0 / 720.0) < 1e-15);
    assert(smoke.angiogenesis.roots_per_event == 1);
    assert(smoke.angiogenesis.lesion_detection_backend ==
           "sparse_coarse_blocks_v1");
    assert(smoke.angiogenesis.lesion_block_edge == 8);
    assert(smoke.angiogenesis.lesion_connectivity == 26);
    assert(smoke.angiogenesis.lesion_core_activation_occupied_fraction == 0.15);
    assert(smoke.angiogenesis.lesion_core_deactivation_occupied_fraction == 0.10);
    assert(smoke.angiogenesis.lesion_minimum_cells_per_core_block == 8);
    assert(smoke.angiogenesis.lesion_halo_blocks == 1);
    assert(smoke.angiogenesis.lesion_refresh_interval_hours == 1.0);
    assert(smoke.angiogenesis.trigger_metric ==
           "per_lesion_biological_cell_volume_voxels3");
    assert(smoke.angiogenesis.trigger_minimum_core_blocks == 4);
    assert(smoke.angiogenesis.seed_process_scope == "per_eligible_lesion");
    assert(smoke.angiogenesis.root_position_policy == "inside_surface_voxel");
    assert(smoke.angiogenesis.surface_min_local_cells == 1);
    assert(smoke.angiogenesis.max_roots_per_lesion == 64);
    assert(smoke.angiogenesis.max_active_tips_per_lesion == 128);
    assert(smoke.angiogenesis.influence_profile == "linear_cutoff");
    assert(smoke.angiogenesis.influence_activation == "immediate");
    assert(smoke.angiogenesis.outward_speed_voxels_per_hour == 2.0);
    assert(smoke.angiogenesis.inward_max_length_voxels == 128);
    assert(smoke.angiogenesis.inward_length_tortuosity_factor == 1.5);
    assert(smoke.angiogenesis.inward_exit_margin_voxels == 16.0);
    assert(smoke.angiogenesis.inward_hard_max_length_voxels == 4096);
    assert(smoke.angiogenesis.inward_far_surface_policy ==
           "continue_to_budget");
    assert(smoke.threads == 4);
    assert(smoke.parallel_mode == "adaptive_cells_and_events_v1");
    assert(smoke.parallel_min_threads == 1);
    assert(smoke.parallel_min_events_per_thread == 1024);
    assert(smoke.parallel_thread_thresholds.size() == 6);
    assert(smoke.parallel_thread_thresholds.front().minimum_cells == 0);
    assert(smoke.parallel_thread_thresholds.back().minimum_cells == 25000);
    assert(smoke.parallel_thread_thresholds.back().max_thread_fraction == 1.0);
    assert(smoke.checkpoint_format == "hdf5_base_v6_slot_journal_v8");
    assert(smoke.storage_mode == "journal_delta_hdf5_v2");
    assert(smoke.vtkhdf_compression_level == 1);
    assert(smoke.hdf5_compression_level == 1);
    assert(smoke.hdf5_chunk_elements == 262144);
    assert(smoke.preview_keyframe_every_hours == 1.0);
    assert(smoke.full_keyframe_every_hours == 24.0);
    assert(smoke.checkpoint_base_every_hours == 168.0);
    assert(smoke.checkpoint_max_delta_chain == 168);
    assert(smoke.delta_full_ratio == 0.70);
    assert(smoke.to_json().find("\"schema_version\":3") != std::string::npos);
    assert(smoke.to_json().find(
        "\"activated_r_migration_beta_scale\":1") != std::string::npos);

    const Model3DConfig defaults;
    assert(defaults.activated_r_migration_beta.scale == 1.0);
    const std::filesystem::path custom_scale_path = write_case(
        "custom_activated_r_scale",
        replace_once(smoke_yaml, "    scale: 1.0\n",
                     "    scale: 7.5\n"));
    const Model3DConfig custom_scale = Model3DConfig::load(custom_scale_path);
    std::filesystem::remove(custom_scale_path);
    assert(custom_scale.activated_r_migration_beta.scale == 7.5);

    const Model3DConfig production = Model3DConfig::load(production_path);
    assert(production.profile == "legacy_2d_mapped_v3");
    assert(production.initialization_mode == "legacy_geometry_fill_v1");
    assert(production.initial_r_cells == 0 && production.initial_K_cells == 0);
    assert(production.initial_radius == 60);
    assert(production.initial_shell_inner_radius == 50);
    assert(production.initial_inner_small_radius == 55);
    assert(production.initial_r_fraction == 0.50);
    assert(production.initial_growth_rate_model == "legacy_truncated_normal_v1");
    assert(production.initial_r_growth_truncated_normal.mean == 1.1832);
    assert(production.initial_r_growth_truncated_normal.standard_deviation == 0.2441);
    assert(production.initial_r_growth_truncated_normal.minimum == 1.0722619);
    assert(production.initial_r_growth_truncated_normal.maximum == 1.3171805);
    assert(production.initial_K_growth_truncated_normal.mean == 0.6832);
    assert(production.initial_K_growth_truncated_normal.standard_deviation == 0.3764);
    assert(production.initial_K_growth_truncated_normal.minimum == 0.33963482);
    assert(production.initial_K_growth_truncated_normal.maximum == 0.99505180);
    assert(production.activated_r_migration_rate_model == "beta");
    assert(production.activated_r_migration_beta.alpha == 0.01);
    assert(production.activated_r_migration_beta.beta == 0.0566666667);
    assert(production.activated_r_migration_beta.scale == 1.0);
    assert(production.activated_r_migration_beta.lower_clamp_enabled);
    assert(production.activated_r_migration_beta.lower_clamp_threshold == 0.5);
    assert(production.activated_r_migration_beta.lower_clamp_value == 0.25);
    assert(production.initial_K_migration_rate_model == "legacy_beta_v1");
    assert(production.initial_K_migration_beta.alpha == 5.0);
    assert(production.initial_K_migration_beta.beta == 5.0);
    assert(production.initial_K_migration_beta.scale == 0.25);
    assert(!production.initial_K_migration_beta.lower_clamp_enabled);
    assert(production.r_to_K_conversion.enabled);
    assert(!production.migration_swap_enabled);
    assert(production.angiogenesis.enabled);
    assert(production.angiogenesis.outward_speed_voxels_per_hour == 2.0);
    Model3DConfig slow_outward_vessel = production;
    slow_outward_vessel.angiogenesis.outward_speed_voxels_per_hour =
        std::sqrt(3.0);
    bool rejected_slow_outward_vessel = false;
    try {
        slow_outward_vessel.validate();
    } catch (const std::invalid_argument&) {
        rejected_slow_outward_vessel = true;
    }
    assert(rejected_slow_outward_vessel);
    Model3DConfig slower_than_inward_vessel = production;
    slower_than_inward_vessel.activated_r_migration_beta.scale = 0.01;
    slower_than_inward_vessel.angiogenesis.outward_speed_voxels_per_hour =
        slower_than_inward_vessel.angiogenesis.inward_speed_voxels_per_hour;
    bool rejected_slower_than_inward_vessel = false;
    try {
        slower_than_inward_vessel.validate();
    } catch (const std::invalid_argument&) {
        rejected_slower_than_inward_vessel = true;
    }
    assert(rejected_slower_than_inward_vessel);
    Model3DConfig faster_outward_vessel = production;
    faster_outward_vessel.angiogenesis.outward_speed_voxels_per_hour = 2.0;
    faster_outward_vessel.validate();
    assert(production.threads == 8);
    assert(production.to_json().find(
        "\"initial_growth_rate_model\":\"legacy_truncated_normal_v1\"") !=
        std::string::npos);
    assert(production.to_json().find(
        "\"initial_K_migration_rate_model\":\"legacy_beta_v1\"") !=
        std::string::npos);
    Model3DConfig changed_rate = production;
    changed_rate.activated_r_migration_beta.scale = 2.0;
    assert(changed_rate.dynamics_json() != production.dynamics_json());
    Model3DConfig changed_storage = production;
    changed_storage.vtkhdf_compression_level = 9;
    changed_storage.hdf5_compression_level = 0;
    changed_storage.full_keyframe_every_hours = 48.0;
    changed_storage.checkpoint_base_every_hours = 336.0;
    changed_storage.checkpoint_max_delta_chain = 336;
    changed_storage.delta_full_ratio = 0.5;
    assert(changed_storage.dynamics_json() == production.dynamics_json());

    const Model3DConfig single_cell_v5 =
        Model3DConfig::load(single_cell_v5_path);
    assert(single_cell_v5.profile ==
           "single_r_stage0_2160h_density_vascular_v5");
    assert(single_cell_v5.activated_r_migration_beta.scale == 3.0);
    assert(single_cell_v5.angiogenesis.outward_speed_voxels_per_hour == 6.0);
    assert(single_cell_v5.angiogenesis.outward_speed_voxels_per_hour >
           std::sqrt(3.0) *
               single_cell_v5.activated_r_migration_beta.scale);
    assert(single_cell_v5.threads == 18);
    assert(single_cell_v5.output_directory ==
           "/Volumes/Work_Active/simulation/ver7/"
           "run_single_r_stage0_2160h_seed1_density_vascular_v5");

    const Model3DConfig requested_single_cell =
        Model3DConfig::load(requested_single_cell_path);
    assert(requested_single_cell.profile == "single_r_stage0_2160h_seed1");
    assert(requested_single_cell.activated_r_migration_beta.scale == 3.0);
    assert(requested_single_cell.angiogenesis.outward_speed_voxels_per_hour ==
           6.0);
    assert(std::abs(
               requested_single_cell.angiogenesis.seed_volume_exponent -
               2.0 / 3.0) < 1e-9);
    assert(requested_single_cell.angiogenesis.seed_minimum_rate_multiplier ==
           0.25);
    assert(requested_single_cell.angiogenesis.inward_far_surface_policy ==
           "continue_to_budget");
    Model3DConfig invalid_inward_budget = requested_single_cell;
    invalid_inward_budget.angiogenesis.inward_hard_max_length_voxels =
        invalid_inward_budget.angiogenesis.inward_max_length_voxels - 1;
    bool rejected_invalid_inward_budget = false;
    try {
        invalid_inward_budget.validate();
    } catch (const std::invalid_argument&) {
        rejected_invalid_inward_budget = true;
    }
    assert(rejected_invalid_inward_budget);
    assert(requested_single_cell.threads == 18);
    assert(requested_single_cell.migration_swap_enabled);
    assert(requested_single_cell.migration_swap_wait_fraction == 0.20);
    assert(requested_single_cell.migration_swap_post_cooldown_fraction ==
           0.20);
    assert(requested_single_cell.scheduler_backend ==
           "deterministic_exact_window_v3");
    assert(requested_single_cell.output_directory ==
           "/Volumes/Work_Active/simulation/ver7/"
           "run_single_r_stage0_2160h_seed1");

    Model3DConfig changed_conversion = production;
    changed_conversion.r_to_K_conversion.probability_per_division += 0.01;
    assert(changed_conversion.dynamics_json() != production.dynamics_json());

    expect_rejected("unknown_root",
                    replace_once(smoke_yaml, "profile: smoke_test_v3\n",
                                 "profile: smoke_test_v3\nunknown_root: 1\n"));
    expect_rejected("schema_v2",
                    replace_once(smoke_yaml, "  version: 3\n",
                                 "  version: 2\n"));
    expect_rejected("missing_activated_r_rate",
                    replace_once(
                        smoke_yaml,
                        "  activated_r_rate:\n    model: beta\n    alpha: 0.01\n    beta: 0.0566666667\n    scale: 1.0\n    lower_clamp:\n      enabled: true\n      threshold: 0.5\n      value: 0.25\n",
                        ""));
    expect_rejected("legacy_initial_r_field",
                    replace_once(smoke_yaml,
                                 "  migration_rate:\n    K:\n",
                                 "  migration_rate:\n    r: 0.25\n    K:\n"));
    expect_rejected("unknown_nested",
                    replace_once(smoke_yaml, "  thin_layer: false\n",
                                 "  thin_layer: false\n  typo_field: 1\n"));
    expect_rejected("duplicate",
                    replace_once(smoke_yaml, "profile: smoke_test_v3\n",
                                 "profile: smoke_test_v3\nprofile: duplicate\n"));
    expect_rejected("missing",
                    replace_once(smoke_yaml, "  threads: 4\n", ""));
    expect_rejected("loose_bool",
                    replace_once(smoke_yaml, "  thin_layer: false\n",
                                 "  thin_layer: yes\n"));
    expect_rejected("quoted_bool",
                    replace_once(smoke_yaml, "  thin_layer: false\n",
                                 "  thin_layer: \"false\"\n"));
    expect_rejected("nan",
                    replace_once(smoke_yaml, "  end_time_hours: 24.0\n",
                                 "  end_time_hours: .nan\n"));
    expect_rejected("infinity",
                    replace_once(smoke_yaml, "  end_time_hours: 24.0\n",
                                 "  end_time_hours: .inf\n"));
    expect_rejected("negative_time",
                    replace_once(smoke_yaml, "  end_time_hours: 24.0\n",
                                 "  end_time_hours: -1.0\n"));
    expect_rejected("wrong_vector_length",
                    replace_once(smoke_yaml,
                                 "  minimum: [-1000000, -1000000, -1000000]\n",
                                 "  minimum: [0, 0]\n"));
    expect_rejected("integer_overflow",
                    replace_once(smoke_yaml, "  threads: 4\n",
                                 "  threads: 999999999999999999999999\n"));
    expect_rejected("parallel_zero_events",
                    replace_once(smoke_yaml, "  min_events_per_thread: 1024\n",
                                 "  min_events_per_thread: 0\n"));
    expect_rejected(
        "zero_swap_wait",
        replace_once(smoke_yaml, "    wait_fraction: 0.20\n",
                     "    wait_fraction: 0.0\n"));
    expect_rejected(
        "unsupported_swap_stage",
        replace_once(smoke_yaml,
                     "    stage_policy: stage1_singleton_v1\n",
                     "    stage_policy: all_stages\n"));
    expect_rejected("parallel_unsorted_thresholds",
                    replace_once(smoke_yaml, "    - minimum_cells: 5000\n",
                                 "    - minimum_cells: 0\n"));
    expect_rejected("parallel_decreasing_fraction",
                    replace_once(smoke_yaml, "      max_thread_fraction: 0.60\n",
                                 "      max_thread_fraction: 0.40\n"));
    expect_rejected("delayed_perfusion",
                    replace_once(smoke_yaml, "    activation: immediate\n",
                                 "    activation: after_outward_connection\n"));
    expect_rejected("multiple_documents", smoke_yaml + "\n---\nprofile: second\n");
    expect_rejected("poisson_roots",
                    replace_once(smoke_yaml, "    roots_per_event: 1\n",
                                 "    roots_per_event: 2\n"));
    expect_rejected("invalid_lesion_hysteresis",
                    replace_once(
                        smoke_yaml,
                        "    core_deactivation_occupied_fraction: 0.10\n",
                        "    core_deactivation_occupied_fraction: 0.20\n"));
    expect_rejected("invalid_lesion_connectivity",
                    replace_once(smoke_yaml, "    connectivity: 26\n",
                                 "    connectivity: 18\n"));
    expect_rejected("overflowing_lesion_block",
                    replace_once(smoke_yaml, "    block_edge: 8\n",
                                 "    block_edge: 2147483647\n"));
    expect_rejected("excessive_lesion_halo",
                    replace_once(smoke_yaml, "    halo_blocks: 1\n",
                                 "    halo_blocks: 1025\n"));
    expect_rejected("zero_lesion_refresh_interval",
                    replace_once(smoke_yaml,
                                 "    refresh_interval_hours: 1.0\n",
                                 "    refresh_interval_hours: 0.0\n"));
    expect_rejected("unsupported_seed_scope",
                    replace_once(smoke_yaml,
                                 "    scope: per_eligible_lesion\n",
                                 "    scope: global_tumour\n"));
    expect_rejected("outside_root_policy",
                    replace_once(smoke_yaml,
                                 "    root_position_policy: inside_surface_voxel\n",
                                 "    root_position_policy: outside_surface_voxel\n"));
    expect_rejected("checkpoint_v1",
                    replace_once(smoke_yaml,
                                 "  checkpoint_format: hdf5_base_v6_slot_journal_v8\n",
                                 "  checkpoint_format: hdf5_v1\n"));
    expect_rejected(
        "invalid_hdf5_compression",
        replace_once(smoke_yaml, "    hdf5_compression_level: 1\n",
                     "    hdf5_compression_level: 10\n"));
    expect_rejected(
        "invalid_delta_ratio",
        replace_once(smoke_yaml, "    delta_full_ratio: 0.70\n",
                     "    delta_full_ratio: 0.0\n"));
    expect_rejected(
        "keyframe_before_frame",
        replace_once(smoke_yaml,
                     "    preview_keyframe_every_hours: 1.0\n",
                     "    preview_keyframe_every_hours: 0.5\n"));
    expect_rejected(
        "unknown_storage_key",
        replace_once(smoke_yaml, "    delta_full_ratio: 0.70\n",
                     "    delta_full_ratio: 0.70\n    typo: 1\n"));
    expect_rejected("complete_density_relief",
                    replace_once(smoke_yaml, "    max_relief_fraction: 0.50\n",
                                 "    max_relief_fraction: 1.0\n"));
    expect_rejected("unsupported_growth_model",
                    replace_once(smoke_yaml, "  growth_rate:\n    model: fixed\n",
                                 "  growth_rate:\n    model: gaussian\n"));
    expect_rejected("fixed_growth_extra_block",
                    replace_once(smoke_yaml, "    fixed:\n      r: 1.0\n      K: 1.0\n",
                                 "    fixed:\n      r: 1.0\n      K: 1.0\n    typo: 1\n"));
    expect_rejected("invalid_growth_sd",
                    replace_once(production_yaml, "        standard_deviation: 0.2441\n",
                                 "        standard_deviation: 0.0\n"));
    expect_rejected("invalid_beta_alpha",
                    replace_once(production_yaml, "    alpha: 0.01\n",
                                 "    alpha: 0.0\n"));
    expect_rejected("invalid_activated_r_scale",
                    replace_once(smoke_yaml, "    scale: 1.0\n",
                                 "    scale: 0.0\n"));
    expect_rejected("unsupported_activated_r_model",
                    replace_once(smoke_yaml,
                                 "  activated_r_rate:\n    model: beta\n",
                                 "  activated_r_rate:\n    model: fixed\n"));
    expect_rejected("disabled_clamp_extra_value",
                    replace_once(production_yaml, "        lower_clamp:\n          enabled: false\n",
                                 "        lower_clamp:\n          enabled: false\n          value: 0.1\n"));
    expect_rejected("invalid_conversion_probability",
                    replace_once(smoke_yaml,
                                 "    probability_per_division: 0.05\n",
                                 "    probability_per_division: 1.01\n"));
    expect_rejected("invalid_conversion_query_block",
                    replace_once(smoke_yaml,
                                 "    query_block_edge: 32\n",
                                 "    query_block_edge: 71\n"));
    expect_rejected("unknown_conversion_key",
                    replace_once(smoke_yaml,
                                 "    probability_per_division: 0.05\n",
                                 "    probability_per_division: 0.05\n    typo: 1\n"));
    expect_rejected("invalid_r_growth_cap",
                    replace_once(smoke_yaml,
                                 "    r_max_inherent_growth_rate: 1.3171805\n",
                                 "    r_max_inherent_growth_rate: 0.0\n"));
    expect_rejected("contradictory_activation_beta_mean",
                    replace_once(smoke_yaml,
                                 "    beta: 0.011666666666666667\n",
                                 "    beta: 0.02\n"));
    expect_rejected("unknown_normal_migration_key",
                    replace_once(smoke_yaml,
                                 "    scale: 0.5\n",
                                 "    scale: 0.5\n    typo: 1\n"));
    expect_rejected("unreachable_outward_connection",
                    replace_once(smoke_yaml,
                                 "    external_connection_distance_voxels: 64.0\n",
                                 "    external_connection_distance_voxels: 129.0\n"));

    const std::string hash_path_yaml = replace_once(
        smoke_yaml, "  directory: atcg3d_smoke_run\n",
        "  directory: \"run # with colon: ok\"\n");
    const std::filesystem::path hash_path = write_case("quoted_path", hash_path_yaml);
    const Model3DConfig quoted_path_config = Model3DConfig::load(hash_path);
    assert(quoted_path_config.output_directory == "run # with colon: ok");
    std::filesystem::remove(hash_path);

    return 0;
}
