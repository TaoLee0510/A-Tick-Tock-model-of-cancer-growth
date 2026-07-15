#include "config/model_config.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace atcg3d {
namespace {

bool approximately_equal(double lhs, double rhs) {
    return std::abs(lhs - rhs) <=
        1e-12 * std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

void require_finite(double value, const char* path) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(path) + " must be finite");
    }
}

void require_probability(double value, const char* path) {
    require_finite(value, path);
    if (value < 0.0 || value > 1.0) {
        throw std::invalid_argument(std::string(path) + " must be between 0 and 1");
    }
}

void validate_truncated_normal(const TruncatedNormalRateConfig& parameters,
                               const char* path) {
    require_finite(parameters.mean, path);
    require_finite(parameters.standard_deviation, path);
    require_finite(parameters.minimum, path);
    require_finite(parameters.maximum, path);
    if (parameters.standard_deviation <= 0.0 ||
        parameters.minimum >= parameters.maximum ||
        parameters.mean < parameters.minimum || parameters.mean > parameters.maximum) {
        throw std::invalid_argument(
            std::string(path) + " must have sd>0 and minimum<=mean<=maximum");
    }
    constexpr double inverse_sqrt_two = 0.70710678118654752440;
    const double lower_z = (parameters.minimum - parameters.mean) /
                           parameters.standard_deviation;
    const double upper_z = (parameters.maximum - parameters.mean) /
                           parameters.standard_deviation;
    const double probability_span =
        0.5 * std::erfc(-upper_z * inverse_sqrt_two) -
        0.5 * std::erfc(-lower_z * inverse_sqrt_two);
    if (!(probability_span > 0.0) || !std::isfinite(probability_span)) {
        throw std::invalid_argument(std::string(path) + " has no numerically resolvable mass");
    }
}

void validate_beta_rate(const BetaRateConfig& parameters, const char* path) {
    require_finite(parameters.alpha, path);
    require_finite(parameters.beta, path);
    require_finite(parameters.scale, path);
    require_finite(parameters.lower_clamp_threshold, path);
    require_finite(parameters.lower_clamp_value, path);
    if (parameters.alpha <= 0.0 || parameters.beta <= 0.0 || parameters.scale <= 0.0) {
        throw std::invalid_argument(std::string(path) + " alpha, beta, and scale must be positive");
    }
    if (parameters.lower_clamp_enabled &&
        (parameters.lower_clamp_threshold < 0.0 ||
         parameters.lower_clamp_threshold > parameters.scale ||
         parameters.lower_clamp_value < 0.0 ||
         parameters.lower_clamp_value > parameters.lower_clamp_threshold)) {
        throw std::invalid_argument(std::string(path) + " lower clamp is invalid");
    }
}

std::string json_escape(const std::string& value) {
    std::ostringstream escaped;
    escaped << std::hex << std::setfill('0');
    for (const unsigned char character : value) {
        switch (character) {
            case '"': escaped << "\\\""; break;
            case '\\': escaped << "\\\\"; break;
            case '\b': escaped << "\\b"; break;
            case '\f': escaped << "\\f"; break;
            case '\n': escaped << "\\n"; break;
            case '\r': escaped << "\\r"; break;
            case '\t': escaped << "\\t"; break;
            default:
                if (character < 0x20U) {
                    escaped << "\\u" << std::setw(4) << static_cast<unsigned int>(character);
                } else {
                    escaped << static_cast<char>(character);
                }
        }
    }
    return escaped.str();
}

const char* json_bool(bool value) {
    return value ? "true" : "false";
}

}  // namespace

void Model3DConfig::validate() const {
    if (schema_name != "atcg3d.model_config" || schema_version != 2) {
        throw std::invalid_argument("unsupported model configuration schema");
    }
    if (profile.empty()) throw std::invalid_argument("profile must not be empty");
    if (legacy_mapping.source_profile.empty() ||
        legacy_mapping.conversion != "density_window_volume_ratio_v1") {
        throw std::invalid_argument("unsupported calibration mapping");
    }
    require_finite(legacy_mapping.density_count_scale_2d_to_3d,
                   "calibration.density_count_scale_2d_to_3d");
    if (legacy_mapping.density_count_scale_2d_to_3d <= 0.0) {
        throw std::invalid_argument("calibration density scale must be positive");
    }
    const double source_values[] = {
        legacy_mapping.source_r_limit, legacy_mapping.source_K_limit,
        legacy_mapping.source_carrying_capacity_r,
        legacy_mapping.source_carrying_capacity_K,
    };
    for (const double value : source_values) {
        if (!std::isfinite(value) || value <= 0.0) {
            throw std::invalid_argument("legacy 2D density source values must be finite and positive");
        }
    }
    if (legacy_mapping.source_inner_radius < 0 ||
        legacy_mapping.source_outer_radius <= legacy_mapping.source_inner_radius) {
        throw std::invalid_argument("legacy 2D source radii are invalid");
    }
    if (!approximately_equal(carrying_capacity_scale_2d_to_3d, 1.0) ||
        !approximately_equal(r_limit, legacy_mapping.source_r_limit *
                                      legacy_mapping.density_count_scale_2d_to_3d) ||
        !approximately_equal(K_limit, legacy_mapping.source_K_limit *
                                      legacy_mapping.density_count_scale_2d_to_3d) ||
        !approximately_equal(carrying_capacity_r,
                             legacy_mapping.source_carrying_capacity_r *
                                 legacy_mapping.density_count_scale_2d_to_3d) ||
        !approximately_equal(carrying_capacity_K,
                             legacy_mapping.source_carrying_capacity_K *
                                 legacy_mapping.density_count_scale_2d_to_3d)) {
        throw std::invalid_argument("effective 3D density values do not match the one-time 2D mapping");
    }

    if (run_mode != "new" && run_mode != "resume") {
        throw std::invalid_argument("run.mode must be new or resume");
    }
    if ((run_mode == "resume") != !resume_checkpoint.empty()) {
        throw std::invalid_argument("run.resume_checkpoint must be set exactly when run.mode is resume");
    }

    if (domain_policy != "expandable_sparse" && domain_policy != "bounded") {
        throw std::invalid_argument("space.domain_policy must be expandable_sparse or bounded");
    }
    if ((domain_policy == "bounded") != bounded_domain) {
        throw std::invalid_argument("derived bounded domain state contradicts space.domain_policy");
    }
    if (chunk_edge < 4 || chunk_edge > 128) {
        throw std::invalid_argument("space.chunk_edge must be in [4,128]");
    }
    if (bounded_domain && !(domain_min.x <= domain_max.x && domain_min.y <= domain_max.y &&
                            domain_min.z <= domain_max.z)) {
        throw std::invalid_argument("domain minimum must not exceed maximum");
    }
    if (thin_layer && bounded_domain && (domain_min.z > 0 || domain_max.z < 0)) {
        throw std::invalid_argument("thin-layer bounded domain must include z=0");
    }

    if (direction_set != "fixed_26_v1") throw std::invalid_argument("unsupported direction.set");
    require_probability(continue_probability, "direction.continue_probability");
    require_probability(direction_density_threshold, "direction.density_threshold");
    require_finite(turn_half_angle_degrees, "direction.turn_half_angle_deg");
    require_finite(direction_density_half_angle_degrees,
                   "direction.density_half_angle_deg");
    require_finite(distance_weight_exponent, "direction.distance_weight_exponent");
    if (turn_half_angle_degrees <= 0.0 || turn_half_angle_degrees > 180.0 ||
        direction_density_half_angle_degrees <= 0.0 ||
        direction_density_half_angle_degrees > 180.0) {
        throw std::invalid_argument("direction angles must be in (0,180]");
    }
    if (direction_density_radius <= 0 || direction_density_radius > 128) {
        throw std::invalid_argument("direction density radius is invalid");
    }
    if (distance_weight_exponent < 0.0) {
        throw std::invalid_argument("distance weight exponent must be non-negative");
    }
    if (migration_activation_window_edge <= 0 || migration_activation_window_edge > 4096 ||
        migration_activation_block_edge <= 0 ||
        migration_activation_block_edge > migration_activation_window_edge) {
        throw std::invalid_argument("migration activation window/block configuration is invalid");
    }
    require_probability(migration_activation_threshold, "migration.activation_threshold");
    validate_beta_rate(normal_r_migration_beta, "migration.normal_r_rate");
    require_finite(migration_activation_duration_alpha,
                   "migration.activation_duration.alpha");
    require_probability(migration_activation_duration_mean_fraction,
                        "migration.activation_duration.mean_fraction");
    require_finite(migration_activation_duration_beta,
                   "migration.activation_duration.beta");
    if (normal_r_migration_beta.lower_clamp_enabled ||
        migration_activation_duration_alpha <= 0.0 ||
        migration_activation_duration_mean_fraction <= 0.0 ||
        migration_activation_duration_mean_fraction >= 1.0 ||
        migration_activation_duration_beta <= 0.0 ||
        !approximately_equal(
            migration_activation_duration_beta,
            migration_activation_duration_alpha *
                (1.0 - migration_activation_duration_mean_fraction) /
                migration_activation_duration_mean_fraction)) {
        throw std::invalid_argument(
            "migration normal-rate or activation-duration beta configuration is contradictory");
    }

    if (density_backend != "block_anchor_v1") throw std::invalid_argument("unsupported density.backend");
    if (density_block_edge <= 0 || density_block_edge > chunk_edge ||
        growth_density_window_edge <= 0 || growth_density_window_edge > 256) {
        throw std::invalid_argument("density block/window configuration is invalid");
    }
    if (r_limit <= 0.0 || K_limit <= 0.0 || carrying_capacity_r <= 0.0 ||
        carrying_capacity_K <= 0.0 || alpha < 0.0 || beta < 0.0) {
        throw std::invalid_argument("biology density parameters are invalid");
    }
    require_finite(alpha, "biology.alpha");
    require_finite(beta, "biology.beta");
    if (r_to_K_conversion.density_window_edge <= 0 ||
        r_to_K_conversion.density_window_edge > 4096 ||
        r_to_K_conversion.query_block_edge <= 0 ||
        r_to_K_conversion.query_block_edge >
            r_to_K_conversion.density_window_edge) {
        throw std::invalid_argument(
            "biology r-to-K conversion density window/block configuration is invalid");
    }
    require_probability(r_to_K_conversion.density_threshold,
                        "biology.r_to_K_conversion.density_threshold");
    require_probability(r_to_K_conversion.probability_per_division,
                        "biology.r_to_K_conversion.probability_per_division");

    if (large_footprint_edge != 2 || small_footprint_voxels != 1) {
        throw std::invalid_argument("unsupported stage footprint configuration");
    }
    if (division_shell_radius != 2) {
        throw std::invalid_argument("division shell radius must be 2");
    }

    if (initialization_mode != "explicit_counts" &&
        initialization_mode != "legacy_geometry_fill_v1") {
        throw std::invalid_argument("unsupported initial.mode");
    }
    if (initial_radius < 0 || initial_shell_thickness < 0 ||
        initial_shell_thickness > initial_radius) {
        throw std::invalid_argument("initial geometry radii are invalid");
    }
    if (initialization_mode == "legacy_geometry_fill_v1" &&
        (initial_shell_inner_radius < 0 || initial_inner_small_radius < 0 ||
         initial_shell_inner_radius >= initial_radius ||
         initial_inner_small_radius > initial_radius ||
         initial_shell_thickness != initial_radius - initial_shell_inner_radius)) {
        throw std::invalid_argument("geometry-fill initialization radii are inconsistent");
    }
    require_probability(initial_r_fraction, "initial.r_fraction");
    require_probability(initial_large_fraction, "initial.large_fraction");
    if (initialization_mode == "legacy_geometry_fill_v1" &&
        (initial_r_cells != 0 || initial_K_cells != 0)) {
        throw std::invalid_argument("geometry-fill initialization must derive counts, not specify them");
    }
    if (initial_r_cells > std::numeric_limits<std::uint64_t>::max() - initial_K_cells ||
        initial_r_cells + initial_K_cells >= static_cast<std::uint64_t>(kEmptySlot)) {
        throw std::invalid_argument("initial cell count exceeds stable-slot capacity");
    }
    if (initial_growth_rate_model == "fixed") {
        if (initial_r_growth_rate <= 0.0 || initial_K_growth_rate <= 0.0) {
            throw std::invalid_argument("fixed initial growth rates must be positive");
        }
    } else if (initial_growth_rate_model == "legacy_truncated_normal_v1") {
        validate_truncated_normal(initial_r_growth_truncated_normal,
                                  "initial.growth_rate.truncated_normal.r");
        validate_truncated_normal(initial_K_growth_truncated_normal,
                                  "initial.growth_rate.truncated_normal.K");
    } else {
        throw std::invalid_argument("unsupported initial.growth_rate.model");
    }
    if (initial_migration_rate_model == "fixed") {
        if (initial_r_migration_rate < 0.0 || initial_K_migration_rate < 0.0) {
            throw std::invalid_argument("fixed initial migration rates must be non-negative");
        }
    } else if (initial_migration_rate_model == "legacy_beta_v1") {
        validate_beta_rate(initial_r_migration_beta, "initial.migration_rate.beta.r");
        validate_beta_rate(initial_K_migration_beta, "initial.migration_rate.beta.K");
    } else {
        throw std::invalid_argument("unsupported initial.migration_rate.model");
    }
    if (death_delay_model != "legacy_geometric_mean_v1" ||
        death_growth_rate_threshold < 0.0 ||
        r_death_delay_hours <= 0.0 || K_death_delay_hours <= 0.0) {
        throw std::invalid_argument("initial biology rates or death delays are invalid");
    }

    require_finite(division_timing.base_cycle_hours, "division.timing.base_cycle_hours");
    require_probability(division_timing.minimum_fraction, "division.timing.minimum_fraction");
    require_probability(division_timing.stochastic_tail_fraction,
                        "division.timing.stochastic_tail_fraction");
    if (division_timing.base_cycle_hours <= 0.0 ||
        !approximately_equal(division_timing.minimum_fraction +
                                 division_timing.stochastic_tail_fraction,
                             1.0) ||
        division_timing.stochastic_time_quantum_hours <= 0.0 ||
        division_timing.retry_delay_hours <= 0.0 ||
        division_timing.inherited_growth_multiplier_min <= 0.0 ||
        division_timing.inherited_growth_multiplier_max <
            division_timing.inherited_growth_multiplier_min ||
        division_timing.r_max_inherent_growth_rate <= 0.0 ||
        division_timing.K_max_inherent_growth_rate <= 0.0) {
        throw std::invalid_argument("division timing configuration is invalid");
    }

    if (end_time_hours < 0.0 || max_events == 0 || threads <= 0 ||
        scheduler_backend != "event_queue_v1" || conflict_bucket_hours < 0.0) {
        throw std::invalid_argument("simulation or scheduler configuration is invalid");
    }
    if (preview_mode != "stable_uid_hash_v1" || full_format != "vtkhdf_points_v1" ||
        checkpoint_format != "hdf5_v2") {
        throw std::invalid_argument("unsupported output strategy");
    }
    if (output_enabled && output_directory.empty()) {
        throw std::invalid_argument("output directory must not be empty");
    }
    if (preview_every_hours < 0.0 || full_every_hours < 0.0 ||
        checkpoint_every_hours < 0.0 || preview_max_cells == 0) {
        throw std::invalid_argument("output intervals or preview size are invalid");
    }
    if (!(display_radius.large > 0.0F && display_radius.small > 0.0F &&
          display_radius.ultrasmall > 0.0F)) {
        throw std::invalid_argument("visualization display radii must be positive");
    }

    if (angiogenesis.trigger_metric != "biological_cell_volume_voxels3") {
        throw std::invalid_argument("unsupported angiogenesis trigger.metric");
    }
    if (angiogenesis.trigger_activation_volume_voxels3 <= 0.0 ||
        angiogenesis.trigger_deactivation_volume_voxels3 < 0.0 ||
        angiogenesis.trigger_deactivation_volume_voxels3 >
            angiogenesis.trigger_activation_volume_voxels3 ||
        angiogenesis.trigger_delay_hours < 0.0 ||
        angiogenesis.stage0_biological_volume_voxels3 <= 0.0 ||
        angiogenesis.stage1_biological_volume_voxels3 <= 0.0 ||
        angiogenesis.stage2_biological_volume_voxels3 <= 0.0) {
        throw std::invalid_argument("angiogenesis trigger/biological volume configuration is invalid");
    }
    if (angiogenesis.seed_process_model != "homogeneous_poisson" ||
        angiogenesis.seed_rate_sites_per_30_days <= 0.0 ||
        !approximately_equal(angiogenesis.seed_rate_sites_per_hour,
                             angiogenesis.seed_rate_sites_per_30_days / 720.0) ||
        angiogenesis.roots_per_event != 1 ||
        angiogenesis.surface_min_separation_voxels < 0 ||
        angiogenesis.surface_max_sampling_attempts == 0 ||
        angiogenesis.max_total_roots == 0 || angiogenesis.max_active_tips < 2) {
        throw std::invalid_argument("angiogenesis homogeneous Poisson seed process is invalid");
    }
    if (angiogenesis.diameter_voxels <= 0.0 ||
        angiogenesis.inward_speed_voxels_per_hour <= 0.0 ||
        angiogenesis.outward_speed_voxels_per_hour <= 0.0 ||
        angiogenesis.inward_max_length_voxels <= 0 ||
        angiogenesis.outward_max_length_voxels <= 0 ||
        angiogenesis.inward_target_tolerance_voxels < 0.0 ||
        angiogenesis.outward_external_connection_distance_voxels <= 0.0) {
        throw std::invalid_argument("angiogenesis vessel geometry/growth configuration is invalid");
    }
    if (angiogenesis.direction_model != "forward_biased_26_v1" ||
        angiogenesis.direction_forward_bias < 0.0 ||
        angiogenesis.direction_half_angle_degrees <= 0.0 ||
        angiogenesis.direction_half_angle_degrees > 180.0 ||
        angiogenesis.direction_turn_half_angle_degrees <= 0.0 ||
        angiogenesis.direction_turn_half_angle_degrees > 180.0 ||
        angiogenesis.direction_distance_weight_exponent < 0.0) {
        throw std::invalid_argument("angiogenesis direction configuration is invalid");
    }
    require_probability(angiogenesis.direction_persistence_probability,
                        "angiogenesis.direction.persistence_probability");
    if (angiogenesis.inward_replacement_policy != "remove_whole_cells" ||
        angiogenesis.outward_occupancy_policy != "permanent_reserved_placeholder" ||
        angiogenesis.vessel_collision_policy != "anastomose_and_stop" ||
        angiogenesis.boundary_policy != "stop") {
        throw std::invalid_argument("unsupported angiogenesis occupancy/collision policy");
    }
    if (angiogenesis.influence_profile != "linear_cutoff" &&
        angiogenesis.influence_profile != "exponential") {
        throw std::invalid_argument("unsupported angiogenesis influence.profile");
    }
    require_probability(angiogenesis.influence_max_relief_fraction,
                        "angiogenesis.influence.max_relief_fraction");
    if (angiogenesis.influence_max_relief_fraction >= 1.0 ||
        angiogenesis.influence_decay_length_voxels <= 0.0 ||
        angiogenesis.influence_cutoff_radius_voxels < angiogenesis.diameter_voxels * 0.5 ||
        angiogenesis.influence_scope != "growth_density" ||
        (angiogenesis.influence_activation != "after_outward_connection" &&
         angiogenesis.influence_activation != "immediate")) {
        throw std::invalid_argument("angiogenesis influence configuration is invalid");
    }
    if (angiogenesis.influence_activation == "after_outward_connection" &&
        angiogenesis.outward_external_connection_distance_voxels >
            static_cast<double>(angiogenesis.outward_max_length_voxels)) {
        throw std::invalid_argument(
            "angiogenesis outward connection distance exceeds outward maximum length");
    }

    const double finite_values[] = {
        continue_probability, turn_half_angle_degrees,
        direction_density_half_angle_degrees, direction_density_threshold,
        distance_weight_exponent, r_limit, K_limit, carrying_capacity_r,
        carrying_capacity_K, alpha, beta, initial_r_growth_rate,
        initial_K_growth_rate, initial_r_migration_rate, initial_K_migration_rate,
        initial_r_growth_truncated_normal.mean,
        initial_r_growth_truncated_normal.standard_deviation,
        initial_r_growth_truncated_normal.minimum,
        initial_r_growth_truncated_normal.maximum,
        initial_K_growth_truncated_normal.mean,
        initial_K_growth_truncated_normal.standard_deviation,
        initial_K_growth_truncated_normal.minimum,
        initial_K_growth_truncated_normal.maximum,
        initial_r_migration_beta.alpha, initial_r_migration_beta.beta,
        initial_r_migration_beta.scale,
        initial_r_migration_beta.lower_clamp_threshold,
        initial_r_migration_beta.lower_clamp_value,
        initial_K_migration_beta.alpha, initial_K_migration_beta.beta,
        initial_K_migration_beta.scale,
        initial_K_migration_beta.lower_clamp_threshold,
        initial_K_migration_beta.lower_clamp_value,
        death_growth_rate_threshold, r_death_delay_hours, K_death_delay_hours,
        r_to_K_conversion.density_threshold,
        r_to_K_conversion.probability_per_division,
        initial_large_fraction,
        initial_r_fraction, division_timing.stochastic_time_quantum_hours,
        division_timing.retry_delay_hours,
        division_timing.inherited_growth_multiplier_min,
        division_timing.inherited_growth_multiplier_max,
        division_timing.r_max_inherent_growth_rate,
        division_timing.K_max_inherent_growth_rate, end_time_hours,
        conflict_bucket_hours, preview_every_hours, full_every_hours,
        checkpoint_every_hours, migration_activation_threshold,
        normal_r_migration_beta.alpha, normal_r_migration_beta.beta,
        normal_r_migration_beta.scale,
        migration_activation_duration_alpha,
        migration_activation_duration_mean_fraction,
        migration_activation_duration_beta,
        angiogenesis.trigger_activation_volume_voxels3,
        angiogenesis.trigger_deactivation_volume_voxels3,
        angiogenesis.trigger_delay_hours,
        angiogenesis.stage0_biological_volume_voxels3,
        angiogenesis.stage1_biological_volume_voxels3,
        angiogenesis.stage2_biological_volume_voxels3,
        angiogenesis.seed_rate_sites_per_30_days,
        angiogenesis.seed_rate_sites_per_hour, angiogenesis.diameter_voxels,
        angiogenesis.inward_speed_voxels_per_hour,
        angiogenesis.outward_speed_voxels_per_hour,
        angiogenesis.inward_target_tolerance_voxels,
        angiogenesis.outward_external_connection_distance_voxels,
        angiogenesis.direction_forward_bias,
        angiogenesis.direction_half_angle_degrees,
        angiogenesis.direction_turn_half_angle_degrees,
        angiogenesis.direction_persistence_probability,
        angiogenesis.direction_distance_weight_exponent,
        angiogenesis.influence_max_relief_fraction,
        angiogenesis.influence_decay_length_voxels,
        angiogenesis.influence_cutoff_radius_voxels,
    };
    for (const double value : finite_values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("configuration contains NaN or infinity");
        }
    }
}

std::string Model3DConfig::to_json() const {
    std::ostringstream out;
    out << std::setprecision(17);
    out << '{'
        << "\"schema_name\":\"" << json_escape(schema_name) << "\","
        << "\"schema_version\":" << schema_version << ','
        << "\"profile\":\"" << json_escape(profile) << "\","
        << "\"calibration_source_profile\":\"" << json_escape(legacy_mapping.source_profile) << "\","
        << "\"calibration_conversion\":\"" << json_escape(legacy_mapping.conversion) << "\","
        << "\"density_count_scale_2d_to_3d\":" << legacy_mapping.density_count_scale_2d_to_3d << ','
        << "\"source_2d_r_limit\":" << legacy_mapping.source_r_limit << ','
        << "\"source_2d_K_limit\":" << legacy_mapping.source_K_limit << ','
        << "\"source_2d_carrying_capacity_r\":" << legacy_mapping.source_carrying_capacity_r << ','
        << "\"source_2d_carrying_capacity_K\":" << legacy_mapping.source_carrying_capacity_K << ','
        << "\"source_2d_outer_radius\":" << legacy_mapping.source_outer_radius << ','
        << "\"source_2d_inner_radius\":" << legacy_mapping.source_inner_radius << ','
        << "\"scaled_down_initialization\":" << json_bool(legacy_mapping.scaled_down_initialization) << ','
        << "\"seed\":" << seed << ','
        << "\"run_mode\":\"" << json_escape(run_mode) << "\","
        << "\"resume_checkpoint\":";
    if (resume_checkpoint.empty()) out << "null";
    else out << '"' << json_escape(resume_checkpoint.string()) << '"';
    out << ','
        << "\"domain_policy\":\"" << json_escape(domain_policy) << "\","
        << "\"bounded_domain\":" << json_bool(bounded_domain) << ','
        << "\"domain_min\":[" << domain_min.x << ',' << domain_min.y << ',' << domain_min.z << "],"
        << "\"domain_max\":[" << domain_max.x << ',' << domain_max.y << ',' << domain_max.z << "],"
        << "\"chunk_edge\":" << chunk_edge << ','
        << "\"thin_layer\":" << json_bool(thin_layer) << ','
        << "\"direction_set\":\"" << json_escape(direction_set) << "\","
        << "\"continue_probability\":" << continue_probability << ','
        << "\"turn_half_angle_degrees\":" << turn_half_angle_degrees << ','
        << "\"direction_density_radius\":" << direction_density_radius << ','
        << "\"direction_density_half_angle_degrees\":" << direction_density_half_angle_degrees << ','
        << "\"direction_density_threshold\":" << direction_density_threshold << ','
        << "\"persistence_uses_density\":" << json_bool(persistence_uses_density) << ','
        << "\"distance_weight_exponent\":" << distance_weight_exponent << ','
        << "\"migration_activation_enabled\":" << json_bool(migration_activation_enabled) << ','
        << "\"migration_activation_window_edge\":" << migration_activation_window_edge << ','
        << "\"migration_activation_block_edge\":" << migration_activation_block_edge << ','
        << "\"migration_activation_threshold\":" << migration_activation_threshold << ','
        << "\"normal_r_migration_beta_alpha\":" << normal_r_migration_beta.alpha << ','
        << "\"normal_r_migration_beta_beta\":" << normal_r_migration_beta.beta << ','
        << "\"normal_r_migration_beta_scale\":" << normal_r_migration_beta.scale << ','
        << "\"migration_activation_duration_alpha\":" << migration_activation_duration_alpha << ','
        << "\"migration_activation_duration_mean_fraction\":" << migration_activation_duration_mean_fraction << ','
        << "\"migration_activation_duration_beta\":" << migration_activation_duration_beta << ','
        << "\"density_backend\":\"" << json_escape(density_backend) << "\","
        << "\"density_block_edge\":" << density_block_edge << ','
        << "\"growth_density_window_edge\":" << growth_density_window_edge << ','
        << "\"carrying_capacity_scale_2d_to_3d\":" << carrying_capacity_scale_2d_to_3d << ','
        << "\"r_limit\":" << r_limit << ','
        << "\"K_limit\":" << K_limit << ','
        << "\"carrying_capacity_r\":" << carrying_capacity_r << ','
        << "\"carrying_capacity_K\":" << carrying_capacity_K << ','
        << "\"alpha\":" << alpha << ','
        << "\"beta\":" << beta << ','
        << "\"large_footprint_edge\":" << large_footprint_edge << ','
        << "\"small_footprint_voxels\":" << small_footprint_voxels << ','
        << "\"division_shell_radius\":" << division_shell_radius << ','
        << "\"allow_shape_reduction\":" << json_bool(allow_shape_reduction) << ','
        << "\"ultrasmall_enabled\":" << json_bool(ultrasmall_enabled) << ','
        << "\"initialization_mode\":\"" << json_escape(initialization_mode) << "\","
        << "\"initial_r_cells\":" << initial_r_cells << ','
        << "\"initial_K_cells\":" << initial_K_cells << ','
        << "\"initial_radius\":" << initial_radius << ','
        << "\"initial_shell_thickness\":" << initial_shell_thickness << ','
        << "\"initial_shell_inner_radius\":" << initial_shell_inner_radius << ','
        << "\"initial_inner_small_radius\":" << initial_inner_small_radius << ','
        << "\"initial_r_fraction\":" << initial_r_fraction << ','
        << "\"initial_large_fraction\":" << initial_large_fraction << ','
        << "\"initial_growth_rate_model\":\"" << json_escape(initial_growth_rate_model) << "\","
        << "\"initial_r_growth_rate\":" << initial_r_growth_rate << ','
        << "\"initial_K_growth_rate\":" << initial_K_growth_rate << ','
        << "\"initial_r_growth_truncated_normal_mean\":" << initial_r_growth_truncated_normal.mean << ','
        << "\"initial_r_growth_truncated_normal_standard_deviation\":" << initial_r_growth_truncated_normal.standard_deviation << ','
        << "\"initial_r_growth_truncated_normal_minimum\":" << initial_r_growth_truncated_normal.minimum << ','
        << "\"initial_r_growth_truncated_normal_maximum\":" << initial_r_growth_truncated_normal.maximum << ','
        << "\"initial_K_growth_truncated_normal_mean\":" << initial_K_growth_truncated_normal.mean << ','
        << "\"initial_K_growth_truncated_normal_standard_deviation\":" << initial_K_growth_truncated_normal.standard_deviation << ','
        << "\"initial_K_growth_truncated_normal_minimum\":" << initial_K_growth_truncated_normal.minimum << ','
        << "\"initial_K_growth_truncated_normal_maximum\":" << initial_K_growth_truncated_normal.maximum << ','
        << "\"initial_migration_rate_model\":\"" << json_escape(initial_migration_rate_model) << "\","
        << "\"initial_r_migration_rate\":" << initial_r_migration_rate << ','
        << "\"initial_K_migration_rate\":" << initial_K_migration_rate << ','
        << "\"initial_r_migration_beta_alpha\":" << initial_r_migration_beta.alpha << ','
        << "\"initial_r_migration_beta_beta\":" << initial_r_migration_beta.beta << ','
        << "\"initial_r_migration_beta_scale\":" << initial_r_migration_beta.scale << ','
        << "\"initial_r_migration_beta_lower_clamp_enabled\":" << json_bool(initial_r_migration_beta.lower_clamp_enabled) << ','
        << "\"initial_r_migration_beta_lower_clamp_threshold\":" << initial_r_migration_beta.lower_clamp_threshold << ','
        << "\"initial_r_migration_beta_lower_clamp_value\":" << initial_r_migration_beta.lower_clamp_value << ','
        << "\"initial_K_migration_beta_alpha\":" << initial_K_migration_beta.alpha << ','
        << "\"initial_K_migration_beta_beta\":" << initial_K_migration_beta.beta << ','
        << "\"initial_K_migration_beta_scale\":" << initial_K_migration_beta.scale << ','
        << "\"initial_K_migration_beta_lower_clamp_enabled\":" << json_bool(initial_K_migration_beta.lower_clamp_enabled) << ','
        << "\"initial_K_migration_beta_lower_clamp_threshold\":" << initial_K_migration_beta.lower_clamp_threshold << ','
        << "\"initial_K_migration_beta_lower_clamp_value\":" << initial_K_migration_beta.lower_clamp_value << ','
        << "\"death_delay_model\":\"" << json_escape(death_delay_model) << "\","
        << "\"death_growth_rate_threshold\":" << death_growth_rate_threshold << ','
        << "\"r_death_delay_hours\":" << r_death_delay_hours << ','
        << "\"K_death_delay_hours\":" << K_death_delay_hours << ','
        << "\"r_to_K_conversion_enabled\":" << json_bool(r_to_K_conversion.enabled) << ','
        << "\"r_to_K_conversion_density_window_edge\":" << r_to_K_conversion.density_window_edge << ','
        << "\"r_to_K_conversion_query_block_edge\":" << r_to_K_conversion.query_block_edge << ','
        << "\"r_to_K_conversion_density_threshold\":" << r_to_K_conversion.density_threshold << ','
        << "\"r_to_K_conversion_probability_per_division\":" << r_to_K_conversion.probability_per_division << ','
        << "\"division_base_cycle_hours\":" << division_timing.base_cycle_hours << ','
        << "\"division_minimum_fraction\":" << division_timing.minimum_fraction << ','
        << "\"division_stochastic_tail_fraction\":" << division_timing.stochastic_tail_fraction << ','
        << "\"division_stochastic_time_quantum_hours\":" << division_timing.stochastic_time_quantum_hours << ','
        << "\"division_retry_delay_hours\":" << division_timing.retry_delay_hours << ','
        << "\"division_inherited_growth_multiplier_min\":" << division_timing.inherited_growth_multiplier_min << ','
        << "\"division_inherited_growth_multiplier_max\":" << division_timing.inherited_growth_multiplier_max << ','
        << "\"division_r_max_inherent_growth_rate\":" << division_timing.r_max_inherent_growth_rate << ','
        << "\"division_K_max_inherent_growth_rate\":" << division_timing.K_max_inherent_growth_rate << ','
        << "\"end_time_hours\":" << end_time_hours << ','
        << "\"max_events\":" << max_events << ','
        << "\"threads\":" << threads << ','
        << "\"scheduler_backend\":\"" << json_escape(scheduler_backend) << "\","
        << "\"conflict_bucket_hours\":" << conflict_bucket_hours << ','
        << "\"output_enabled\":" << json_bool(output_enabled) << ','
        << "\"preview_mode\":\"" << json_escape(preview_mode) << "\","
        << "\"full_format\":\"" << json_escape(full_format) << "\","
        << "\"checkpoint_format\":\"" << json_escape(checkpoint_format) << "\","
        << "\"output_directory\":\"" << json_escape(output_directory.string()) << "\","
        << "\"preview_every_hours\":" << preview_every_hours << ','
        << "\"full_every_hours\":" << full_every_hours << ','
        << "\"checkpoint_every_hours\":" << checkpoint_every_hours << ','
        << "\"preview_max_cells\":" << preview_max_cells << ','
        << "\"preview_seed\":" << preview_seed << ','
        << "\"display_radius_large\":" << display_radius.large << ','
        << "\"display_radius_small\":" << display_radius.small << ','
        << "\"display_radius_ultrasmall\":" << display_radius.ultrasmall << ','
        << "\"angiogenesis_enabled\":" << json_bool(angiogenesis.enabled) << ','
        << "\"angiogenesis_trigger_metric\":\"" << json_escape(angiogenesis.trigger_metric) << "\","
        << "\"angiogenesis_trigger_activation_volume_voxels3\":" << angiogenesis.trigger_activation_volume_voxels3 << ','
        << "\"angiogenesis_trigger_deactivation_volume_voxels3\":" << angiogenesis.trigger_deactivation_volume_voxels3 << ','
        << "\"angiogenesis_trigger_delay_hours\":" << angiogenesis.trigger_delay_hours << ','
        << "\"angiogenesis_stage0_biological_volume_voxels3\":" << angiogenesis.stage0_biological_volume_voxels3 << ','
        << "\"angiogenesis_stage1_biological_volume_voxels3\":" << angiogenesis.stage1_biological_volume_voxels3 << ','
        << "\"angiogenesis_stage2_biological_volume_voxels3\":" << angiogenesis.stage2_biological_volume_voxels3 << ','
        << "\"angiogenesis_seed_process_model\":\"" << json_escape(angiogenesis.seed_process_model) << "\","
        << "\"angiogenesis_seed_rate_sites_per_30_days\":" << angiogenesis.seed_rate_sites_per_30_days << ','
        << "\"angiogenesis_seed_rate_sites_per_hour\":" << angiogenesis.seed_rate_sites_per_hour << ','
        << "\"angiogenesis_roots_per_event\":" << angiogenesis.roots_per_event << ','
        << "\"angiogenesis_surface_min_separation_voxels\":" << angiogenesis.surface_min_separation_voxels << ','
        << "\"angiogenesis_surface_max_sampling_attempts\":" << angiogenesis.surface_max_sampling_attempts << ','
        << "\"angiogenesis_max_total_roots\":" << angiogenesis.max_total_roots << ','
        << "\"angiogenesis_max_active_tips\":" << angiogenesis.max_active_tips << ','
        << "\"angiogenesis_diameter_voxels\":" << angiogenesis.diameter_voxels << ','
        << "\"angiogenesis_inward_speed_voxels_per_hour\":" << angiogenesis.inward_speed_voxels_per_hour << ','
        << "\"angiogenesis_outward_speed_voxels_per_hour\":" << angiogenesis.outward_speed_voxels_per_hour << ','
        << "\"angiogenesis_inward_max_length_voxels\":" << angiogenesis.inward_max_length_voxels << ','
        << "\"angiogenesis_outward_max_length_voxels\":" << angiogenesis.outward_max_length_voxels << ','
        << "\"angiogenesis_inward_target_tolerance_voxels\":" << angiogenesis.inward_target_tolerance_voxels << ','
        << "\"angiogenesis_outward_external_connection_distance_voxels\":" << angiogenesis.outward_external_connection_distance_voxels << ','
        << "\"angiogenesis_direction_model\":\"" << json_escape(angiogenesis.direction_model) << "\","
        << "\"angiogenesis_direction_forward_bias\":" << angiogenesis.direction_forward_bias << ','
        << "\"angiogenesis_direction_half_angle_degrees\":" << angiogenesis.direction_half_angle_degrees << ','
        << "\"angiogenesis_direction_turn_half_angle_degrees\":" << angiogenesis.direction_turn_half_angle_degrees << ','
        << "\"angiogenesis_direction_persistence_probability\":" << angiogenesis.direction_persistence_probability << ','
        << "\"angiogenesis_direction_distance_weight_exponent\":" << angiogenesis.direction_distance_weight_exponent << ','
        << "\"angiogenesis_inward_replacement_policy\":\"" << json_escape(angiogenesis.inward_replacement_policy) << "\","
        << "\"angiogenesis_outward_occupancy_policy\":\"" << json_escape(angiogenesis.outward_occupancy_policy) << "\","
        << "\"angiogenesis_vessel_collision_policy\":\"" << json_escape(angiogenesis.vessel_collision_policy) << "\","
        << "\"angiogenesis_boundary_policy\":\"" << json_escape(angiogenesis.boundary_policy) << "\","
        << "\"angiogenesis_influence_profile\":\"" << json_escape(angiogenesis.influence_profile) << "\","
        << "\"angiogenesis_influence_max_relief_fraction\":" << angiogenesis.influence_max_relief_fraction << ','
        << "\"angiogenesis_influence_decay_length_voxels\":" << angiogenesis.influence_decay_length_voxels << ','
        << "\"angiogenesis_influence_cutoff_radius_voxels\":" << angiogenesis.influence_cutoff_radius_voxels << ','
        << "\"angiogenesis_influence_scope\":\"" << json_escape(angiogenesis.influence_scope) << "\","
        << "\"angiogenesis_influence_activation\":\"" << json_escape(angiogenesis.influence_activation) << "\""
        << "}\n";
    return out.str();
}

std::string Model3DConfig::dynamics_json() const {
    Model3DConfig normalized = *this;
    normalized.run_mode = "new";
    normalized.resume_checkpoint.clear();
    normalized.end_time_hours = 0.0;
    normalized.max_events = 1;
    normalized.threads = 1;
    normalized.output_enabled = false;
    normalized.preview_mode = "stable_uid_hash_v1";
    normalized.full_format = "vtkhdf_points_v1";
    normalized.checkpoint_format = "hdf5_v2";
    normalized.output_directory.clear();
    normalized.preview_every_hours = 0.0;
    normalized.full_every_hours = 0.0;
    normalized.checkpoint_every_hours = 0.0;
    normalized.preview_max_cells = 1;
    normalized.preview_seed = 0;
    normalized.display_radius = {};
    return normalized.to_json();
}

}  // namespace atcg3d
