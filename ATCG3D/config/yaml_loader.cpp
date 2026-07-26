#include "config/model_config.hpp"

#include <yaml-cpp/yaml.h>

#include <charconv>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace atcg3d {
namespace {

[[noreturn]] void config_error(const YAML::Node& node,
                               const std::string& path,
                               const std::string& message) {
    const YAML::Mark mark = node.Mark();
    std::string location;
    if (!mark.is_null()) {
        location = " (line " + std::to_string(mark.line + 1) + ", column " +
                   std::to_string(mark.column + 1) + ')';
    }
    throw std::invalid_argument(path + location + ": " + message);
}

void check_map(const YAML::Node& node,
               const std::string& path,
               std::initializer_list<std::string_view> allowed_keys) {
    if (!node || !node.IsMap()) {
        config_error(node, path, "must be a mapping");
    }
    std::unordered_set<std::string> allowed;
    for (const std::string_view key : allowed_keys) allowed.emplace(key);
    std::unordered_set<std::string> seen;
    for (const auto& entry : node) {
        if (!entry.first.IsScalar()) {
            config_error(entry.first, path, "mapping keys must be scalar strings");
        }
        const std::string key = entry.first.Scalar();
        if (!seen.insert(key).second) {
            config_error(entry.first, path + '.' + key, "duplicate key");
        }
        if (!allowed.contains(key)) {
            config_error(entry.first, path + '.' + key, "unknown key");
        }
    }
}

YAML::Node required(const YAML::Node& map,
                    const std::string& key,
                    const std::string& parent_path) {
    const YAML::Node value = map[key];
    if (!value || !value.IsDefined()) {
        config_error(map, parent_path + '.' + key, "required key is missing");
    }
    return value;
}

std::string scalar_text(const YAML::Node& node, const std::string& path) {
    if (!node || !node.IsScalar()) config_error(node, path, "must be a scalar");
    return node.Scalar();
}

std::string plain_scalar_text(const YAML::Node& node, const std::string& path) {
    const std::string value = scalar_text(node, path);
    const std::string tag = node.Tag();
    if (tag == "!" || (!tag.empty() && tag != "?")) {
        config_error(node, path, "must use an unquoted, untagged YAML scalar");
    }
    return value;
}

std::string strict_string(const YAML::Node& node, const std::string& path) {
    return scalar_text(node, path);
}

bool strict_bool(const YAML::Node& node, const std::string& path) {
    const std::string value = plain_scalar_text(node, path);
    if (value == "true") return true;
    if (value == "false") return false;
    config_error(node, path, "must be the YAML boolean true or false");
}

template <class Integer>
Integer strict_integer(const YAML::Node& node, const std::string& path) {
    const std::string value = plain_scalar_text(node, path);
    if (value.empty() || value.front() == '+' || value.find_first_of("_xXoO") != std::string::npos) {
        config_error(node, path, "must be a decimal integer");
    }
    Integer result{};
    const char* begin = value.data();
    const char* end = begin + value.size();
    const auto [pointer, error] = std::from_chars(begin, end, result, 10);
    if (error != std::errc{} || pointer != end) {
        config_error(node, path, "integer is invalid or outside the supported range");
    }
    return result;
}

double strict_double(const YAML::Node& node, const std::string& path) {
    const std::string value = plain_scalar_text(node, path);
    if (value.empty() || value.find('_') != std::string::npos) {
        config_error(node, path, "must be a finite floating point value");
    }
    std::size_t consumed = 0;
    double result = 0.0;
    try {
        result = std::stod(value, &consumed);
    } catch (const std::exception&) {
        config_error(node, path, "must be a finite floating point value");
    }
    if (consumed != value.size() || !std::isfinite(result)) {
        config_error(node, path, "must be a finite floating point value");
    }
    return result;
}

float strict_float(const YAML::Node& node, const std::string& path) {
    const double value = strict_double(node, path);
    if (value < -static_cast<double>(std::numeric_limits<float>::max()) ||
        value > static_cast<double>(std::numeric_limits<float>::max())) {
        config_error(node, path, "is outside the float32 range");
    }
    return static_cast<float>(value);
}

Vec3i strict_vec3i(const YAML::Node& node, const std::string& path) {
    if (!node || !node.IsSequence() || node.size() != 3) {
        config_error(node, path, "must be a sequence of exactly three int32 values");
    }
    return {
        strict_integer<std::int32_t>(node[0], path + "[0]"),
        strict_integer<std::int32_t>(node[1], path + "[1]"),
        strict_integer<std::int32_t>(node[2], path + "[2]"),
    };
}

std::filesystem::path optional_path(const YAML::Node& node,
                                    const std::string& path,
                                    const std::filesystem::path& base_directory) {
    if (!node || node.IsNull()) return {};
    std::filesystem::path value(strict_string(node, path));
    if (value.empty()) return {};
    if (value.is_relative()) value = base_directory / value;
    return value.lexically_normal();
}

YAML::Node checked_section(const YAML::Node& parent,
                           const std::string& key,
                           const std::string& parent_path,
                           std::initializer_list<std::string_view> allowed_keys) {
    YAML::Node section = required(parent, key, parent_path);
    check_map(section, parent_path + '.' + key, allowed_keys);
    return section;
}

TruncatedNormalRateConfig parse_truncated_normal_rate(const YAML::Node& parent,
                                                       const std::string& key,
                                                       const std::string& parent_path) {
    const YAML::Node node = checked_section(
        parent, key, parent_path,
        {"mean", "standard_deviation", "minimum", "maximum"});
    const std::string path = parent_path + '.' + key;
    TruncatedNormalRateConfig result;
    result.mean = strict_double(required(node, "mean", path), path + ".mean");
    result.standard_deviation = strict_double(
        required(node, "standard_deviation", path), path + ".standard_deviation");
    result.minimum = strict_double(required(node, "minimum", path), path + ".minimum");
    result.maximum = strict_double(required(node, "maximum", path), path + ".maximum");
    return result;
}

BetaRateConfig parse_beta_rate_fields(const YAML::Node& node,
                                      const std::string& path) {
    BetaRateConfig result;
    result.alpha = strict_double(required(node, "alpha", path), path + ".alpha");
    result.beta = strict_double(required(node, "beta", path), path + ".beta");
    result.scale = strict_double(required(node, "scale", path), path + ".scale");

    const YAML::Node clamp = required(node, "lower_clamp", path);
    if (!clamp || !clamp.IsMap()) {
        config_error(clamp, path + ".lower_clamp", "must be a mapping");
    }
    result.lower_clamp_enabled = strict_bool(
        required(clamp, "enabled", path + ".lower_clamp"),
        path + ".lower_clamp.enabled");
    if (result.lower_clamp_enabled) {
        check_map(clamp, path + ".lower_clamp", {"enabled", "threshold", "value"});
        result.lower_clamp_threshold = strict_double(
            required(clamp, "threshold", path + ".lower_clamp"),
            path + ".lower_clamp.threshold");
        result.lower_clamp_value = strict_double(
            required(clamp, "value", path + ".lower_clamp"),
            path + ".lower_clamp.value");
    } else {
        check_map(clamp, path + ".lower_clamp", {"enabled"});
    }
    return result;
}

BetaRateConfig parse_beta_rate(const YAML::Node& parent,
                               const std::string& key,
                               const std::string& parent_path) {
    const YAML::Node node = checked_section(
        parent, key, parent_path, {"alpha", "beta", "scale", "lower_clamp"});
    return parse_beta_rate_fields(node, parent_path + '.' + key);
}

}  // namespace

Model3DConfig Model3DConfig::load(const std::filesystem::path& path) {
    std::vector<YAML::Node> documents;
    try {
        documents = YAML::LoadAllFromFile(path.string());
    } catch (const YAML::Exception& error) {
        const YAML::Mark mark = error.mark;
        std::string location;
        if (!mark.is_null()) {
            location = ":" + std::to_string(mark.line + 1) + ':' +
                       std::to_string(mark.column + 1);
        }
        throw std::invalid_argument("unable to parse YAML configuration " + path.string() +
                                    location + ": " + error.msg);
    }
    if (documents.size() != 1 || !documents.front()) {
        throw std::invalid_argument("configuration must contain exactly one YAML document: " +
                                    path.string());
    }

    const YAML::Node root = documents.front();
    check_map(root, "$", {"schema", "profile", "calibration", "rng", "run", "space",
                           "direction", "migration", "density", "biology", "stage",
                           "division", "initial", "simulation", "parallel", "scheduler",
                           "angiogenesis", "output", "control"});

    Model3DConfig config;
    const YAML::Node schema = checked_section(root, "schema", "$", {"name", "version"});
    config.schema_name = strict_string(required(schema, "name", "$.schema"), "schema.name");
    config.schema_version = strict_integer<std::uint32_t>(
        required(schema, "version", "$.schema"), "schema.version");
    config.profile = strict_string(required(root, "profile", "$"), "profile");

    const YAML::Node calibration = checked_section(
        root, "calibration", "$", {"source_profile", "conversion",
                                     "density_count_scale_2d_to_3d",
                                     "scaled_down_initialization", "source_2d"});
    config.legacy_mapping.source_profile = strict_string(
        required(calibration, "source_profile", "$.calibration"),
        "calibration.source_profile");
    config.legacy_mapping.conversion = strict_string(
        required(calibration, "conversion", "$.calibration"), "calibration.conversion");
    config.legacy_mapping.density_count_scale_2d_to_3d = strict_double(
        required(calibration, "density_count_scale_2d_to_3d", "$.calibration"),
        "calibration.density_count_scale_2d_to_3d");
    config.legacy_mapping.scaled_down_initialization = strict_bool(
        required(calibration, "scaled_down_initialization", "$.calibration"),
        "calibration.scaled_down_initialization");
    const YAML::Node source_2d = checked_section(
        calibration, "source_2d", "$.calibration",
        {"r_limit", "K_limit", "carrying_capacity_r", "carrying_capacity_K",
         "outer_radius", "inner_radius"});
    config.legacy_mapping.source_r_limit = strict_double(
        required(source_2d, "r_limit", "$.calibration.source_2d"),
        "calibration.source_2d.r_limit");
    config.legacy_mapping.source_K_limit = strict_double(
        required(source_2d, "K_limit", "$.calibration.source_2d"),
        "calibration.source_2d.K_limit");
    config.legacy_mapping.source_carrying_capacity_r = strict_double(
        required(source_2d, "carrying_capacity_r", "$.calibration.source_2d"),
        "calibration.source_2d.carrying_capacity_r");
    config.legacy_mapping.source_carrying_capacity_K = strict_double(
        required(source_2d, "carrying_capacity_K", "$.calibration.source_2d"),
        "calibration.source_2d.carrying_capacity_K");
    config.legacy_mapping.source_outer_radius = strict_integer<int>(
        required(source_2d, "outer_radius", "$.calibration.source_2d"),
        "calibration.source_2d.outer_radius");
    config.legacy_mapping.source_inner_radius = strict_integer<int>(
        required(source_2d, "inner_radius", "$.calibration.source_2d"),
        "calibration.source_2d.inner_radius");

    // Convert once at startup. Existing biology code retains a multiplier for
    // ABI compatibility, but it is fixed to one in every YAML-loaded config.
    const double scale = config.legacy_mapping.density_count_scale_2d_to_3d;
    config.carrying_capacity_scale_2d_to_3d = 1.0;
    config.r_limit = config.legacy_mapping.source_r_limit * scale;
    config.K_limit = config.legacy_mapping.source_K_limit * scale;
    config.carrying_capacity_r = config.legacy_mapping.source_carrying_capacity_r * scale;
    config.carrying_capacity_K = config.legacy_mapping.source_carrying_capacity_K * scale;

    const YAML::Node rng = checked_section(root, "rng", "$", {"seed"});
    config.seed = strict_integer<std::uint64_t>(required(rng, "seed", "$.rng"), "rng.seed");

    const YAML::Node run = checked_section(root, "run", "$", {"mode", "resume_checkpoint"});
    config.run_mode = strict_string(required(run, "mode", "$.run"), "run.mode");
    config.resume_checkpoint = optional_path(
        required(run, "resume_checkpoint", "$.run"), "run.resume_checkpoint",
        std::filesystem::absolute(path).parent_path());

    const YAML::Node space = checked_section(
        root, "space", "$", {"domain_policy", "chunk_edge", "thin_layer", "minimum", "maximum"});
    config.domain_policy = strict_string(required(space, "domain_policy", "$.space"),
                                         "space.domain_policy");
    config.bounded_domain = config.domain_policy == "bounded";
    config.chunk_edge = strict_integer<int>(required(space, "chunk_edge", "$.space"),
                                            "space.chunk_edge");
    config.thin_layer = strict_bool(required(space, "thin_layer", "$.space"),
                                    "space.thin_layer");
    config.domain_min = strict_vec3i(required(space, "minimum", "$.space"), "space.minimum");
    config.domain_max = strict_vec3i(required(space, "maximum", "$.space"), "space.maximum");

    const YAML::Node direction = checked_section(
        root, "direction", "$", {"set", "continue_probability", "turn_half_angle_deg",
                                   "density_radius", "density_half_angle_deg",
                                   "density_threshold", "persistence_uses_density",
                                   "distance_weight_exponent"});
    config.direction_set = strict_string(required(direction, "set", "$.direction"), "direction.set");
    config.continue_probability = strict_double(
        required(direction, "continue_probability", "$.direction"),
        "direction.continue_probability");
    config.turn_half_angle_degrees = strict_double(
        required(direction, "turn_half_angle_deg", "$.direction"),
        "direction.turn_half_angle_deg");
    config.direction_density_radius = strict_integer<int>(
        required(direction, "density_radius", "$.direction"), "direction.density_radius");
    config.direction_density_half_angle_degrees = strict_double(
        required(direction, "density_half_angle_deg", "$.direction"),
        "direction.density_half_angle_deg");
    config.direction_density_threshold = strict_double(
        required(direction, "density_threshold", "$.direction"),
        "direction.density_threshold");
    config.persistence_uses_density = strict_bool(
        required(direction, "persistence_uses_density", "$.direction"),
        "direction.persistence_uses_density");
    config.distance_weight_exponent = strict_double(
        required(direction, "distance_weight_exponent", "$.direction"),
        "direction.distance_weight_exponent");

    const YAML::Node migration = checked_section(
        root, "migration", "$", {"activation_enabled", "activation_window_edge",
                                   "activation_block_edge", "activation_threshold",
                                   "normal_r_rate", "activated_r_rate",
                                   "activation_duration",
                                   "crowding_exchange"});
    config.migration_activation_enabled = strict_bool(
        required(migration, "activation_enabled", "$.migration"),
        "migration.activation_enabled");
    config.migration_activation_window_edge = strict_integer<int>(
        required(migration, "activation_window_edge", "$.migration"),
        "migration.activation_window_edge");
    config.migration_activation_block_edge = strict_integer<int>(
        required(migration, "activation_block_edge", "$.migration"),
        "migration.activation_block_edge");
    config.migration_activation_threshold = strict_double(
        required(migration, "activation_threshold", "$.migration"),
        "migration.activation_threshold");
    const YAML::Node normal_r_rate = checked_section(
        migration, "normal_r_rate", "$.migration",
        {"model", "alpha", "beta", "scale"});
    if (strict_string(required(normal_r_rate, "model", "$.migration.normal_r_rate"),
                      "migration.normal_r_rate.model") != "beta") {
        config_error(required(normal_r_rate, "model", "$.migration.normal_r_rate"),
                     "migration.normal_r_rate.model", "unsupported model");
    }
    config.normal_r_migration_beta.alpha = strict_double(
        required(normal_r_rate, "alpha", "$.migration.normal_r_rate"),
        "migration.normal_r_rate.alpha");
    config.normal_r_migration_beta.beta = strict_double(
        required(normal_r_rate, "beta", "$.migration.normal_r_rate"),
        "migration.normal_r_rate.beta");
    config.normal_r_migration_beta.scale = strict_double(
        required(normal_r_rate, "scale", "$.migration.normal_r_rate"),
        "migration.normal_r_rate.scale");
    config.normal_r_migration_beta.lower_clamp_enabled = false;
    config.normal_r_migration_beta.lower_clamp_threshold = 0.0;
    config.normal_r_migration_beta.lower_clamp_value = 0.0;

    const YAML::Node activated_r_rate = checked_section(
        migration, "activated_r_rate", "$.migration",
        {"model", "alpha", "beta", "scale", "lower_clamp"});
    config.activated_r_migration_rate_model = strict_string(
        required(activated_r_rate, "model", "$.migration.activated_r_rate"),
        "migration.activated_r_rate.model");
    if (config.activated_r_migration_rate_model != "beta") {
        config_error(required(activated_r_rate, "model",
                              "$.migration.activated_r_rate"),
                     "migration.activated_r_rate.model", "unsupported model");
    }
    config.activated_r_migration_beta = parse_beta_rate_fields(
        activated_r_rate, "migration.activated_r_rate");

    const YAML::Node activation_duration = checked_section(
        migration, "activation_duration", "$.migration",
        {"model", "alpha", "mean_fraction", "beta"});
    if (strict_string(required(activation_duration, "model",
                               "$.migration.activation_duration"),
                      "migration.activation_duration.model") !=
        "beta_remaining_division_cycle") {
        config_error(required(activation_duration, "model",
                              "$.migration.activation_duration"),
                     "migration.activation_duration.model", "unsupported model");
    }
    config.migration_activation_duration_alpha = strict_double(
        required(activation_duration, "alpha", "$.migration.activation_duration"),
        "migration.activation_duration.alpha");
    config.migration_activation_duration_mean_fraction = strict_double(
        required(activation_duration, "mean_fraction",
                 "$.migration.activation_duration"),
        "migration.activation_duration.mean_fraction");
    config.migration_activation_duration_beta = strict_double(
        required(activation_duration, "beta", "$.migration.activation_duration"),
        "migration.activation_duration.beta");
    const YAML::Node crowding_exchange = checked_section(
        migration, "crowding_exchange", "$.migration",
        {"enabled", "stage_policy", "wait_fraction",
         "post_exchange_cooldown_fraction"});
    config.migration_swap_enabled = strict_bool(
        required(crowding_exchange, "enabled",
                 "$.migration.crowding_exchange"),
        "migration.crowding_exchange.enabled");
    config.migration_swap_stage_policy = strict_string(
        required(crowding_exchange, "stage_policy",
                 "$.migration.crowding_exchange"),
        "migration.crowding_exchange.stage_policy");
    config.migration_swap_wait_fraction = strict_double(
        required(crowding_exchange, "wait_fraction",
                 "$.migration.crowding_exchange"),
        "migration.crowding_exchange.wait_fraction");
    config.migration_swap_post_cooldown_fraction = strict_double(
        required(crowding_exchange, "post_exchange_cooldown_fraction",
                 "$.migration.crowding_exchange"),
        "migration.crowding_exchange.post_exchange_cooldown_fraction");

    const YAML::Node density = checked_section(
        root, "density", "$", {"backend", "block_edge", "growth_window_edge"});
    config.density_backend = strict_string(required(density, "backend", "$.density"),
                                           "density.backend");
    config.density_block_edge = strict_integer<int>(
        required(density, "block_edge", "$.density"), "density.block_edge");
    config.growth_density_window_edge = strict_integer<int>(
        required(density, "growth_window_edge", "$.density"),
        "density.growth_window_edge");

    const YAML::Node biology = checked_section(
        root, "biology", "$", {"alpha", "beta", "death_delay_model",
                                  "death_growth_rate_threshold",
                                  "r_death_delay_hours", "K_death_delay_hours",
                                  "r_to_K_conversion"});
    config.alpha = strict_double(required(biology, "alpha", "$.biology"), "biology.alpha");
    config.beta = strict_double(required(biology, "beta", "$.biology"), "biology.beta");
    config.death_delay_model = strict_string(
        required(biology, "death_delay_model", "$.biology"),
        "biology.death_delay_model");
    config.death_growth_rate_threshold = strict_double(
        required(biology, "death_growth_rate_threshold", "$.biology"),
        "biology.death_growth_rate_threshold");
    config.r_death_delay_hours = strict_double(
        required(biology, "r_death_delay_hours", "$.biology"),
        "biology.r_death_delay_hours");
    config.K_death_delay_hours = strict_double(
        required(biology, "K_death_delay_hours", "$.biology"),
        "biology.K_death_delay_hours");
    const YAML::Node conversion = checked_section(
        biology, "r_to_K_conversion", "$.biology",
        {"enabled", "density_window_edge", "query_block_edge",
         "density_threshold", "probability_per_division"});
    config.r_to_K_conversion.enabled = strict_bool(
        required(conversion, "enabled", "$.biology.r_to_K_conversion"),
        "biology.r_to_K_conversion.enabled");
    config.r_to_K_conversion.density_window_edge = strict_integer<int>(
        required(conversion, "density_window_edge", "$.biology.r_to_K_conversion"),
        "biology.r_to_K_conversion.density_window_edge");
    config.r_to_K_conversion.query_block_edge = strict_integer<int>(
        required(conversion, "query_block_edge", "$.biology.r_to_K_conversion"),
        "biology.r_to_K_conversion.query_block_edge");
    config.r_to_K_conversion.density_threshold = strict_double(
        required(conversion, "density_threshold", "$.biology.r_to_K_conversion"),
        "biology.r_to_K_conversion.density_threshold");
    config.r_to_K_conversion.probability_per_division = strict_double(
        required(conversion, "probability_per_division",
                 "$.biology.r_to_K_conversion"),
        "biology.r_to_K_conversion.probability_per_division");

    const YAML::Node stage = checked_section(
        root, "stage", "$", {"large_footprint_edge", "small_footprint_voxels",
                               "ultrasmall_enabled"});
    config.large_footprint_edge = strict_integer<int>(
        required(stage, "large_footprint_edge", "$.stage"), "stage.large_footprint_edge");
    config.small_footprint_voxels = strict_integer<int>(
        required(stage, "small_footprint_voxels", "$.stage"),
        "stage.small_footprint_voxels");
    config.ultrasmall_enabled = strict_bool(
        required(stage, "ultrasmall_enabled", "$.stage"), "stage.ultrasmall_enabled");

    const YAML::Node division = checked_section(
        root, "division", "$", {"shell_radius", "allow_shape_reduction", "timing"});
    config.division_shell_radius = strict_integer<int>(
        required(division, "shell_radius", "$.division"), "division.shell_radius");
    config.allow_shape_reduction = strict_bool(
        required(division, "allow_shape_reduction", "$.division"),
        "division.allow_shape_reduction");
    const YAML::Node timing = checked_section(
        division, "timing", "$.division",
        {"base_cycle_hours", "minimum_fraction", "stochastic_tail_fraction",
         "stochastic_time_quantum_hours", "retry_delay_hours",
         "inherited_growth_multiplier_min", "inherited_growth_multiplier_max",
         "r_max_inherent_growth_rate", "K_max_inherent_growth_rate"});
    config.division_timing.base_cycle_hours = strict_double(
        required(timing, "base_cycle_hours", "$.division.timing"),
        "division.timing.base_cycle_hours");
    config.division_timing.minimum_fraction = strict_double(
        required(timing, "minimum_fraction", "$.division.timing"),
        "division.timing.minimum_fraction");
    config.division_timing.stochastic_tail_fraction = strict_double(
        required(timing, "stochastic_tail_fraction", "$.division.timing"),
        "division.timing.stochastic_tail_fraction");
    config.division_timing.stochastic_time_quantum_hours = strict_double(
        required(timing, "stochastic_time_quantum_hours", "$.division.timing"),
        "division.timing.stochastic_time_quantum_hours");
    config.division_timing.retry_delay_hours = strict_double(
        required(timing, "retry_delay_hours", "$.division.timing"),
        "division.timing.retry_delay_hours");
    config.division_timing.inherited_growth_multiplier_min = strict_double(
        required(timing, "inherited_growth_multiplier_min", "$.division.timing"),
        "division.timing.inherited_growth_multiplier_min");
    config.division_timing.inherited_growth_multiplier_max = strict_double(
        required(timing, "inherited_growth_multiplier_max", "$.division.timing"),
        "division.timing.inherited_growth_multiplier_max");
    config.division_timing.r_max_inherent_growth_rate = strict_double(
        required(timing, "r_max_inherent_growth_rate", "$.division.timing"),
        "division.timing.r_max_inherent_growth_rate");
    config.division_timing.K_max_inherent_growth_rate = strict_double(
        required(timing, "K_max_inherent_growth_rate", "$.division.timing"),
        "division.timing.K_max_inherent_growth_rate");

    const YAML::Node initial = checked_section(
        root, "initial", "$", {"mode", "explicit_counts", "geometry", "r_fraction",
                                 "large_fraction", "growth_rate", "migration_rate"});
    config.initialization_mode = strict_string(required(initial, "mode", "$.initial"),
                                               "initial.mode");
    const YAML::Node explicit_counts = checked_section(
        initial, "explicit_counts", "$.initial", {"r_cells", "K_cells"});
    config.initial_r_cells = strict_integer<std::uint64_t>(
        required(explicit_counts, "r_cells", "$.initial.explicit_counts"),
        "initial.explicit_counts.r_cells");
    config.initial_K_cells = strict_integer<std::uint64_t>(
        required(explicit_counts, "K_cells", "$.initial.explicit_counts"),
        "initial.explicit_counts.K_cells");
    const YAML::Node geometry = checked_section(
        initial, "geometry", "$.initial",
        {"outer_radius", "shell_inner_radius", "inner_small_radius"});
    config.initial_radius = strict_integer<int>(
        required(geometry, "outer_radius", "$.initial.geometry"),
        "initial.geometry.outer_radius");
    config.initial_shell_inner_radius = strict_integer<int>(
        required(geometry, "shell_inner_radius", "$.initial.geometry"),
        "initial.geometry.shell_inner_radius");
    config.initial_inner_small_radius = strict_integer<int>(
        required(geometry, "inner_small_radius", "$.initial.geometry"),
        "initial.geometry.inner_small_radius");
    config.initial_shell_thickness = config.initial_radius - config.initial_shell_inner_radius;
    config.initial_r_fraction = strict_double(
        required(initial, "r_fraction", "$.initial"), "initial.r_fraction");
    config.initial_large_fraction = strict_double(
        required(initial, "large_fraction", "$.initial"), "initial.large_fraction");
    const YAML::Node growth_rate = required(initial, "growth_rate", "$.initial");
    if (!growth_rate || !growth_rate.IsMap()) {
        config_error(growth_rate, "$.initial.growth_rate", "must be a mapping");
    }
    config.initial_growth_rate_model = strict_string(
        required(growth_rate, "model", "$.initial.growth_rate"),
        "initial.growth_rate.model");
    if (config.initial_growth_rate_model == "fixed") {
        check_map(growth_rate, "$.initial.growth_rate", {"model", "fixed"});
        const YAML::Node fixed = checked_section(
            growth_rate, "fixed", "$.initial.growth_rate", {"r", "K"});
        config.initial_r_growth_rate = strict_double(
            required(fixed, "r", "$.initial.growth_rate.fixed"),
            "initial.growth_rate.fixed.r");
        config.initial_K_growth_rate = strict_double(
            required(fixed, "K", "$.initial.growth_rate.fixed"),
            "initial.growth_rate.fixed.K");
    } else if (config.initial_growth_rate_model == "legacy_truncated_normal_v1") {
        check_map(growth_rate, "$.initial.growth_rate", {"model", "truncated_normal"});
        const YAML::Node truncated_normal = checked_section(
            growth_rate, "truncated_normal", "$.initial.growth_rate", {"r", "K"});
        config.initial_r_growth_truncated_normal = parse_truncated_normal_rate(
            truncated_normal, "r", "$.initial.growth_rate.truncated_normal");
        config.initial_K_growth_truncated_normal = parse_truncated_normal_rate(
            truncated_normal, "K", "$.initial.growth_rate.truncated_normal");
    } else {
        config_error(required(growth_rate, "model", "$.initial.growth_rate"),
                     "initial.growth_rate.model", "unsupported model");
    }

    const YAML::Node migration_rate = checked_section(
        initial, "migration_rate", "$.initial", {"K"});
    const YAML::Node K_migration_rate = required(
        migration_rate, "K", "$.initial.migration_rate");
    if (!K_migration_rate || !K_migration_rate.IsMap()) {
        config_error(K_migration_rate, "$.initial.migration_rate.K",
                     "must be a mapping");
    }
    config.initial_K_migration_rate_model = strict_string(
        required(K_migration_rate, "model", "$.initial.migration_rate.K"),
        "initial.migration_rate.K.model");
    if (config.initial_K_migration_rate_model == "fixed") {
        check_map(K_migration_rate, "$.initial.migration_rate.K",
                  {"model", "value"});
        config.initial_K_migration_rate = strict_double(
            required(K_migration_rate, "value", "$.initial.migration_rate.K"),
            "initial.migration_rate.K.value");
    } else if (config.initial_K_migration_rate_model == "legacy_beta_v1") {
        check_map(K_migration_rate, "$.initial.migration_rate.K",
                  {"model", "beta"});
        config.initial_K_migration_beta = parse_beta_rate(
            K_migration_rate, "beta", "$.initial.migration_rate.K");
    } else {
        config_error(required(K_migration_rate, "model",
                              "$.initial.migration_rate.K"),
                     "initial.migration_rate.K.model", "unsupported model");
    }

    const YAML::Node simulation = checked_section(
        root, "simulation", "$", {"end_time_hours", "max_events", "threads"});
    config.end_time_hours = strict_double(
        required(simulation, "end_time_hours", "$.simulation"),
        "simulation.end_time_hours");
    config.max_events = strict_integer<std::uint64_t>(
        required(simulation, "max_events", "$.simulation"), "simulation.max_events");
    config.threads = strict_integer<int>(required(simulation, "threads", "$.simulation"),
                                         "simulation.threads");

    const YAML::Node parallel = root["parallel"];
    if (parallel && parallel.IsDefined()) {
        check_map(parallel, "$.parallel",
                  {"mode", "min_threads", "min_events_per_thread",
                   "min_refresh_items_per_thread", "cell_thresholds"});
        config.parallel_mode = strict_string(
            required(parallel, "mode", "$.parallel"), "parallel.mode");
        config.parallel_min_threads = strict_integer<int>(
            required(parallel, "min_threads", "$.parallel"),
            "parallel.min_threads");
        config.parallel_min_events_per_thread = strict_integer<std::uint64_t>(
            required(parallel, "min_events_per_thread", "$.parallel"),
            "parallel.min_events_per_thread");
        config.parallel_min_refresh_items_per_thread =
            strict_integer<std::uint64_t>(
                required(parallel, "min_refresh_items_per_thread", "$.parallel"),
                "parallel.min_refresh_items_per_thread");
        const YAML::Node thresholds = required(
            parallel, "cell_thresholds", "$.parallel");
        if (!thresholds.IsSequence() || thresholds.size() == 0) {
            config_error(thresholds, "parallel.cell_thresholds",
                         "must be a non-empty sequence");
        }
        config.parallel_thread_thresholds.clear();
        config.parallel_thread_thresholds.reserve(thresholds.size());
        for (std::size_t index = 0; index < thresholds.size(); ++index) {
            const YAML::Node threshold = thresholds[index];
            const std::string path =
                "parallel.cell_thresholds[" + std::to_string(index) + ']';
            check_map(threshold, path, {"minimum_cells", "max_thread_fraction"});
            config.parallel_thread_thresholds.push_back({
                strict_integer<std::uint64_t>(
                    required(threshold, "minimum_cells", path),
                    path + ".minimum_cells"),
                strict_double(required(threshold, "max_thread_fraction", path),
                              path + ".max_thread_fraction")});
        }
    }

    const YAML::Node scheduler = checked_section(
        root, "scheduler", "$", {"backend",
                                   "conflict_bucket_hours",
                                   "proposal_window_hours",
                                   "proposal_window_max_events",
                                   "proposal_min_events_per_thread",
                                   "proposal_dependency_block_edge"});
    config.scheduler_backend = strict_string(
        required(scheduler, "backend", "$.scheduler"), "scheduler.backend");
    config.conflict_bucket_hours = strict_double(
        required(scheduler, "conflict_bucket_hours", "$.scheduler"),
        "scheduler.conflict_bucket_hours");
    config.proposal_window_hours = strict_double(
        required(scheduler, "proposal_window_hours", "$.scheduler"),
        "scheduler.proposal_window_hours");
    config.proposal_window_max_events = strict_integer<std::uint64_t>(
        required(scheduler, "proposal_window_max_events", "$.scheduler"),
        "scheduler.proposal_window_max_events");
    config.proposal_min_events_per_thread = strict_integer<std::uint64_t>(
        required(scheduler, "proposal_min_events_per_thread", "$.scheduler"),
        "scheduler.proposal_min_events_per_thread");
    config.proposal_dependency_block_edge = strict_integer<int>(
        required(scheduler, "proposal_dependency_block_edge", "$.scheduler"),
        "scheduler.proposal_dependency_block_edge");

    const YAML::Node angiogenesis = checked_section(
        root, "angiogenesis", "$", {"enabled", "lesion_detection", "trigger",
                                     "seed_process", "vessel", "direction",
                                     "influence"});
    config.angiogenesis.enabled = strict_bool(
        required(angiogenesis, "enabled", "$.angiogenesis"), "angiogenesis.enabled");
    const YAML::Node lesion_detection = checked_section(
        angiogenesis, "lesion_detection", "$.angiogenesis",
        {"backend", "block_edge", "connectivity",
         "core_activation_occupied_fraction",
         "core_deactivation_occupied_fraction",
         "minimum_cells_per_core_block",
         "minimum_biological_volume_per_core_block", "halo_blocks",
         "refresh_interval_hours"});
    config.angiogenesis.lesion_detection_backend = strict_string(
        required(lesion_detection, "backend", "$.angiogenesis.lesion_detection"),
        "angiogenesis.lesion_detection.backend");
    config.angiogenesis.lesion_block_edge = strict_integer<int>(
        required(lesion_detection, "block_edge", "$.angiogenesis.lesion_detection"),
        "angiogenesis.lesion_detection.block_edge");
    config.angiogenesis.lesion_connectivity = strict_integer<int>(
        required(lesion_detection, "connectivity", "$.angiogenesis.lesion_detection"),
        "angiogenesis.lesion_detection.connectivity");
    config.angiogenesis.lesion_core_activation_occupied_fraction = strict_double(
        required(lesion_detection, "core_activation_occupied_fraction",
                 "$.angiogenesis.lesion_detection"),
        "angiogenesis.lesion_detection.core_activation_occupied_fraction");
    config.angiogenesis.lesion_core_deactivation_occupied_fraction = strict_double(
        required(lesion_detection, "core_deactivation_occupied_fraction",
                 "$.angiogenesis.lesion_detection"),
        "angiogenesis.lesion_detection.core_deactivation_occupied_fraction");
    config.angiogenesis.lesion_minimum_cells_per_core_block =
        strict_integer<std::uint64_t>(
            required(lesion_detection, "minimum_cells_per_core_block",
                     "$.angiogenesis.lesion_detection"),
            "angiogenesis.lesion_detection.minimum_cells_per_core_block");
    config.angiogenesis.lesion_minimum_biological_volume_per_core_block =
        strict_double(
            required(lesion_detection,
                     "minimum_biological_volume_per_core_block",
                     "$.angiogenesis.lesion_detection"),
            "angiogenesis.lesion_detection.minimum_biological_volume_per_core_block");
    config.angiogenesis.lesion_halo_blocks = strict_integer<int>(
        required(lesion_detection, "halo_blocks",
                 "$.angiogenesis.lesion_detection"),
        "angiogenesis.lesion_detection.halo_blocks");
    config.angiogenesis.lesion_refresh_interval_hours = strict_double(
        required(lesion_detection, "refresh_interval_hours",
                 "$.angiogenesis.lesion_detection"),
        "angiogenesis.lesion_detection.refresh_interval_hours");
    const YAML::Node trigger = checked_section(
        angiogenesis, "trigger", "$.angiogenesis",
        {"metric", "activation_volume_voxels3", "deactivation_volume_voxels3",
         "delay_hours", "minimum_core_blocks", "stage_biological_volume_voxels3"});
    config.angiogenesis.trigger_metric = strict_string(
        required(trigger, "metric", "$.angiogenesis.trigger"),
        "angiogenesis.trigger.metric");
    config.angiogenesis.trigger_activation_volume_voxels3 = strict_double(
        required(trigger, "activation_volume_voxels3", "$.angiogenesis.trigger"),
        "angiogenesis.trigger.activation_volume_voxels3");
    config.angiogenesis.trigger_deactivation_volume_voxels3 = strict_double(
        required(trigger, "deactivation_volume_voxels3", "$.angiogenesis.trigger"),
        "angiogenesis.trigger.deactivation_volume_voxels3");
    config.angiogenesis.trigger_delay_hours = strict_double(
        required(trigger, "delay_hours", "$.angiogenesis.trigger"),
        "angiogenesis.trigger.delay_hours");
    config.angiogenesis.trigger_minimum_core_blocks = strict_integer<std::uint64_t>(
        required(trigger, "minimum_core_blocks", "$.angiogenesis.trigger"),
        "angiogenesis.trigger.minimum_core_blocks");
    const YAML::Node stage_volume = checked_section(
        trigger, "stage_biological_volume_voxels3", "$.angiogenesis.trigger",
        {"stage0", "stage1", "stage2"});
    config.angiogenesis.stage0_biological_volume_voxels3 = strict_double(
        required(stage_volume, "stage0", "$.angiogenesis.trigger.stage_biological_volume_voxels3"),
        "angiogenesis.trigger.stage_biological_volume_voxels3.stage0");
    config.angiogenesis.stage1_biological_volume_voxels3 = strict_double(
        required(stage_volume, "stage1", "$.angiogenesis.trigger.stage_biological_volume_voxels3"),
        "angiogenesis.trigger.stage_biological_volume_voxels3.stage1");
    config.angiogenesis.stage2_biological_volume_voxels3 = strict_double(
        required(stage_volume, "stage2", "$.angiogenesis.trigger.stage_biological_volume_voxels3"),
        "angiogenesis.trigger.stage_biological_volume_voxels3.stage2");

    const YAML::Node seed_process = checked_section(
        angiogenesis, "seed_process", "$.angiogenesis",
        {"model", "scope", "rate_sites_per_30_days", "roots_per_event",
         "density_stress_on_fraction", "density_stress_full_fraction",
         "density_stress_exponent", "volume_reference_voxels3",
         "volume_exponent", "minimum_rate_multiplier",
         "maximum_rate_multiplier",
         "surface_min_separation_voxels", "surface_max_sampling_attempts",
         "surface_min_local_cells", "root_position_policy",
         "max_total_roots", "max_active_tips", "max_roots_per_lesion",
         "max_active_tips_per_lesion"});
    config.angiogenesis.seed_process_model = strict_string(
        required(seed_process, "model", "$.angiogenesis.seed_process"),
        "angiogenesis.seed_process.model");
    config.angiogenesis.seed_process_scope = strict_string(
        required(seed_process, "scope", "$.angiogenesis.seed_process"),
        "angiogenesis.seed_process.scope");
    config.angiogenesis.seed_rate_sites_per_30_days = strict_double(
        required(seed_process, "rate_sites_per_30_days", "$.angiogenesis.seed_process"),
        "angiogenesis.seed_process.rate_sites_per_30_days");
    config.angiogenesis.seed_rate_sites_per_hour =
        config.angiogenesis.seed_rate_sites_per_30_days / 720.0;
    config.angiogenesis.seed_density_stress_on_fraction = strict_double(
        required(seed_process, "density_stress_on_fraction",
                 "$.angiogenesis.seed_process"),
        "angiogenesis.seed_process.density_stress_on_fraction");
    config.angiogenesis.seed_density_stress_full_fraction = strict_double(
        required(seed_process, "density_stress_full_fraction",
                 "$.angiogenesis.seed_process"),
        "angiogenesis.seed_process.density_stress_full_fraction");
    config.angiogenesis.seed_density_stress_exponent = strict_double(
        required(seed_process, "density_stress_exponent",
                 "$.angiogenesis.seed_process"),
        "angiogenesis.seed_process.density_stress_exponent");
    config.angiogenesis.seed_volume_reference_voxels3 = strict_double(
        required(seed_process, "volume_reference_voxels3",
                 "$.angiogenesis.seed_process"),
        "angiogenesis.seed_process.volume_reference_voxels3");
    config.angiogenesis.seed_volume_exponent = strict_double(
        required(seed_process, "volume_exponent", "$.angiogenesis.seed_process"),
        "angiogenesis.seed_process.volume_exponent");
    config.angiogenesis.seed_minimum_rate_multiplier = strict_double(
        required(seed_process, "minimum_rate_multiplier",
                 "$.angiogenesis.seed_process"),
        "angiogenesis.seed_process.minimum_rate_multiplier");
    config.angiogenesis.seed_maximum_rate_multiplier = strict_double(
        required(seed_process, "maximum_rate_multiplier",
                 "$.angiogenesis.seed_process"),
        "angiogenesis.seed_process.maximum_rate_multiplier");
    config.angiogenesis.roots_per_event = strict_integer<std::uint32_t>(
        required(seed_process, "roots_per_event", "$.angiogenesis.seed_process"),
        "angiogenesis.seed_process.roots_per_event");
    config.angiogenesis.surface_min_separation_voxels = strict_integer<int>(
        required(seed_process, "surface_min_separation_voxels", "$.angiogenesis.seed_process"),
        "angiogenesis.seed_process.surface_min_separation_voxels");
    config.angiogenesis.surface_max_sampling_attempts = strict_integer<std::uint32_t>(
        required(seed_process, "surface_max_sampling_attempts", "$.angiogenesis.seed_process"),
        "angiogenesis.seed_process.surface_max_sampling_attempts");
    config.angiogenesis.surface_min_local_cells = strict_integer<std::uint64_t>(
        required(seed_process, "surface_min_local_cells", "$.angiogenesis.seed_process"),
        "angiogenesis.seed_process.surface_min_local_cells");
    config.angiogenesis.root_position_policy = strict_string(
        required(seed_process, "root_position_policy", "$.angiogenesis.seed_process"),
        "angiogenesis.seed_process.root_position_policy");
    config.angiogenesis.max_total_roots = strict_integer<std::uint64_t>(
        required(seed_process, "max_total_roots", "$.angiogenesis.seed_process"),
        "angiogenesis.seed_process.max_total_roots");
    config.angiogenesis.max_active_tips = strict_integer<std::uint64_t>(
        required(seed_process, "max_active_tips", "$.angiogenesis.seed_process"),
        "angiogenesis.seed_process.max_active_tips");
    config.angiogenesis.max_roots_per_lesion = strict_integer<std::uint64_t>(
        required(seed_process, "max_roots_per_lesion", "$.angiogenesis.seed_process"),
        "angiogenesis.seed_process.max_roots_per_lesion");
    config.angiogenesis.max_active_tips_per_lesion = strict_integer<std::uint64_t>(
        required(seed_process, "max_active_tips_per_lesion",
                 "$.angiogenesis.seed_process"),
        "angiogenesis.seed_process.max_active_tips_per_lesion");

    const YAML::Node vessel = checked_section(
        angiogenesis, "vessel", "$.angiogenesis",
        {"diameter_voxels", "inward", "outward", "contact",
         "blocked_policy", "blocked_retry_interval_hours",
         "collision_policy", "boundary_policy"});
    config.angiogenesis.diameter_voxels = strict_double(
        required(vessel, "diameter_voxels", "$.angiogenesis.vessel"),
        "angiogenesis.vessel.diameter_voxels");
    const YAML::Node inward = checked_section(
        vessel, "inward", "$.angiogenesis.vessel",
        {"speed_voxels_per_hour", "max_length_voxels",
         "length_tortuosity_factor", "exit_margin_voxels",
         "hard_max_length_voxels", "target_tolerance_voxels",
         "replacement_policy", "path_policy", "far_surface_policy"});
    config.angiogenesis.inward_speed_voxels_per_hour = strict_double(
        required(inward, "speed_voxels_per_hour", "$.angiogenesis.vessel.inward"),
        "angiogenesis.vessel.inward.speed_voxels_per_hour");
    config.angiogenesis.inward_max_length_voxels = strict_integer<int>(
        required(inward, "max_length_voxels", "$.angiogenesis.vessel.inward"),
        "angiogenesis.vessel.inward.max_length_voxels");
    config.angiogenesis.inward_length_tortuosity_factor = strict_double(
        required(inward, "length_tortuosity_factor",
                 "$.angiogenesis.vessel.inward"),
        "angiogenesis.vessel.inward.length_tortuosity_factor");
    config.angiogenesis.inward_exit_margin_voxels = strict_double(
        required(inward, "exit_margin_voxels", "$.angiogenesis.vessel.inward"),
        "angiogenesis.vessel.inward.exit_margin_voxels");
    config.angiogenesis.inward_hard_max_length_voxels =
        strict_integer<int>(
            required(inward, "hard_max_length_voxels",
                     "$.angiogenesis.vessel.inward"),
            "angiogenesis.vessel.inward.hard_max_length_voxels");
    config.angiogenesis.inward_target_tolerance_voxels = strict_double(
        required(inward, "target_tolerance_voxels", "$.angiogenesis.vessel.inward"),
        "angiogenesis.vessel.inward.target_tolerance_voxels");
    config.angiogenesis.inward_replacement_policy = strict_string(
        required(inward, "replacement_policy", "$.angiogenesis.vessel.inward"),
        "angiogenesis.vessel.inward.replacement_policy");
    config.angiogenesis.inward_path_policy = strict_string(
        required(inward, "path_policy", "$.angiogenesis.vessel.inward"),
        "angiogenesis.vessel.inward.path_policy");
    config.angiogenesis.inward_far_surface_policy = strict_string(
        required(inward, "far_surface_policy", "$.angiogenesis.vessel.inward"),
        "angiogenesis.vessel.inward.far_surface_policy");
    const YAML::Node outward = checked_section(
        vessel, "outward", "$.angiogenesis.vessel",
        {"speed_voxels_per_hour", "max_length_voxels",
         "external_connection_distance_voxels", "occupancy_policy"});
    config.angiogenesis.outward_speed_voxels_per_hour = strict_double(
        required(outward, "speed_voxels_per_hour", "$.angiogenesis.vessel.outward"),
        "angiogenesis.vessel.outward.speed_voxels_per_hour");
    config.angiogenesis.outward_max_length_voxels = strict_integer<int>(
        required(outward, "max_length_voxels", "$.angiogenesis.vessel.outward"),
        "angiogenesis.vessel.outward.max_length_voxels");
    config.angiogenesis.outward_external_connection_distance_voxels = strict_double(
        required(outward, "external_connection_distance_voxels", "$.angiogenesis.vessel.outward"),
        "angiogenesis.vessel.outward.external_connection_distance_voxels");
    config.angiogenesis.outward_occupancy_policy = strict_string(
        required(outward, "occupancy_policy", "$.angiogenesis.vessel.outward"),
        "angiogenesis.vessel.outward.occupancy_policy");
    const YAML::Node contact = checked_section(
        vessel, "contact", "$.angiogenesis.vessel",
        {"outward_same_lesion", "outward_other_lesion",
         "inward_other_lesion"});
    config.angiogenesis.outward_same_lesion_contact_policy = strict_string(
        required(contact, "outward_same_lesion", "$.angiogenesis.vessel.contact"),
        "angiogenesis.vessel.contact.outward_same_lesion");
    config.angiogenesis.outward_other_lesion_contact_policy = strict_string(
        required(contact, "outward_other_lesion", "$.angiogenesis.vessel.contact"),
        "angiogenesis.vessel.contact.outward_other_lesion");
    config.angiogenesis.inward_other_lesion_contact_policy = strict_string(
        required(contact, "inward_other_lesion", "$.angiogenesis.vessel.contact"),
        "angiogenesis.vessel.contact.inward_other_lesion");
    config.angiogenesis.vessel_blocked_policy = strict_string(
        required(vessel, "blocked_policy", "$.angiogenesis.vessel"),
        "angiogenesis.vessel.blocked_policy");
    config.angiogenesis.vessel_blocked_retry_interval_hours = strict_double(
        required(vessel, "blocked_retry_interval_hours", "$.angiogenesis.vessel"),
        "angiogenesis.vessel.blocked_retry_interval_hours");
    config.angiogenesis.vessel_collision_policy = strict_string(
        required(vessel, "collision_policy", "$.angiogenesis.vessel"),
        "angiogenesis.vessel.collision_policy");
    config.angiogenesis.boundary_policy = strict_string(
        required(vessel, "boundary_policy", "$.angiogenesis.vessel"),
        "angiogenesis.vessel.boundary_policy");

    const YAML::Node vessel_direction = checked_section(
        angiogenesis, "direction", "$.angiogenesis",
        {"model", "forward_bias", "forward_half_angle_deg", "turn_half_angle_deg",
         "persistence_probability", "distance_weight_exponent"});
    config.angiogenesis.direction_model = strict_string(
        required(vessel_direction, "model", "$.angiogenesis.direction"),
        "angiogenesis.direction.model");
    config.angiogenesis.direction_forward_bias = strict_double(
        required(vessel_direction, "forward_bias", "$.angiogenesis.direction"),
        "angiogenesis.direction.forward_bias");
    config.angiogenesis.direction_half_angle_degrees = strict_double(
        required(vessel_direction, "forward_half_angle_deg", "$.angiogenesis.direction"),
        "angiogenesis.direction.forward_half_angle_deg");
    config.angiogenesis.direction_turn_half_angle_degrees = strict_double(
        required(vessel_direction, "turn_half_angle_deg", "$.angiogenesis.direction"),
        "angiogenesis.direction.turn_half_angle_deg");
    config.angiogenesis.direction_persistence_probability = strict_double(
        required(vessel_direction, "persistence_probability", "$.angiogenesis.direction"),
        "angiogenesis.direction.persistence_probability");
    config.angiogenesis.direction_distance_weight_exponent = strict_double(
        required(vessel_direction, "distance_weight_exponent", "$.angiogenesis.direction"),
        "angiogenesis.direction.distance_weight_exponent");

    const YAML::Node influence = checked_section(
        angiogenesis, "influence", "$.angiogenesis",
        {"profile", "max_relief_fraction", "decay_length_voxels", "cutoff_radius_voxels",
         "scope", "activation"});
    config.angiogenesis.influence_profile = strict_string(
        required(influence, "profile", "$.angiogenesis.influence"),
        "angiogenesis.influence.profile");
    config.angiogenesis.influence_max_relief_fraction = strict_double(
        required(influence, "max_relief_fraction", "$.angiogenesis.influence"),
        "angiogenesis.influence.max_relief_fraction");
    config.angiogenesis.influence_decay_length_voxels = strict_double(
        required(influence, "decay_length_voxels", "$.angiogenesis.influence"),
        "angiogenesis.influence.decay_length_voxels");
    config.angiogenesis.influence_cutoff_radius_voxels = strict_double(
        required(influence, "cutoff_radius_voxels", "$.angiogenesis.influence"),
        "angiogenesis.influence.cutoff_radius_voxels");
    config.angiogenesis.influence_scope = strict_string(
        required(influence, "scope", "$.angiogenesis.influence"),
        "angiogenesis.influence.scope");
    config.angiogenesis.influence_activation = strict_string(
        required(influence, "activation", "$.angiogenesis.influence"),
        "angiogenesis.influence.activation");

    const YAML::Node output = checked_section(
        root, "output", "$", {"enabled", "preview_mode", "full_format", "checkpoint_format",
                                "directory", "preview_every_hours", "full_every_hours",
                                "checkpoint_every_hours", "storage", "async_enabled",
                                "async_queue_depth", "async_max_pending_bytes",
                                "preview_overflow_policy", "preview_max_cells",
                                "preview_seed", "display_radius", "sampling_preset",
                                "on_demand", "live_preview"});
    config.output_enabled = strict_bool(required(output, "enabled", "$.output"), "output.enabled");
    if (const YAML::Node value = output["sampling_preset"];
        value && value.IsDefined()) {
        config.output_sampling_preset =
            strict_string(value, "output.sampling_preset");
    }
    config.preview_mode = strict_string(required(output, "preview_mode", "$.output"),
                                        "output.preview_mode");
    config.full_format = strict_string(required(output, "full_format", "$.output"),
                                       "output.full_format");
    config.checkpoint_format = strict_string(
        required(output, "checkpoint_format", "$.output"), "output.checkpoint_format");
    config.output_directory = strict_string(required(output, "directory", "$.output"),
                                            "output.directory");
    config.preview_every_hours = strict_double(
        required(output, "preview_every_hours", "$.output"), "output.preview_every_hours");
    config.full_every_hours = strict_double(
        required(output, "full_every_hours", "$.output"), "output.full_every_hours");
    config.checkpoint_every_hours = strict_double(
        required(output, "checkpoint_every_hours", "$.output"),
        "output.checkpoint_every_hours");
    const YAML::Node storage = checked_section(
        output, "storage", "$.output",
        {"mode", "vtkhdf_compression_level", "hdf5_compression_level",
         "hdf5_chunk_elements", "preview_keyframe_every_hours",
         "full_keyframe_every_hours", "checkpoint_base_every_hours",
         "checkpoint_max_delta_chain", "delta_full_ratio"});
    config.storage_mode = strict_string(
        required(storage, "mode", "$.output.storage"), "output.storage.mode");
    config.vtkhdf_compression_level = strict_integer<int>(
        required(storage, "vtkhdf_compression_level", "$.output.storage"),
        "output.storage.vtkhdf_compression_level");
    config.hdf5_compression_level = strict_integer<int>(
        required(storage, "hdf5_compression_level", "$.output.storage"),
        "output.storage.hdf5_compression_level");
    config.hdf5_chunk_elements = strict_integer<std::uint64_t>(
        required(storage, "hdf5_chunk_elements", "$.output.storage"),
        "output.storage.hdf5_chunk_elements");
    config.preview_keyframe_every_hours = strict_double(
        required(storage, "preview_keyframe_every_hours", "$.output.storage"),
        "output.storage.preview_keyframe_every_hours");
    config.full_keyframe_every_hours = strict_double(
        required(storage, "full_keyframe_every_hours", "$.output.storage"),
        "output.storage.full_keyframe_every_hours");
    config.checkpoint_base_every_hours = strict_double(
        required(storage, "checkpoint_base_every_hours", "$.output.storage"),
        "output.storage.checkpoint_base_every_hours");
    config.checkpoint_max_delta_chain = strict_integer<std::uint64_t>(
        required(storage, "checkpoint_max_delta_chain", "$.output.storage"),
        "output.storage.checkpoint_max_delta_chain");
    config.delta_full_ratio = strict_double(
        required(storage, "delta_full_ratio", "$.output.storage"),
        "output.storage.delta_full_ratio");
    config.output_async_enabled = strict_bool(
        required(output, "async_enabled", "$.output"),
        "output.async_enabled");
    config.output_async_queue_depth = strict_integer<std::uint64_t>(
        required(output, "async_queue_depth", "$.output"),
        "output.async_queue_depth");
    if (const YAML::Node value = output["async_max_pending_bytes"];
        value && value.IsDefined()) {
        config.output_async_max_pending_bytes =
            strict_integer<std::uint64_t>(
                value, "output.async_max_pending_bytes");
    }
    if (const YAML::Node value = output["preview_overflow_policy"];
        value && value.IsDefined()) {
        config.preview_overflow_policy =
            strict_string(value, "output.preview_overflow_policy");
    }
    if (const YAML::Node on_demand = output["on_demand"];
        on_demand && on_demand.IsDefined()) {
        check_map(on_demand, "$.output.on_demand",
                  {"preview", "checkpoint", "full"});
        config.output_on_demand_preview = strict_bool(
            required(on_demand, "preview", "$.output.on_demand"),
            "output.on_demand.preview");
        config.output_on_demand_checkpoint = strict_bool(
            required(on_demand, "checkpoint", "$.output.on_demand"),
            "output.on_demand.checkpoint");
        config.output_on_demand_full = strict_bool(
            required(on_demand, "full", "$.output.on_demand"),
            "output.on_demand.full");
    }
    if (const YAML::Node live = output["live_preview"];
        live && live.IsDefined()) {
        check_map(live, "$.output.live_preview",
                  {"enabled_when_viewer_attached",
                   "wall_interval_seconds", "persist"});
        config.live_preview_when_attached = strict_bool(
            required(live, "enabled_when_viewer_attached",
                     "$.output.live_preview"),
            "output.live_preview.enabled_when_viewer_attached");
        config.live_preview_wall_interval_seconds = strict_double(
            required(live, "wall_interval_seconds",
                     "$.output.live_preview"),
            "output.live_preview.wall_interval_seconds");
        config.live_preview_persist = strict_bool(
            required(live, "persist", "$.output.live_preview"),
            "output.live_preview.persist");
    }
    config.preview_max_cells = strict_integer<std::uint64_t>(
        required(output, "preview_max_cells", "$.output"), "output.preview_max_cells");
    config.preview_seed = strict_integer<std::uint64_t>(
        required(output, "preview_seed", "$.output"), "output.preview_seed");
    const YAML::Node display_radius = checked_section(
        output, "display_radius", "$.output", {"large", "small", "ultrasmall"});
    config.display_radius.large = strict_float(
        required(display_radius, "large", "$.output.display_radius"),
        "output.display_radius.large");
    config.display_radius.small = strict_float(
        required(display_radius, "small", "$.output.display_radius"),
        "output.display_radius.small");
    config.display_radius.ultrasmall = strict_float(
        required(display_radius, "ultrasmall", "$.output.display_radius"),
        "output.display_radius.ultrasmall");

    if (const YAML::Node control = root["control"];
        control && control.IsDefined()) {
        check_map(control, "$.control",
                  {"enabled", "status_wall_interval_seconds",
                   "poll_wall_interval_seconds"});
        config.control_enabled = strict_bool(
            required(control, "enabled", "$.control"),
            "control.enabled");
        config.control_status_wall_interval_seconds = strict_double(
            required(control, "status_wall_interval_seconds",
                     "$.control"),
            "control.status_wall_interval_seconds");
        config.control_poll_wall_interval_seconds = strict_double(
            required(control, "poll_wall_interval_seconds",
                     "$.control"),
            "control.poll_wall_interval_seconds");
    }

    config.validate();
    return config;
}

}  // namespace atcg3d
