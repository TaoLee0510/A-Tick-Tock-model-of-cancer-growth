#include "config/structured_config.hpp"
#include "model/division_renewal.hpp"
#include "model/beta_duration.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

namespace atcg3d::structured_pde {
namespace {

[[noreturn]] void fail(const std::string& path, const std::string& message) {
    throw std::invalid_argument("structured PDE config " + path + ": " + message);
}

YAML::Node required(const YAML::Node& parent,
                    const char* key,
                    const std::string& path) {
    const YAML::Node result = parent[key];
    if (!result) fail(path + "." + key, "is required");
    return result;
}

void mapping(const YAML::Node& node,
             const std::string& path,
             std::initializer_list<const char*> allowed) {
    if (!node.IsMap()) fail(path, "must be a mapping");
    std::set<std::string> accepted;
    for (const char* key : allowed) accepted.emplace(key);
    for (const auto& item : node) {
        const std::string key = item.first.as<std::string>();
        if (!accepted.contains(key)) fail(path + "." + key, "is unknown");
    }
}

std::string text(const YAML::Node& node, const std::string& path) {
    if (!node.IsScalar()) fail(path, "must be a scalar string");
    try {
        return node.as<std::string>();
    } catch (const std::exception&) {
        fail(path, "must be a string");
    }
}

double number(const YAML::Node& node, const std::string& path) {
    if (!node.IsScalar()) fail(path, "must be a finite number");
    double result{};
    try {
        result = node.as<double>();
    } catch (const std::exception&) {
        fail(path, "must be a finite number");
    }
    if (!std::isfinite(result)) fail(path, "must be finite");
    return result;
}

int integer(const YAML::Node& node, const std::string& path) {
    if (!node.IsScalar()) fail(path, "must be an integer");
    try {
        return node.as<int>();
    } catch (const std::exception&) {
        fail(path, "must be an integer");
    }
}

bool boolean(const YAML::Node& node, const std::string& path) {
    if (!node.IsScalar()) fail(path, "must be a boolean");
    try {
        return node.as<bool>();
    } catch (const std::exception&) {
        fail(path, "must be a boolean");
    }
}

bool close(double lhs, double rhs) noexcept {
    return std::abs(lhs - rhs) <=
        1.0e-12 * std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

std::uint64_t mix(std::uint64_t state, std::uint64_t value) noexcept {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    value ^= value >> 31U;
    return state ^ (value + (state << 6U) + (state >> 2U));
}

void hash_text(std::uint64_t& state, const std::string& value) noexcept {
    for (const unsigned char character : value) state = mix(state, character);
    state = mix(state, value.size());
}

std::string escaped(const std::string& value) {
    std::string result;
    for (const char character : value) {
        if (character == '\\' || character == '"') result.push_back('\\');
        result.push_back(character);
    }
    return result;
}

}  // namespace

StructuredPdeConfig3D StructuredPdeConfig3D::load(
    const std::filesystem::path& path) {
    YAML::Node root;
    try {
        root = YAML::LoadFile(path.string());
    } catch (const std::exception& error) {
        throw std::invalid_argument(
            "unable to load structured PDE config " + path.string() + ": " +
            error.what());
    }
    mapping(root, "$", {"schema", "profile", "continuum_config",
                         "structured_migration", "output", "storage", "division_clock"});
    const YAML::Node schema = required(root, "schema", "$");
    mapping(schema, "$.schema", {"name", "version"});
    if (text(required(schema, "name", "$.schema"), "$.schema.name") !=
        "atcg3d.structured_pde_config") {
        fail("$.schema.name", "must equal atcg3d.structured_pde_config");
    }

    StructuredPdeConfig3D result;
    result.schema_version = integer(
        required(schema, "version", "$.schema"), "$.schema.version");
    result.profile = text(required(root, "profile", "$"), "$.profile");
    result.source_path = std::filesystem::absolute(path).lexically_normal();
    const auto relative = std::filesystem::path(text(
        required(root, "continuum_config", "$"), "$.continuum_config"));
    result.continuum_config_path = relative.is_absolute()
        ? relative.lexically_normal()
        : (result.source_path.parent_path() / relative).lexically_normal();
    result.continuum =
        continuum::ContinuumModelConfig3D::load(result.continuum_config_path);
    if (const YAML::Node output = root["output"]) {
        mapping(output, "$.output", {"directory"});
        result.continuum.output.directory = text(
            required(output, "directory", "$.output"), "$.output.directory");
        result.continuum.validate();
    }

    const YAML::Node migration =
        required(root, "structured_migration", "$");
    mapping(migration, "$.structured_migration",
            {"model", "activation_density", "activation_clock",
             "activation_time_bin_width_hours", "activation_maximum_hours",
             "activation_rate_model", "activation_rate_bins",
             "activation_stop", "direction_transport",
             "direction_density_window_edge", "direction_nutrient_window_edge",
             "chemotaxis_strength", "zero_gradient_tolerance",
             "reactivation_cooldown_hours",
             "reactivation_density_threshold",
             "crowding_exchange", "vessel_exclusion",
             "maximum_move_probability_per_substep", "minimum_density"});
    result.migration.model = text(
        required(migration, "model", "$.structured_migration"),
        "$.structured_migration.model");
    result.migration.activation_density = text(
        required(migration, "activation_density", "$.structured_migration"),
        "$.structured_migration.activation_density");
    result.migration.activation_clock = text(
        required(migration, "activation_clock", "$.structured_migration"),
        "$.structured_migration.activation_clock");
    if (migration["activation_time_bin_width_hours"] || migration["activation_maximum_hours"]) {
        if (result.schema_version < 11)
            fail("$.structured_migration.activation_clock", "duration grid requires schema v11");
        result.migration.activation_time_bin_width_hours = number(
            required(migration, "activation_time_bin_width_hours", "$.structured_migration"),
            "$.structured_migration.activation_time_bin_width_hours");
        result.migration.activation_maximum_hours = number(
            required(migration, "activation_maximum_hours", "$.structured_migration"),
            "$.structured_migration.activation_maximum_hours");
    }
    if (migration["activation_rate_model"] || migration["activation_rate_bins"]) {
        if (result.schema_version < 12)
            fail("$.structured_migration.activation_rate_model", "requires schema v12");
        result.migration.activation_rate_model = text(
            required(migration, "activation_rate_model", "$.structured_migration"),
            "$.structured_migration.activation_rate_model");
        result.migration.activation_rate_bins = integer(
            required(migration, "activation_rate_bins", "$.structured_migration"),
            "$.structured_migration.activation_rate_bins");
    }
    if (result.schema_version >= 5) {
        result.migration.activation_stop = text(
            required(migration, "activation_stop", "$.structured_migration"),
            "$.structured_migration.activation_stop");
    }
    result.migration.direction_transport = text(
        required(migration, "direction_transport", "$.structured_migration"),
        "$.structured_migration.direction_transport");
    if (result.schema_version >= 3 && result.schema_version < 5) {
        result.migration.direction_density_window_edge = integer(
            required(migration, "direction_density_window_edge",
                     "$.structured_migration"),
            "$.structured_migration.direction_density_window_edge");
    }
    if (result.schema_version >= 5) {
        result.migration.direction_nutrient_window_edge = integer(
            required(migration, "direction_nutrient_window_edge",
                     "$.structured_migration"),
            "$.structured_migration.direction_nutrient_window_edge");
        result.migration.chemotaxis_strength = number(
            required(migration, "chemotaxis_strength",
                     "$.structured_migration"),
            "$.structured_migration.chemotaxis_strength");
        result.migration.zero_gradient_tolerance = number(
            required(migration, "zero_gradient_tolerance",
                     "$.structured_migration"),
            "$.structured_migration.zero_gradient_tolerance");
        result.migration.reactivation_cooldown_hours = number(
            required(migration, "reactivation_cooldown_hours",
                     "$.structured_migration"),
            "$.structured_migration.reactivation_cooldown_hours");
        result.migration.reactivation_density_threshold = number(
            required(migration, "reactivation_density_threshold",
                     "$.structured_migration"),
            "$.structured_migration.reactivation_density_threshold");
    }
    if (result.schema_version >= 3) {
        result.migration.crowding_exchange = text(
            required(migration, "crowding_exchange",
                     "$.structured_migration"),
            "$.structured_migration.crowding_exchange");
        result.migration.vessel_exclusion = boolean(
            required(migration, "vessel_exclusion",
                     "$.structured_migration"),
            "$.structured_migration.vessel_exclusion");
    }
    result.migration.maximum_move_probability_per_substep = number(
        required(migration, "maximum_move_probability_per_substep",
                 "$.structured_migration"),
        "$.structured_migration.maximum_move_probability_per_substep");
    result.migration.minimum_density = number(
        required(migration, "minimum_density", "$.structured_migration"),
        "$.structured_migration.minimum_density");
    if (root["storage"]) {
        if (result.schema_version < 9) fail("$.storage", "requires schema v9");
        const auto storage = root["storage"];
        mapping(storage, "$.storage", {"model", "maximum_active_voxels"});
        result.storage_model = text(required(storage, "model", "$.storage"),
            "$.storage.model");
        result.maximum_active_voxels = integer(
            required(storage, "maximum_active_voxels", "$.storage"),
            "$.storage.maximum_active_voxels");
    }
    if (root["division_clock"]) {
        if (result.schema_version < 10) fail("$.division_clock", "requires schema v10");
        const auto clock = root["division_clock"];
        mapping(clock, "$.division_clock", {"model", "work_bin_width", "maximum_work"});
        result.division_clock_model = text(required(clock, "model", "$.division_clock"), "$.division_clock.model");
        result.division_work_bin_width = number(required(clock, "work_bin_width", "$.division_clock"), "$.division_clock.work_bin_width");
        result.division_maximum_work = number(required(clock, "maximum_work", "$.division_clock"), "$.division_clock.maximum_work");
    }
    result.validate();
    return result;
}

void StructuredPdeConfig3D::validate() const {
    if ((schema_version < 1 || schema_version > 13) || profile.empty()) {
        throw std::invalid_argument("structured PDE schema/profile is invalid");
    }
    if (schema_version < 9 && storage_model != "dense_v1") {
        throw std::invalid_argument("sparse storage requires structured v9");
    }
    if ((storage_model != "dense_v1" && storage_model != "sparse_zero_pages_v1") ||
        maximum_active_voxels < 1 || maximum_active_voxels > 1000000) {
        throw std::invalid_argument("invalid structured storage model/budget");
    }
    if (storage_model == "sparse_zero_pages_v1" &&
        (!continuum.base.thin_layer || continuum.angiogenesis.model != "disabled")) {
        throw std::invalid_argument(
            "sparse zero pages v1 requires thin layer and static vasculature");
    }
    continuum.validate();
    if ((division_clock_model != "mean_rate_v1" &&
         division_clock_model != "transported_shifted_geometric_v1") ||
        (schema_version < 10 && division_clock_model != "mean_rate_v1") ||
        !(division_work_bin_width > 0.0) || !std::isfinite(division_work_bin_width) ||
        !(division_maximum_work > division_work_bin_width) ||
        !std::isfinite(division_maximum_work) ||
        division_maximum_work / division_work_bin_width > 4096) {
        throw std::invalid_argument("invalid structured division clock model/grid");
    }
    if (schema_version >= 11 && (!(migration.activation_time_bin_width_hours > 0.0) ||
        !std::isfinite(migration.activation_time_bin_width_hours) ||
        !(migration.activation_maximum_hours > migration.activation_time_bin_width_hours) ||
        !std::isfinite(migration.activation_maximum_hours) ||
        migration.activation_maximum_hours / migration.activation_time_bin_width_hours > 4096))
        throw std::invalid_argument("invalid activation duration grid");
    if (migration.activation_rate_model != "phenotype_mean_v1" &&
        migration.activation_rate_model != "beta_rate_distribution_v2")
        throw std::invalid_argument("unsupported activation rate closure");
    if (schema_version >= 12 && (migration.activation_rate_bins < 4 || migration.activation_rate_bins > 64))
        throw std::invalid_argument("activation rate bins must be between 4 and 64");
    if (migration.activation_rate_model == "beta_rate_distribution_v2" &&
        (schema_version < 12 || migration.activation_clock != "beta_duration_distribution_v2" ||
         continuum.base.normal_r_migration_beta.lower_clamp_enabled ||
         !(continuum.base.normal_r_migration_beta.scale > 0.0) ||
         migration.activation_rate_bins < 4 || migration.activation_rate_bins > 64))
        throw std::invalid_argument("rate distribution requires v12, duration distributions and an unclamped positive beta rate");
    if (migration.activation_clock == "beta_duration_distribution_v2") {
        if (schema_version < 11) throw std::invalid_argument("activation duration distribution requires schema v11");
        const auto& base = continuum.base;
        const double minimum_rate = base.initial_growth_rate_model == "fixed"
            ? base.initial_r_growth_rate : base.initial_r_growth_truncated_normal.minimum;
        (void)beta_duration_kernel(base.migration_activation_duration_alpha,
            base.migration_activation_duration_beta, base.division_timing.base_cycle_hours / minimum_rate,
            migration.activation_time_bin_width_hours, migration.activation_maximum_hours);
    }
    if (division_clock_model == "transported_shifted_geometric_v1") {
        const auto& base = continuum.base;
        std::array<double, 2> inherent{base.initial_r_growth_rate, base.initial_K_growth_rate};
        if (base.initial_growth_rate_model != "fixed") {
            const auto& r = base.initial_r_growth_truncated_normal;
            const auto& K = base.initial_K_growth_truncated_normal;
            inherent = {std::clamp(r.mean, r.minimum, r.maximum),
                        std::clamp(K.mean, K.minimum, K.maximum)};
        }
        (void)DivisionRenewal3D(base.division_timing, division_work_bin_width,
                              division_maximum_work, inherent);
    }
    if (migration.model != "abm_activation_clock_discrete_velocity_v1" ||
        migration.activation_density != "abm_anchor_box_v1" ||
        (migration.activation_clock != "beta_mean_remaining_cycle_v1" &&
         migration.activation_clock != "beta_duration_distribution_v2") ||
        (migration.direction_transport != "fixed_direction_jump_v1" &&
         migration.direction_transport !=
             "guided_fixed_direction_jump_v2" &&
         migration.direction_transport !=
             "guided_fixed_direction_jump_exchange_v3" &&
         migration.direction_transport !=
             "nutrient_gradient_fixed_direction_jump_exchange_v4" &&
         migration.direction_transport !=
             "nutrient_gradient_feasible_direction_jump_exchange_v5")) {
        throw std::invalid_argument("unsupported structured migration model");
    }
    if (!(migration.maximum_move_probability_per_substep > 0.0) ||
        migration.maximum_move_probability_per_substep > 0.25 ||
        !(migration.minimum_density > 0.0)) {
        throw std::invalid_argument("structured migration numerics are invalid");
    }
    const auto& base = continuum.base;
    if (migration.direction_transport == "nutrient_gradient_feasible_direction_jump_exchange_v5" &&
        base.turn_half_angle_degrees > 45.0)
        throw std::invalid_argument("feasible-direction subset closure requires a turn cone at most 45 degrees");
    if (!base.migration_activation_enabled ||
        base.activated_r_migration_rate_model != "normal_multiplier" ||
        !close(base.activated_r_normal_multiplier,
               continuum.migration.activated_r_mobility_multiplier) ||
        base.direction_set != "fixed_26_v1") {
        throw std::invalid_argument(
            "structured PDE requires the shared ABM normal-multiplier contract");
    }
    if (schema_version == 2 &&
        (continuum.schema_version != 2 ||
         migration.direction_transport !=
             "guided_fixed_direction_jump_v2" ||
         base.direction_guidance_model !=
             "low_density_high_resource_v1" ||
         continuum.nutrient.consumption_model != "per_cell_ratio_v2")) {
        throw std::invalid_argument(
            "structured PDE v2 requires per-cell resource consumption and "
            "low-density/high-resource guidance");
    }
    if (schema_version >= 3 && schema_version < 5 &&
        (continuum.schema_version != 2 ||
         migration.direction_transport !=
             "guided_fixed_direction_jump_exchange_v3" ||
         migration.direction_density_window_edge != 70 ||
         migration.crowding_exchange !=
             "active_r_K_stage1_conservative_v1" ||
         !migration.vessel_exclusion ||
         base.direction_guidance_model !=
             "low_density_high_resource_v1" ||
         continuum.nutrient.consumption_model != "per_cell_ratio_v2" ||
         base.r_to_K_conversion.density_window_edge !=
             base.migration_activation_window_edge ||
         base.r_to_K_conversion.query_block_edge !=
             base.migration_activation_block_edge ||
         base.r_to_K_conversion.density_window_edge != 70)) {
        throw std::invalid_argument(
            "structured PDE v3+ requires a 70-voxel shared density window, "
            "conservative active-r/K exchange, hard vessel exclusion, and "
            "v2 resource guidance");
    }
    if (schema_version >= 4 &&
        (continuum.reaction.model !=
             "abm_work_clock_neighbor_availability_v2" ||
         !close(continuum.reaction.small_daughter_vacancy_exponent,
                base.thin_layer ? 8.0 : 26.0))) {
        throw std::invalid_argument(
            "structured PDE v4 requires the ABM work clock and the exact "
            "small-cell neighbour count");
    }
    const bool v5_boundary =
        continuum.nutrient.boundary_mode ==
            "planar_edges_and_vessels_dirichlet_v1" ||
        continuum.nutrient.boundary_mode ==
            "planar_edges_dirichlet_v1" ||
        continuum.nutrient.boundary_mode == "vessels_dirichlet_v1";
    const bool v6_boundary =
        continuum.nutrient.boundary_mode ==
            "moving_tumor_front_dirichlet_v2" ||
        continuum.nutrient.boundary_mode ==
            "moving_tumor_front_and_vessels_dirichlet_v2";
    if (schema_version >= 5 &&
        ((schema_version == 5 && continuum.schema_version != 3) ||
         (schema_version == 6 && continuum.schema_version != 4) ||
         (schema_version == 7 && continuum.schema_version != 5) ||
         (schema_version >= 8 && continuum.schema_version != 6) ||
         migration.activation_stop != (schema_version >= 7
             ? "cohort_clock_refractory_hysteresis_v3"
             : "clock_expiry_refractory_hysteresis_v2") ||
         (migration.direction_transport !=
             "nutrient_gradient_fixed_direction_jump_exchange_v4" &&
          migration.direction_transport !=
             "nutrient_gradient_feasible_direction_jump_exchange_v5") ||
         (base.direction_guidance_model != "low_density_high_resource_v1" &&
          base.direction_guidance_model != "low_density_high_resource_bounded_v2" &&
          base.direction_guidance_model != "nutrient_gradient_shared_resource_v3" &&
        base.direction_guidance_model != "nutrient_gradient_shared_resource_v4") ||
         migration.direction_nutrient_window_edge != 70 ||
         !(migration.chemotaxis_strength > 0.0) ||
         !(migration.zero_gradient_tolerance > 0.0) ||
         migration.reactivation_cooldown_hours < 0.0 ||
         !(migration.reactivation_density_threshold >= 0.0) ||
         !(migration.reactivation_density_threshold <
             base.migration_activation_threshold) ||
         migration.crowding_exchange !=
             "active_r_K_stage1_conservative_v1" ||
         !migration.vessel_exclusion ||
         continuum.nutrient.model != "transient_shared_resource_v2" ||
         continuum.nutrient.solver !=
             "transient_explicit_dirichlet_sources_v2" ||
         (schema_version == 5 && !v5_boundary) ||
         (schema_version == 6 && !v6_boundary) ||
         (schema_version >= 7 && !v5_boundary && !v6_boundary) ||
         continuum.nutrient.consumption_model != "per_cell_ratio_v2" ||
         !close(continuum.nutrient.r_consumption_rate_per_hour,
                continuum.nutrient.K_consumption_rate_per_hour) ||
         !close(continuum.nutrient.maximum_capacity_multiplier, 1.0))) {
        throw std::invalid_argument(
            "structured PDE v5+ requires the versioned transient resource "
            "boundary, "
            "equal per-cell demand, fixed capacity, pure nutrient-gradient "
            "guidance, and clock-expiry hysteresis");
    }
}

std::uint64_t StructuredPdeConfig3D::dynamics_fingerprint() const {
    std::uint64_t state = mix(
        0x4154434753504445ULL, continuum.dynamics_fingerprint());
    state = mix(state, static_cast<std::uint64_t>(schema_version));
    hash_text(state, migration.model);
    hash_text(state, migration.activation_density);
    hash_text(state, migration.activation_clock);
    if (schema_version >= 12) {
        hash_text(state, migration.activation_rate_model);
        state = mix(state, migration.activation_rate_bins);
    }
    if (schema_version >= 11) {
        state = mix(state, std::bit_cast<std::uint64_t>(migration.activation_time_bin_width_hours));
        state = mix(state, std::bit_cast<std::uint64_t>(migration.activation_maximum_hours));
    }
    hash_text(state, migration.activation_stop);
    hash_text(state, migration.direction_transport);
    state = mix(state, std::bit_cast<std::uint64_t>(
        migration.maximum_move_probability_per_substep));
    state = mix(state,
        std::bit_cast<std::uint64_t>(migration.minimum_density));
    if (schema_version >= 3) {
        state = mix(state, static_cast<std::uint64_t>(schema_version >= 5
            ? migration.direction_nutrient_window_edge
            : migration.direction_density_window_edge));
        hash_text(state, migration.crowding_exchange);
        state = mix(state, migration.vessel_exclusion ? 1U : 0U);
    }
    if (schema_version >= 5) {
        state = mix(state, std::bit_cast<std::uint64_t>(
            migration.chemotaxis_strength));
        state = mix(state, std::bit_cast<std::uint64_t>(
            migration.zero_gradient_tolerance));
        state = mix(state, std::bit_cast<std::uint64_t>(
            migration.reactivation_cooldown_hours));
        state = mix(state, std::bit_cast<std::uint64_t>(
            migration.reactivation_density_threshold));
    }
    if (schema_version >= 9) {
        hash_text(state, storage_model);
        state = mix(state, maximum_active_voxels);
    }
    if (schema_version >= 10) {
        hash_text(state, division_clock_model);
        state = mix(state, std::bit_cast<std::uint64_t>(division_work_bin_width));
        state = mix(state, std::bit_cast<std::uint64_t>(division_maximum_work));
    }
    return state;
}

std::string StructuredPdeConfig3D::to_json() const {
    std::ostringstream stream;
    stream << std::setprecision(17)
           << "{\"schema\":\"atcg3d.structured_pde_config\","
           << "\"version\":" << schema_version
           << ",\"profile\":\"" << escaped(profile) << "\","
           << "\"continuum_config\":" << continuum.to_json() << ','
           << "\"structured_migration\":{"
           << "\"model\":\"" << escaped(migration.model) << "\","
           << "\"activation_density\":\""
           << escaped(migration.activation_density) << "\","
           << "\"activation_clock\":\""
           << escaped(migration.activation_clock) << "\","
           << "\"activation_stop\":\""
           << escaped(migration.activation_stop) << "\","
           << "\"direction_transport\":\""
           << escaped(migration.direction_transport) << "\","
           << "\"direction_density_window_edge\":"
           << migration.direction_density_window_edge << ','
           << "\"direction_nutrient_window_edge\":"
           << migration.direction_nutrient_window_edge << ','
           << "\"chemotaxis_strength\":"
           << migration.chemotaxis_strength << ','
           << "\"zero_gradient_tolerance\":"
           << migration.zero_gradient_tolerance << ','
           << "\"reactivation_cooldown_hours\":"
           << migration.reactivation_cooldown_hours << ','
           << "\"reactivation_density_threshold\":"
           << migration.reactivation_density_threshold << ','
           << "\"crowding_exchange\":\""
           << escaped(migration.crowding_exchange) << "\","
           << "\"vessel_exclusion\":"
           << (migration.vessel_exclusion ? "true" : "false") << ','
           << "\"maximum_move_probability_per_substep\":"
           << migration.maximum_move_probability_per_substep << ','
           << "\"minimum_density\":" << migration.minimum_density
           << "}";
    if (schema_version >= 12) {
        stream << ",\"activation_rate_distribution\":{\"model\":\"" << migration.activation_rate_model
               << "\",\"bins\":" << migration.activation_rate_bins << '}';
    }
    if (schema_version >= 11) {
        stream << ",\"activation_duration_grid\":{\"width_hours\":" << migration.activation_time_bin_width_hours
               << ",\"maximum_hours\":" << migration.activation_maximum_hours << '}';
    }
    if (schema_version >= 9) {
        stream << ",\"storage\":{\"model\":\"" << storage_model
               << "\",\"maximum_active_voxels\":" << maximum_active_voxels << '}';
    }
    if (schema_version >= 10) {
        stream << ",\"division_clock\":{\"model\":\"" << division_clock_model
               << "\",\"work_bin_width\":" << division_work_bin_width
               << ",\"maximum_work\":" << division_maximum_work << '}';
    }
    stream << '}';
    return stream.str();
}

}  // namespace atcg3d::structured_pde
