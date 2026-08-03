#include "config/continuum_config.hpp"

#include <bit>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unordered_set>

#include <yaml-cpp/yaml.h>

namespace atcg3d::continuum {
namespace {

[[noreturn]] void fail(const std::string& path, const std::string& message) {
    throw std::invalid_argument(path + ": " + message);
}

YAML::Node required(const YAML::Node& parent,
                    const char* key,
                    const std::string& path) {
    const YAML::Node node = parent[key];
    if (!node) fail(path + "." + key, "missing required field");
    return node;
}

void mapping(const YAML::Node& node,
             const std::string& path,
             std::initializer_list<std::string_view> allowed_keys) {
    if (!node || !node.IsMap()) fail(path, "must be a mapping");
    std::unordered_set<std::string> allowed;
    for (const auto key : allowed_keys) allowed.emplace(key);
    std::unordered_set<std::string> observed;
    for (const auto& entry : node) {
        if (!entry.first.IsScalar()) fail(path, "contains a non-scalar key");
        const std::string key = entry.first.Scalar();
        if (!observed.insert(key).second) fail(path + "." + key, "duplicate field");
        if (!allowed.contains(key)) fail(path + "." + key, "unknown field");
    }
}

std::string text_value(const YAML::Node& node, const std::string& path) {
    if (!node.IsScalar()) fail(path, "must be a scalar string");
    return node.Scalar();
}

double number(const YAML::Node& node, const std::string& path) {
    if (!node.IsScalar()) fail(path, "must be a scalar number");
    double value{};
    try {
        value = node.as<double>();
    } catch (const std::exception&) {
        fail(path, "must be a number");
    }
    if (!std::isfinite(value)) fail(path, "must be finite");
    return value;
}

int integer(const YAML::Node& node, const std::string& path) {
    if (!node.IsScalar()) fail(path, "must be a scalar integer");
    long long value{};
    try {
        value = node.as<long long>();
    } catch (const std::exception&) {
        fail(path, "must be an integer");
    }
    if (value < std::numeric_limits<int>::min() ||
        value > std::numeric_limits<int>::max()) fail(path, "integer is out of range");
    return static_cast<int>(value);
}

bool boolean(const YAML::Node& node, const std::string& path) {
    if (!node.IsScalar()) fail(path, "must be a scalar boolean");
    try {
        return node.as<bool>();
    } catch (const std::exception&) {
        fail(path, "must be true or false");
    }
}

template <class T, std::size_t N, class Reader>
std::array<T, N> fixed_sequence(const YAML::Node& node,
                                const std::string& path,
                                Reader reader) {
    if (!node.IsSequence() || node.size() != N) {
        fail(path, "must be a sequence of length " + std::to_string(N));
    }
    std::array<T, N> result{};
    for (std::size_t index = 0; index < N; ++index) {
        result[index] = reader(node[index], path + "[" + std::to_string(index) + "]");
    }
    return result;
}

std::filesystem::path resolved_path(const YAML::Node& node,
                                    const std::string& path,
                                    const std::filesystem::path& source,
                                    bool nullable) {
    if (node.IsNull()) {
        if (nullable) return {};
        fail(path, "must not be null");
    }
    std::filesystem::path configured = text_value(node, path);
    if (configured.empty()) {
        if (nullable) return {};
        fail(path, "must not be empty");
    }
    return configured.is_absolute()
        ? configured.lexically_normal()
        : (source.parent_path() / configured).lexically_normal();
}

std::string escaped(std::string_view value) {
    std::string result;
    for (const char character : value) {
        switch (character) {
            case '\\': result += "\\\\"; break;
            case '"': result += "\\\""; break;
            case '\n': result += "\\n"; break;
            case '\r': result += "\\r"; break;
            case '\t': result += "\\t"; break;
            default: result += character; break;
        }
    }
    return result;
}

std::uint64_t mix(std::uint64_t state, std::uint64_t value) noexcept {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    value ^= value >> 31U;
    return state ^ (value + (state << 6U) + (state >> 2U));
}

void hash_text(std::uint64_t& state, std::string_view value) noexcept {
    for (const unsigned char character : value) state = mix(state, character);
}

void hash_double(std::uint64_t& state, double value) noexcept {
    state = mix(state, std::bit_cast<std::uint64_t>(value));
}

}  // namespace

ContinuumModelConfig3D ContinuumModelConfig3D::load(
    const std::filesystem::path& path) {
    YAML::Node root;
    try {
        root = YAML::LoadFile(path.string());
    } catch (const std::exception& error) {
        throw std::invalid_argument(
            "unable to load continuum config " + path.string() + ": " + error.what());
    }
    mapping(root, "$", {"schema", "profile", "base_config", "run",
                          "initialization", "grid", "numerics", "migration",
                          "reaction", "nutrient", "vascular", "output"});

    ContinuumModelConfig3D result;
    result.source_path = std::filesystem::absolute(path).lexically_normal();
    const YAML::Node schema = required(root, "schema", "$");
    mapping(schema, "$.schema", {"name", "version"});
    result.schema_name = text_value(required(schema, "name", "$.schema"),
                                    "$.schema.name");
    result.schema_version = integer(required(schema, "version", "$.schema"),
                                    "$.schema.version");
    if (result.schema_name != "atcg3d.continuum_model_config" ||
        result.schema_version != 1) fail("$.schema", "unsupported schema identity");
    result.profile = text_value(required(root, "profile", "$"), "$.profile");
    if (result.profile.empty()) fail("$.profile", "must not be empty");

    const std::filesystem::path configured_base = text_value(
        required(root, "base_config", "$"), "$.base_config");
    result.base_config_path = configured_base.is_absolute()
        ? configured_base.lexically_normal()
        : (result.source_path.parent_path() / configured_base).lexically_normal();
    if (!configured_base.is_absolute() &&
        !std::filesystem::exists(result.base_config_path)) {
        const auto colocated =
            (result.source_path.parent_path() / configured_base.filename())
                .lexically_normal();
        if (std::filesystem::exists(colocated)) result.base_config_path = colocated;
    }
    result.base = Model3DConfig::load(result.base_config_path);

    const YAML::Node run = required(root, "run", "$");
    mapping(run, "$.run", {"mode", "resume_checkpoint"});
    result.run_mode = text_value(required(run, "mode", "$.run"), "$.run.mode");
    result.resume_checkpoint = resolved_path(
        required(run, "resume_checkpoint", "$.run"), "$.run.resume_checkpoint",
        result.source_path, true);

    const YAML::Node initialization = required(root, "initialization", "$");
    mapping(initialization, "$.initialization", {"mode", "abm_checkpoint"});
    result.initialization_mode = text_value(
        required(initialization, "mode", "$.initialization"),
        "$.initialization.mode");
    result.abm_checkpoint = resolved_path(
        required(initialization, "abm_checkpoint", "$.initialization"),
        "$.initialization.abm_checkpoint", result.source_path, true);

    const YAML::Node grid = required(root, "grid", "$");
    mapping(grid, "$.grid", {"shape", "origin", "spacing_voxels"});
    result.grid.shape = fixed_sequence<int, 3>(
        required(grid, "shape", "$.grid"), "$.grid.shape", integer);
    result.grid.origin = fixed_sequence<double, 3>(
        required(grid, "origin", "$.grid"), "$.grid.origin", number);
    result.grid.spacing_voxels = number(
        required(grid, "spacing_voxels", "$.grid"), "$.grid.spacing_voxels");

    const YAML::Node numerics = required(root, "numerics", "$");
    mapping(numerics, "$.numerics",
            {"start_time_hours", "end_time_hours", "time_step_hours"});
    result.start_time_hours = number(
        required(numerics, "start_time_hours", "$.numerics"),
        "$.numerics.start_time_hours");
    result.end_time_hours = number(
        required(numerics, "end_time_hours", "$.numerics"),
        "$.numerics.end_time_hours");
    result.time_step_hours = number(
        required(numerics, "time_step_hours", "$.numerics"),
        "$.numerics.time_step_hours");

    const YAML::Node migration = required(root, "migration", "$");
    mapping(migration, "$.migration",
            {"mapping", "diffusion_scale", "large_mobility_multiplier",
             "activated_r_mobility_multiplier", "crowding_exponent"});
    result.migration.mapping = text_value(
        required(migration, "mapping", "$.migration"), "$.migration.mapping");
    result.migration.diffusion_scale = number(
        required(migration, "diffusion_scale", "$.migration"),
        "$.migration.diffusion_scale");
    result.migration.large_mobility_multiplier = number(
        required(migration, "large_mobility_multiplier", "$.migration"),
        "$.migration.large_mobility_multiplier");
    result.migration.activated_r_mobility_multiplier = number(
        required(migration, "activated_r_mobility_multiplier", "$.migration"),
        "$.migration.activated_r_mobility_multiplier");
    result.migration.crowding_exponent = number(
        required(migration, "crowding_exponent", "$.migration"),
        "$.migration.crowding_exponent");

    const YAML::Node reaction = required(root, "reaction", "$");
    mapping(reaction, "$.reaction",
            {"model", "large_daughter_vacancy_exponent",
             "small_daughter_vacancy_exponent",
             "failed_r_division_death_fraction", "maximum_occupied_fraction"});
    result.reaction.model = text_value(
        required(reaction, "model", "$.reaction"), "$.reaction.model");
    result.reaction.large_daughter_vacancy_exponent = number(
        required(reaction, "large_daughter_vacancy_exponent", "$.reaction"),
        "$.reaction.large_daughter_vacancy_exponent");
    result.reaction.small_daughter_vacancy_exponent = number(
        required(reaction, "small_daughter_vacancy_exponent", "$.reaction"),
        "$.reaction.small_daughter_vacancy_exponent");
    result.reaction.failed_r_division_death_fraction = number(
        required(reaction, "failed_r_division_death_fraction", "$.reaction"),
        "$.reaction.failed_r_division_death_fraction");
    result.reaction.maximum_occupied_fraction = number(
        required(reaction, "maximum_occupied_fraction", "$.reaction"),
        "$.reaction.maximum_occupied_fraction");

    const YAML::Node nutrient = required(root, "nutrient", "$");
    mapping(nutrient, "$.nutrient",
            {"model", "solver", "diffusion_voxels2_per_hour", "decay_per_hour",
             "vessel_exchange_per_hour", "vessel_value", "consumption", "capacity",
             "refresh_every_hours", "solver_iterations", "relaxation"});
    auto& n = result.nutrient;
    n.model = text_value(required(nutrient, "model", "$.nutrient"),
                         "$.nutrient.model");
    n.solver = text_value(required(nutrient, "solver", "$.nutrient"),
                          "$.nutrient.solver");
    n.diffusion_voxels2_per_hour = number(
        required(nutrient, "diffusion_voxels2_per_hour", "$.nutrient"),
        "$.nutrient.diffusion_voxels2_per_hour");
    n.decay_per_hour = number(required(nutrient, "decay_per_hour", "$.nutrient"),
                              "$.nutrient.decay_per_hour");
    n.vessel_exchange_per_hour = number(
        required(nutrient, "vessel_exchange_per_hour", "$.nutrient"),
        "$.nutrient.vessel_exchange_per_hour");
    n.vessel_value = number(required(nutrient, "vessel_value", "$.nutrient"),
                            "$.nutrient.vessel_value");
    n.refresh_every_hours = number(
        required(nutrient, "refresh_every_hours", "$.nutrient"),
        "$.nutrient.refresh_every_hours");
    n.solver_iterations = integer(
        required(nutrient, "solver_iterations", "$.nutrient"),
        "$.nutrient.solver_iterations");
    n.relaxation = number(required(nutrient, "relaxation", "$.nutrient"),
                          "$.nutrient.relaxation");
    const YAML::Node consumption = required(nutrient, "consumption", "$.nutrient");
    mapping(consumption, "$.nutrient.consumption",
            {"r_per_occupied_voxel_hour", "K_per_occupied_voxel_hour",
             "r_half_saturation", "K_half_saturation"});
    n.r_consumption_per_occupied_voxel_hour = number(
        required(consumption, "r_per_occupied_voxel_hour", "$.nutrient.consumption"),
        "$.nutrient.consumption.r_per_occupied_voxel_hour");
    n.K_consumption_per_occupied_voxel_hour = number(
        required(consumption, "K_per_occupied_voxel_hour", "$.nutrient.consumption"),
        "$.nutrient.consumption.K_per_occupied_voxel_hour");
    n.r_consumption_half_saturation = number(
        required(consumption, "r_half_saturation", "$.nutrient.consumption"),
        "$.nutrient.consumption.r_half_saturation");
    n.K_consumption_half_saturation = number(
        required(consumption, "K_half_saturation", "$.nutrient.consumption"),
        "$.nutrient.consumption.K_half_saturation");
    const YAML::Node capacity = required(nutrient, "capacity", "$.nutrient");
    mapping(capacity, "$.nutrient.capacity",
            {"half_saturation", "maximum_multiplier"});
    n.capacity_half_saturation = number(
        required(capacity, "half_saturation", "$.nutrient.capacity"),
        "$.nutrient.capacity.half_saturation");
    n.maximum_capacity_multiplier = number(
        required(capacity, "maximum_multiplier", "$.nutrient.capacity"),
        "$.nutrient.capacity.maximum_multiplier");

    const YAML::Node vascular = required(root, "vascular", "$");
    mapping(vascular, "$.vascular",
            {"source_mode", "synthetic_axis", "synthetic_center",
             "synthetic_radius_voxels"});
    result.vascular.source_mode = text_value(
        required(vascular, "source_mode", "$.vascular"), "$.vascular.source_mode");
    result.vascular.synthetic_axis = text_value(
        required(vascular, "synthetic_axis", "$.vascular"),
        "$.vascular.synthetic_axis");
    result.vascular.synthetic_center = fixed_sequence<double, 3>(
        required(vascular, "synthetic_center", "$.vascular"),
        "$.vascular.synthetic_center", number);
    result.vascular.synthetic_radius_voxels = number(
        required(vascular, "synthetic_radius_voxels", "$.vascular"),
        "$.vascular.synthetic_radius_voxels");

    const YAML::Node output = required(root, "output", "$");
    mapping(output, "$.output",
            {"enabled", "directory", "metrics_every_hours", "field_every_hours",
             "radial_profile_every_hours", "checkpoint_every_hours"});
    result.output.enabled = boolean(required(output, "enabled", "$.output"),
                                    "$.output.enabled");
    result.output.directory = text_value(
        required(output, "directory", "$.output"), "$.output.directory");
    if (result.output.directory.empty()) fail("$.output.directory", "must not be empty");
    result.output.metrics_every_hours = number(
        required(output, "metrics_every_hours", "$.output"),
        "$.output.metrics_every_hours");
    result.output.field_every_hours = number(
        required(output, "field_every_hours", "$.output"),
        "$.output.field_every_hours");
    result.output.radial_profile_every_hours = number(
        required(output, "radial_profile_every_hours", "$.output"),
        "$.output.radial_profile_every_hours");
    result.output.checkpoint_every_hours = number(
        required(output, "checkpoint_every_hours", "$.output"),
        "$.output.checkpoint_every_hours");

    result.validate();
    return result;
}

void ContinuumModelConfig3D::validate() const {
    if (run_mode != "new" && run_mode != "resume") {
        throw std::invalid_argument("continuum run mode must be new or resume");
    }
    if ((run_mode == "resume") != !resume_checkpoint.empty()) {
        throw std::invalid_argument("continuum resume mode/path are inconsistent");
    }
    if (initialization_mode != "base_model" &&
        initialization_mode != "abm_checkpoint") {
        throw std::invalid_argument("unsupported continuum initialization mode");
    }
    if ((initialization_mode == "abm_checkpoint") != !abm_checkpoint.empty()) {
        throw std::invalid_argument("ABM checkpoint initialization/path are inconsistent");
    }
    std::uint64_t cells = 1;
    for (const int extent : grid.shape) {
        if (extent <= 0 || extent > 4096 ||
            cells > 500000000ULL / static_cast<std::uint64_t>(extent)) {
            throw std::invalid_argument("continuum grid dimensions are invalid or too large");
        }
        cells *= static_cast<std::uint64_t>(extent);
    }
    if (base.thin_layer != (grid.shape[2] == 1)) {
        throw std::invalid_argument(
            "continuum grid z extent must be one exactly for a thin-layer base model");
    }
    if (!(grid.spacing_voxels > 0.0) || start_time_hours < 0.0 ||
        !(end_time_hours > start_time_hours) || !(time_step_hours > 0.0) ||
        time_step_hours > end_time_hours - start_time_hours) {
        throw std::invalid_argument("continuum grid/time parameters are invalid");
    }
    if (migration.mapping != "fixed_26_from_base_means_v1" ||
        !(migration.diffusion_scale > 0.0) ||
        !(migration.large_mobility_multiplier > 0.0) ||
        !(migration.activated_r_mobility_multiplier >= 1.0) ||
        migration.crowding_exponent < 0.0) {
        throw std::invalid_argument("continuum migration parameters are invalid");
    }
    if (reaction.model != "abm_mean_clock_volume_filling_v1" ||
        reaction.large_daughter_vacancy_exponent < 0.0 ||
        reaction.small_daughter_vacancy_exponent < 0.0 ||
        reaction.failed_r_division_death_fraction < 0.0 ||
        reaction.failed_r_division_death_fraction > 1.0 ||
        !(reaction.maximum_occupied_fraction > 0.0) ||
        reaction.maximum_occupied_fraction > 1.0) {
        throw std::invalid_argument("continuum reaction parameters are invalid");
    }
    if (nutrient.model != "effective_resource_surplus_v1" ||
        nutrient.solver != "deterministic_quasi_steady_jacobi_v1" ||
        !(nutrient.diffusion_voxels2_per_hour > 0.0) ||
        nutrient.decay_per_hour < 0.0 ||
        nutrient.vessel_exchange_per_hour < 0.0 ||
        !(nutrient.vessel_value > 0.0) ||
        nutrient.r_consumption_per_occupied_voxel_hour < 0.0 ||
        nutrient.K_consumption_per_occupied_voxel_hour < 0.0 ||
        !(nutrient.r_consumption_half_saturation > 0.0) ||
        !(nutrient.K_consumption_half_saturation > 0.0) ||
        !(nutrient.capacity_half_saturation > 0.0) ||
        !(nutrient.maximum_capacity_multiplier >= 1.0) ||
        !(nutrient.refresh_every_hours > 0.0) ||
        nutrient.solver_iterations <= 0 || nutrient.solver_iterations > 100000 ||
        !(nutrient.relaxation > 0.0) || nutrient.relaxation > 1.0) {
        throw std::invalid_argument("continuum nutrient parameters are invalid");
    }
    if (vascular.source_mode != "abm_perfusion" &&
        vascular.source_mode != "synthetic_central_line" &&
        vascular.source_mode != "abm_plus_synthetic_line") {
        throw std::invalid_argument("unsupported continuum vascular source mode");
    }
    if (vascular.synthetic_axis != "x" && vascular.synthetic_axis != "y" &&
        vascular.synthetic_axis != "z") {
        throw std::invalid_argument("synthetic vascular axis must be x, y, or z");
    }
    if (!(vascular.synthetic_radius_voxels > 0.0) ||
        output.metrics_every_hours < 0.0 || output.field_every_hours < 0.0 ||
        output.radial_profile_every_hours < 0.0 ||
        output.checkpoint_every_hours < 0.0 ||
        (output.enabled && output.directory.empty())) {
        throw std::invalid_argument("continuum vascular/output parameters are invalid");
    }
}

std::uint64_t ContinuumModelConfig3D::dynamics_fingerprint() const {
    std::uint64_t state = 0x4154434743504445ULL;
    hash_text(state, base.dynamics_json());
    hash_text(state, initialization_mode);
    for (const int extent : grid.shape) state = mix(state, extent);
    for (const double value : grid.origin) hash_double(state, value);
    hash_double(state, grid.spacing_voxels);
    hash_double(state, time_step_hours);
    hash_text(state, migration.mapping);
    for (const double value : {migration.diffusion_scale,
             migration.large_mobility_multiplier,
             migration.activated_r_mobility_multiplier,
             migration.crowding_exponent,
             reaction.large_daughter_vacancy_exponent,
             reaction.small_daughter_vacancy_exponent,
             reaction.failed_r_division_death_fraction,
             reaction.maximum_occupied_fraction,
             nutrient.diffusion_voxels2_per_hour, nutrient.decay_per_hour,
             nutrient.vessel_exchange_per_hour, nutrient.vessel_value,
             nutrient.r_consumption_per_occupied_voxel_hour,
             nutrient.K_consumption_per_occupied_voxel_hour,
             nutrient.r_consumption_half_saturation,
             nutrient.K_consumption_half_saturation,
             nutrient.capacity_half_saturation,
             nutrient.maximum_capacity_multiplier,
             nutrient.refresh_every_hours, nutrient.relaxation,
             vascular.synthetic_radius_voxels}) hash_double(state, value);
    state = mix(state, nutrient.solver_iterations);
    hash_text(state, reaction.model);
    hash_text(state, nutrient.model);
    hash_text(state, nutrient.solver);
    hash_text(state, vascular.source_mode);
    hash_text(state, vascular.synthetic_axis);
    for (const double value : vascular.synthetic_center) hash_double(state, value);
    return state;
}

std::string ContinuumModelConfig3D::to_json() const {
    std::ostringstream out;
    out << std::setprecision(17)
        << '{'
        << "\"schema_name\":\"" << escaped(schema_name) << "\","
        << "\"schema_version\":" << schema_version << ','
        << "\"profile\":\"" << escaped(profile) << "\","
        << "\"base_config_path\":\""
        << escaped(base_config_path.generic_string()) << "\","
        << "\"run\":{\"mode\":\"" << escaped(run_mode)
        << "\",\"resume_checkpoint\":\""
        << escaped(resume_checkpoint.generic_string()) << "\"},"
        << "\"initialization\":{\"mode\":\""
        << escaped(initialization_mode) << "\",\"abm_checkpoint\":\""
        << escaped(abm_checkpoint.generic_string()) << "\"},"
        << "\"grid\":{\"shape\":[" << grid.shape[0] << ',' << grid.shape[1]
        << ',' << grid.shape[2] << "],\"origin\":[" << grid.origin[0] << ','
        << grid.origin[1] << ',' << grid.origin[2] << "],\"spacing_voxels\":"
        << grid.spacing_voxels << "},"
        << "\"numerics\":{\"start_time_hours\":" << start_time_hours
        << ",\"end_time_hours\":" << end_time_hours
        << ",\"time_step_hours\":" << time_step_hours << "},"
        << "\"migration\":{\"mapping\":\"" << escaped(migration.mapping)
        << "\",\"diffusion_scale\":" << migration.diffusion_scale
        << ",\"large_mobility_multiplier\":"
        << migration.large_mobility_multiplier
        << ",\"activated_r_mobility_multiplier\":"
        << migration.activated_r_mobility_multiplier
        << ",\"crowding_exponent\":" << migration.crowding_exponent << "},"
        << "\"reaction\":{\"model\":\"" << escaped(reaction.model)
        << "\",\"large_daughter_vacancy_exponent\":"
        << reaction.large_daughter_vacancy_exponent
        << ",\"small_daughter_vacancy_exponent\":"
        << reaction.small_daughter_vacancy_exponent
        << ",\"failed_r_division_death_fraction\":"
        << reaction.failed_r_division_death_fraction
        << ",\"maximum_occupied_fraction\":"
        << reaction.maximum_occupied_fraction << "},"
        << "\"nutrient\":{\"model\":\"" << escaped(nutrient.model)
        << "\",\"solver\":\"" << escaped(nutrient.solver)
        << "\",\"diffusion_voxels2_per_hour\":"
        << nutrient.diffusion_voxels2_per_hour
        << ",\"decay_per_hour\":" << nutrient.decay_per_hour
        << ",\"vessel_exchange_per_hour\":"
        << nutrient.vessel_exchange_per_hour
        << ",\"vessel_value\":" << nutrient.vessel_value
        << ",\"r_consumption_per_occupied_voxel_hour\":"
        << nutrient.r_consumption_per_occupied_voxel_hour
        << ",\"K_consumption_per_occupied_voxel_hour\":"
        << nutrient.K_consumption_per_occupied_voxel_hour
        << ",\"r_consumption_half_saturation\":"
        << nutrient.r_consumption_half_saturation
        << ",\"K_consumption_half_saturation\":"
        << nutrient.K_consumption_half_saturation
        << ",\"capacity_half_saturation\":"
        << nutrient.capacity_half_saturation
        << ",\"maximum_capacity_multiplier\":"
        << nutrient.maximum_capacity_multiplier
        << ",\"refresh_every_hours\":" << nutrient.refresh_every_hours
        << ",\"solver_iterations\":" << nutrient.solver_iterations
        << ",\"relaxation\":" << nutrient.relaxation << "},"
        << "\"vascular\":{\"source_mode\":\""
        << escaped(vascular.source_mode) << "\",\"synthetic_axis\":\""
        << escaped(vascular.synthetic_axis) << "\",\"synthetic_center\":["
        << vascular.synthetic_center[0] << ',' << vascular.synthetic_center[1]
        << ',' << vascular.synthetic_center[2]
        << "],\"synthetic_radius_voxels\":"
        << vascular.synthetic_radius_voxels << "},"
        << "\"output\":{\"enabled\":"
        << (output.enabled ? "true" : "false") << ",\"directory\":\""
        << escaped(output.directory.generic_string())
        << "\",\"metrics_every_hours\":" << output.metrics_every_hours
        << ",\"field_every_hours\":" << output.field_every_hours
        << ",\"radial_profile_every_hours\":"
        << output.radial_profile_every_hours
        << ",\"checkpoint_every_hours\":"
        << output.checkpoint_every_hours << "},"
        << "\"dynamics_fingerprint\":" << dynamics_fingerprint() << ','
        << "\"base\":" << base.to_json()
        << '}';
    return out.str();
}

}  // namespace atcg3d::continuum
