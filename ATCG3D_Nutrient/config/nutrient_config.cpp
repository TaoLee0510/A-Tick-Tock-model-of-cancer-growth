#include "config/nutrient_config.hpp"

#include <bit>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unordered_set>

#include <yaml-cpp/yaml.h>

namespace atcg3d::nutrient {
namespace {

[[noreturn]] void config_error(const std::string& path,
                               const std::string& message) {
    throw std::invalid_argument(path + ": " + message);
}

YAML::Node require(const YAML::Node& parent,
                   const char* key,
                   const std::string& path) {
    const YAML::Node value = parent[key];
    if (!value) config_error(path + "." + key, "missing required field");
    return value;
}

void check_mapping(const YAML::Node& node,
                   const std::string& path,
                   std::initializer_list<std::string_view> keys) {
    if (!node || !node.IsMap()) config_error(path, "must be a mapping");
    std::unordered_set<std::string> allowed;
    for (const std::string_view key : keys) allowed.emplace(key);
    std::unordered_set<std::string> observed;
    for (const auto& entry : node) {
        if (!entry.first.IsScalar()) config_error(path, "contains a non-scalar key");
        const std::string key = entry.first.Scalar();
        if (!observed.insert(key).second) {
            config_error(path + "." + key, "duplicate field");
        }
        if (!allowed.contains(key)) {
            config_error(path + "." + key, "unknown field");
        }
    }
}

std::string scalar_string(const YAML::Node& node, const std::string& path) {
    if (!node.IsScalar()) config_error(path, "must be a scalar string");
    return node.Scalar();
}

double finite_double(const YAML::Node& node, const std::string& path) {
    if (!node.IsScalar()) config_error(path, "must be a scalar number");
    double value{};
    try {
        value = node.as<double>();
    } catch (const std::exception&) {
        config_error(path, "must be a number");
    }
    if (!std::isfinite(value)) config_error(path, "must be finite");
    return value;
}

int strict_int(const YAML::Node& node, const std::string& path) {
    if (!node.IsScalar()) config_error(path, "must be a scalar integer");
    long long value{};
    try {
        value = node.as<long long>();
    } catch (const std::exception&) {
        config_error(path, "must be an integer");
    }
    if (value < std::numeric_limits<int>::min() ||
        value > std::numeric_limits<int>::max()) {
        config_error(path, "integer is out of range");
    }
    return static_cast<int>(value);
}

std::string json_escape(std::string_view value) {
    std::string result;
    result.reserve(value.size());
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

}  // namespace

void NutrientFieldConfig3D::validate() const {
    if (model != "effective_resource_surplus_v1") {
        throw std::invalid_argument("unsupported nutrient field model: " + model);
    }
    if (solver != "deterministic_quasi_steady_jacobi_v1") {
        throw std::invalid_argument("unsupported nutrient solver: " + solver);
    }
    if (block_edge <= 0 || block_edge > 64 || halo_voxels < 1 ||
        halo_voxels > 4096) {
        throw std::invalid_argument("nutrient block/halo geometry is invalid");
    }
    if (!(diffusion_voxels2_per_hour > 0.0) || decay_per_hour < 0.0 ||
        vessel_exchange_per_hour < 0.0 || !(vessel_value > 0.0) ||
        r_consumption_per_voxel_hour < 0.0 ||
        K_consumption_per_voxel_hour < 0.0 ||
        !(r_consumption_half_saturation > 0.0) ||
        !(K_consumption_half_saturation > 0.0) ||
        !(capacity_half_saturation > 0.0) ||
        !(maximum_capacity_multiplier >= 1.0) ||
        !(refresh_every_hours > 0.0) || solver_iterations <= 0 ||
        solver_iterations > 100000 || !(relaxation > 0.0) ||
        relaxation > 1.0 || metrics_every_hours < 0.0 ||
        field_snapshot_every_hours < 0.0) {
        throw std::invalid_argument("nutrient numerical/biological parameters are invalid");
    }
}

std::uint64_t NutrientFieldConfig3D::fingerprint() const noexcept {
    std::uint64_t value = 0x415443474e555431ULL;
    for (const unsigned char character : model) value = mix(value, character);
    for (const unsigned char character : solver) value = mix(value, character);
    value = mix(value, static_cast<std::uint64_t>(block_edge));
    value = mix(value, static_cast<std::uint64_t>(halo_voxels));
    for (const double number : {
             diffusion_voxels2_per_hour, decay_per_hour,
             vessel_exchange_per_hour, vessel_value,
             r_consumption_per_voxel_hour, K_consumption_per_voxel_hour,
             r_consumption_half_saturation, K_consumption_half_saturation,
             capacity_half_saturation, maximum_capacity_multiplier,
             refresh_every_hours, relaxation, metrics_every_hours,
             field_snapshot_every_hours}) {
        value = mix(value, std::bit_cast<std::uint64_t>(number));
    }
    value = mix(value, static_cast<std::uint64_t>(solver_iterations));
    return value;
}

std::string NutrientFieldConfig3D::to_json() const {
    std::ostringstream out;
    out << std::setprecision(17)
        << '{'
        << "\"model\":\"" << json_escape(model) << "\","
        << "\"solver\":\"" << json_escape(solver) << "\","
        << "\"block_edge\":" << block_edge << ','
        << "\"halo_voxels\":" << halo_voxels << ','
        << "\"diffusion_voxels2_per_hour\":" << diffusion_voxels2_per_hour << ','
        << "\"decay_per_hour\":" << decay_per_hour << ','
        << "\"vessel_exchange_per_hour\":" << vessel_exchange_per_hour << ','
        << "\"vessel_value\":" << vessel_value << ','
        << "\"r_consumption_per_voxel_hour\":" << r_consumption_per_voxel_hour << ','
        << "\"K_consumption_per_voxel_hour\":" << K_consumption_per_voxel_hour << ','
        << "\"r_consumption_half_saturation\":" << r_consumption_half_saturation << ','
        << "\"K_consumption_half_saturation\":" << K_consumption_half_saturation << ','
        << "\"capacity_half_saturation\":" << capacity_half_saturation << ','
        << "\"maximum_capacity_multiplier\":" << maximum_capacity_multiplier << ','
        << "\"refresh_every_hours\":" << refresh_every_hours << ','
        << "\"solver_iterations\":" << solver_iterations << ','
        << "\"relaxation\":" << relaxation << ','
        << "\"metrics_every_hours\":" << metrics_every_hours << ','
        << "\"field_snapshot_every_hours\":" << field_snapshot_every_hours
        << '}';
    return out.str();
}

NutrientModelConfig3D NutrientModelConfig3D::load(
    const std::filesystem::path& path) {
    YAML::Node root;
    try {
        root = YAML::LoadFile(path.string());
    } catch (const std::exception& error) {
        throw std::invalid_argument(
            "unable to load nutrient config " + path.string() + ": " + error.what());
    }
    check_mapping(root, "$", {"schema", "profile", "base_config", "nutrient"});
    const YAML::Node schema = require(root, "schema", "$");
    check_mapping(schema, "$.schema", {"name", "version"});

    NutrientModelConfig3D result;
    result.source_path = std::filesystem::absolute(path).lexically_normal();
    result.schema_name = scalar_string(
        require(schema, "name", "$.schema"), "$.schema.name");
    result.schema_version = strict_int(
        require(schema, "version", "$.schema"), "$.schema.version");
    if (result.schema_name != "atcg3d.nutrient_model_config" ||
        result.schema_version != 1) {
        config_error("$.schema", "unsupported nutrient schema identity");
    }
    result.profile = scalar_string(require(root, "profile", "$"), "$.profile");
    if (result.profile.empty()) config_error("$.profile", "must not be empty");
    const std::filesystem::path configured_base = scalar_string(
        require(root, "base_config", "$"), "$.base_config");
    result.base_config_path = configured_base.is_absolute()
        ? configured_base.lexically_normal()
        : (result.source_path.parent_path() / configured_base).lexically_normal();
    if (!configured_base.is_absolute() &&
        !std::filesystem::exists(result.base_config_path)) {
        const auto colocated_base =
            (result.source_path.parent_path() / configured_base.filename())
                .lexically_normal();
        if (std::filesystem::exists(colocated_base)) {
            result.base_config_path = colocated_base;
        }
    }
    result.base = Model3DConfig::load(result.base_config_path);

    const YAML::Node nutrient = require(root, "nutrient", "$");
    check_mapping(
        nutrient, "$.nutrient",
        {"model", "solver", "block_edge", "halo_voxels", "diffusion_voxels2_per_hour",
         "decay_per_hour", "vessel_exchange_per_hour", "vessel_value", "consumption",
         "capacity", "refresh_every_hours", "solver_iterations", "relaxation", "output"});
    auto& field = result.nutrient;
    field.model = scalar_string(require(nutrient, "model", "$.nutrient"),
                                "$.nutrient.model");
    field.solver = scalar_string(require(nutrient, "solver", "$.nutrient"),
                                 "$.nutrient.solver");
    field.block_edge = strict_int(require(nutrient, "block_edge", "$.nutrient"),
                                  "$.nutrient.block_edge");
    field.halo_voxels = strict_int(require(nutrient, "halo_voxels", "$.nutrient"),
                                   "$.nutrient.halo_voxels");
    field.diffusion_voxels2_per_hour = finite_double(
        require(nutrient, "diffusion_voxels2_per_hour", "$.nutrient"),
        "$.nutrient.diffusion_voxels2_per_hour");
    field.decay_per_hour = finite_double(
        require(nutrient, "decay_per_hour", "$.nutrient"),
        "$.nutrient.decay_per_hour");
    field.vessel_exchange_per_hour = finite_double(
        require(nutrient, "vessel_exchange_per_hour", "$.nutrient"),
        "$.nutrient.vessel_exchange_per_hour");
    field.vessel_value = finite_double(
        require(nutrient, "vessel_value", "$.nutrient"),
        "$.nutrient.vessel_value");
    field.refresh_every_hours = finite_double(
        require(nutrient, "refresh_every_hours", "$.nutrient"),
        "$.nutrient.refresh_every_hours");
    field.solver_iterations = strict_int(
        require(nutrient, "solver_iterations", "$.nutrient"),
        "$.nutrient.solver_iterations");
    field.relaxation = finite_double(
        require(nutrient, "relaxation", "$.nutrient"),
        "$.nutrient.relaxation");

    const YAML::Node consumption = require(nutrient, "consumption", "$.nutrient");
    check_mapping(consumption, "$.nutrient.consumption",
                  {"r_per_voxel_hour", "K_per_voxel_hour",
                   "r_half_saturation", "K_half_saturation"});
    field.r_consumption_per_voxel_hour = finite_double(
        require(consumption, "r_per_voxel_hour", "$.nutrient.consumption"),
        "$.nutrient.consumption.r_per_voxel_hour");
    field.K_consumption_per_voxel_hour = finite_double(
        require(consumption, "K_per_voxel_hour", "$.nutrient.consumption"),
        "$.nutrient.consumption.K_per_voxel_hour");
    field.r_consumption_half_saturation = finite_double(
        require(consumption, "r_half_saturation", "$.nutrient.consumption"),
        "$.nutrient.consumption.r_half_saturation");
    field.K_consumption_half_saturation = finite_double(
        require(consumption, "K_half_saturation", "$.nutrient.consumption"),
        "$.nutrient.consumption.K_half_saturation");

    const YAML::Node capacity = require(nutrient, "capacity", "$.nutrient");
    check_mapping(capacity, "$.nutrient.capacity",
                  {"half_saturation", "maximum_multiplier"});
    field.capacity_half_saturation = finite_double(
        require(capacity, "half_saturation", "$.nutrient.capacity"),
        "$.nutrient.capacity.half_saturation");
    field.maximum_capacity_multiplier = finite_double(
        require(capacity, "maximum_multiplier", "$.nutrient.capacity"),
        "$.nutrient.capacity.maximum_multiplier");

    const YAML::Node output = require(nutrient, "output", "$.nutrient");
    check_mapping(output, "$.nutrient.output",
                  {"metrics_every_hours", "field_snapshot_every_hours"});
    field.metrics_every_hours = finite_double(
        require(output, "metrics_every_hours", "$.nutrient.output"),
        "$.nutrient.output.metrics_every_hours");
    field.field_snapshot_every_hours = finite_double(
        require(output, "field_snapshot_every_hours", "$.nutrient.output"),
        "$.nutrient.output.field_snapshot_every_hours");
    field.validate();
    return result;
}

std::string NutrientModelConfig3D::to_json() const {
    std::ostringstream out;
    out << '{'
        << "\"schema_name\":\"" << json_escape(schema_name) << "\","
        << "\"schema_version\":" << schema_version << ','
        << "\"profile\":\"" << json_escape(profile) << "\","
        << "\"base_config_path\":\""
        << json_escape(base_config_path.generic_string()) << "\","
        << "\"base\":" << base.to_json() << ','
        << "\"nutrient\":" << nutrient.to_json()
        << '}';
    return out.str();
}

}  // namespace atcg3d::nutrient
