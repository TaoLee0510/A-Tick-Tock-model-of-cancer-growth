#include "config/model_config.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string_view>

namespace atcg3d {
namespace {

std::string trim(std::string value) {
    const auto is_space = [](unsigned char character) { return std::isspace(character) != 0; };
    value.erase(value.begin(), std::find_if_not(value.begin(), value.end(), is_space));
    value.erase(std::find_if_not(value.rbegin(), value.rend(), is_space).base(), value.end());
    return value;
}

template <class Integer>
Integer parse_integer(const std::string& value, const std::string& key) {
    Integer parsed{};
    const char* begin = value.data();
    const char* end = begin + value.size();
    const auto [pointer, error] = std::from_chars(begin, end, parsed);
    if (error != std::errc{} || pointer != end) {
        throw std::invalid_argument("invalid integer for " + key + ": " + value);
    }
    return parsed;
}

double parse_double(const std::string& value, const std::string& key) {
    std::size_t consumed = 0;
    double parsed = 0.0;
    try {
        parsed = std::stod(value, &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument("invalid floating point value for " + key + ": " + value);
    }
    if (consumed != value.size() || !std::isfinite(parsed)) {
        throw std::invalid_argument("invalid finite value for " + key + ": " + value);
    }
    return parsed;
}

bool parse_bool(const std::string& value, const std::string& key) {
    if (value == "true" || value == "1") {
        return true;
    }
    if (value == "false" || value == "0") {
        return false;
    }
    throw std::invalid_argument("invalid boolean for " + key + ": " + value);
}

std::string unquote(const std::string& value) {
    if (value.size() >= 2 && ((value.front() == '"' && value.back() == '"') ||
                              (value.front() == '\'' && value.back() == '\''))) {
        return value.substr(1, value.size() - 2);
    }
    return value;
}

void set_value(Model3DConfig& config, const std::string& key, const std::string& raw_value) {
    const std::string value = trim(raw_value);
    if (key == "schema.version") config.schema_version = parse_integer<std::uint32_t>(value, key);
    else if (key == "profile") config.profile = unquote(value);
    else if (key == "rng.seed") config.seed = parse_integer<std::uint64_t>(value, key);
    else if (key == "space.domain_policy") {
        config.domain_policy = unquote(value);
        if (config.domain_policy == "bounded") config.bounded_domain = true;
        else if (config.domain_policy == "expandable_sparse") config.bounded_domain = false;
    }
    else if (key == "space.chunk_edge") config.chunk_edge = parse_integer<int>(value, key);
    else if (key == "space.min_x") config.domain_min.x = parse_integer<std::int32_t>(value, key);
    else if (key == "space.min_y") config.domain_min.y = parse_integer<std::int32_t>(value, key);
    else if (key == "space.min_z") config.domain_min.z = parse_integer<std::int32_t>(value, key);
    else if (key == "space.max_x") config.domain_max.x = parse_integer<std::int32_t>(value, key);
    else if (key == "space.max_y") config.domain_max.y = parse_integer<std::int32_t>(value, key);
    else if (key == "space.max_z") config.domain_max.z = parse_integer<std::int32_t>(value, key);
    else if (key == "space.thin_layer") config.thin_layer = parse_bool(value, key);
    else if (key == "direction.set") config.direction_set = unquote(value);
    else if (key == "direction.continue_probability") config.continue_probability = parse_double(value, key);
    else if (key == "direction.turn_half_angle_deg") config.turn_half_angle_degrees = parse_double(value, key);
    else if (key == "direction.density_radius") config.direction_density_radius = parse_integer<int>(value, key);
    else if (key == "direction.density_half_angle_deg") config.direction_density_half_angle_degrees = parse_double(value, key);
    else if (key == "direction.density_threshold") config.direction_density_threshold = parse_double(value, key);
    else if (key == "direction.persistence_uses_density") config.persistence_uses_density = parse_bool(value, key);
    else if (key == "direction.distance_weight_exponent") config.distance_weight_exponent = parse_double(value, key);
    else if (key == "migration.activation_enabled") config.migration_activation_enabled = parse_bool(value, key);
    else if (key == "migration.activation_window_edge") config.migration_activation_window_edge = parse_integer<int>(value, key);
    else if (key == "migration.activation_block_edge") config.migration_activation_block_edge = parse_integer<int>(value, key);
    else if (key == "migration.activation_threshold") config.migration_activation_threshold = parse_double(value, key);
    else if (key == "density.backend") config.density_backend = unquote(value);
    else if (key == "density.block_edge") config.density_block_edge = parse_integer<int>(value, key);
    else if (key == "density.growth_window_edge") config.growth_density_window_edge = parse_integer<int>(value, key);
    else if (key == "density.capacity_scale_2d_to_3d") config.carrying_capacity_scale_2d_to_3d = parse_double(value, key);
    else if (key == "biology.r_limit") config.r_limit = parse_double(value, key);
    else if (key == "biology.K_limit") config.K_limit = parse_double(value, key);
    else if (key == "biology.carrying_capacity_r") config.carrying_capacity_r = parse_double(value, key);
    else if (key == "biology.carrying_capacity_K") config.carrying_capacity_K = parse_double(value, key);
    else if (key == "biology.alpha") config.alpha = parse_double(value, key);
    else if (key == "biology.beta") config.beta = parse_double(value, key);
    else if (key == "stage.large_footprint_edge") config.large_footprint_edge = parse_integer<int>(value, key);
    else if (key == "stage.small_footprint_voxels") config.small_footprint_voxels = parse_integer<int>(value, key);
    else if (key == "division.shell_radius") config.division_shell_radius = parse_integer<int>(value, key);
    else if (key == "division.allow_shape_reduction") config.allow_shape_reduction = parse_bool(value, key);
    else if (key == "stage.ultrasmall_enabled") config.ultrasmall_enabled = parse_bool(value, key);
    else if (key == "initial.r_cells") config.initial_r_cells = parse_integer<std::uint64_t>(value, key);
    else if (key == "initial.K_cells") config.initial_K_cells = parse_integer<std::uint64_t>(value, key);
    else if (key == "initial.radius") config.initial_radius = parse_integer<int>(value, key);
    else if (key == "initial.shell_thickness") config.initial_shell_thickness = parse_integer<int>(value, key);
    else if (key == "initial.r_growth_rate") config.initial_r_growth_rate = parse_double(value, key);
    else if (key == "initial.K_growth_rate") config.initial_K_growth_rate = parse_double(value, key);
    else if (key == "initial.r_migration_rate") config.initial_r_migration_rate = parse_double(value, key);
    else if (key == "initial.K_migration_rate") config.initial_K_migration_rate = parse_double(value, key);
    else if (key == "biology.r_death_delay_hours") config.r_death_delay_hours = parse_double(value, key);
    else if (key == "biology.K_death_delay_hours") config.K_death_delay_hours = parse_double(value, key);
    else if (key == "simulation.end_time_hours") config.end_time_hours = parse_double(value, key);
    else if (key == "simulation.max_events") config.max_events = parse_integer<std::uint64_t>(value, key);
    else if (key == "simulation.threads") config.threads = parse_integer<int>(value, key);
    else if (key == "scheduler.backend") config.scheduler_backend = unquote(value);
    else if (key == "scheduler.conflict_bucket_hours") config.conflict_bucket_hours = parse_double(value, key);
    else if (key == "output.enabled") config.output_enabled = parse_bool(value, key);
    else if (key == "output.preview_mode") config.preview_mode = unquote(value);
    else if (key == "output.full_format") config.full_format = unquote(value);
    else if (key == "output.checkpoint_format") config.checkpoint_format = unquote(value);
    else if (key == "output.directory") config.output_directory = unquote(value);
    else if (key == "output.preview_every_hours") config.preview_every_hours = parse_double(value, key);
    else if (key == "output.full_every_hours") config.full_every_hours = parse_double(value, key);
    else if (key == "output.checkpoint_every_hours") config.checkpoint_every_hours = parse_double(value, key);
    else if (key == "output.preview_max_cells") config.preview_max_cells = parse_integer<std::uint64_t>(value, key);
    else if (key == "output.preview_seed") config.preview_seed = parse_integer<std::uint64_t>(value, key);
    else throw std::invalid_argument("unknown configuration key: " + key);
}

std::string json_escape(const std::string& value) {
    std::string result;
    result.reserve(value.size());
    for (const char character : value) {
        if (character == '\\' || character == '"') {
            result.push_back('\\');
        }
        result.push_back(character);
    }
    return result;
}

}  // namespace

Model3DConfig Model3DConfig::load(const std::filesystem::path& path) {
    std::ifstream stream(path);
    if (!stream) {
        throw std::runtime_error("unable to open configuration: " + path.string());
    }
    Model3DConfig config;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(stream, line)) {
        ++line_number;
        const std::size_t comment = line.find('#');
        if (comment != std::string::npos) {
            line.erase(comment);
        }
        line = trim(line);
        if (line.empty()) {
            continue;
        }
        const std::size_t equals = line.find('=');
        if (equals == std::string::npos) {
            throw std::invalid_argument("configuration line " + std::to_string(line_number) +
                                        " does not contain '='");
        }
        try {
            set_value(config, trim(line.substr(0, equals)), trim(line.substr(equals + 1)));
        } catch (const std::exception& error) {
            throw std::invalid_argument("configuration line " + std::to_string(line_number) +
                                        ": " + error.what());
        }
    }
    config.validate();
    return config;
}

void Model3DConfig::apply_override(const std::string& assignment) {
    const std::size_t equals = assignment.find('=');
    if (equals == std::string::npos) {
        throw std::invalid_argument("override must be key=value: " + assignment);
    }
    Model3DConfig candidate = *this;
    set_value(candidate, trim(assignment.substr(0, equals)), trim(assignment.substr(equals + 1)));
    candidate.validate();
    *this = std::move(candidate);
}

void Model3DConfig::validate() const {
    const auto probability = [](double value, const char* name) {
        if (!std::isfinite(value) || value < 0.0 || value > 1.0) {
            throw std::invalid_argument(std::string(name) + " must be between 0 and 1");
        }
    };
    if (schema_version != 1) throw std::invalid_argument("unsupported schema.version");
    if (profile.empty()) throw std::invalid_argument("profile must not be empty");
    if (domain_policy != "expandable_sparse" && domain_policy != "bounded") {
        throw std::invalid_argument("space.domain_policy must be expandable_sparse or bounded");
    }
    if ((domain_policy == "bounded") != bounded_domain) {
        throw std::invalid_argument("space.domain_policy contradicts space.bounded");
    }
    if (chunk_edge < 4 || chunk_edge > 128) throw std::invalid_argument("space.chunk_edge must be in [4,128]");
    if (bounded_domain && !(domain_min.x <= domain_max.x && domain_min.y <= domain_max.y && domain_min.z <= domain_max.z)) {
        throw std::invalid_argument("domain minimum must not exceed maximum");
    }
    if (thin_layer && bounded_domain && (domain_min.z > 0 || domain_max.z < 0)) {
        throw std::invalid_argument("thin-layer bounded domain must include z=0");
    }
    if (direction_set != "fixed_26_v1") throw std::invalid_argument("unsupported direction.set");
    probability(continue_probability, "direction.continue_probability");
    probability(direction_density_threshold, "direction.density_threshold");
    if (turn_half_angle_degrees <= 0.0 || turn_half_angle_degrees > 180.0) throw std::invalid_argument("invalid turn angle");
    if (direction_density_half_angle_degrees <= 0.0 || direction_density_half_angle_degrees > 180.0) throw std::invalid_argument("invalid density cone angle");
    if (direction_density_radius <= 0 || direction_density_radius > 128) throw std::invalid_argument("invalid density radius");
    if (distance_weight_exponent < 0.0) throw std::invalid_argument("distance weight exponent must be non-negative");
    if (migration_activation_window_edge <= 0 || migration_activation_window_edge > 4096) {
        throw std::invalid_argument("invalid migration activation window");
    }
    if (migration_activation_block_edge <= 0 ||
        migration_activation_block_edge > migration_activation_window_edge) {
        throw std::invalid_argument("invalid migration activation block edge");
    }
    probability(migration_activation_threshold, "migration.activation_threshold");
    if (density_backend != "block_anchor_v1") throw std::invalid_argument("unsupported density.backend");
    if (density_block_edge <= 0 || density_block_edge > chunk_edge) throw std::invalid_argument("invalid density block edge");
    if (growth_density_window_edge <= 0 || growth_density_window_edge > 256) throw std::invalid_argument("invalid growth density window");
    if (carrying_capacity_scale_2d_to_3d <= 0.0 || carrying_capacity_r <= 0.0 || carrying_capacity_K <= 0.0) throw std::invalid_argument("carrying capacities must be positive");
    if (r_limit < 0.0 || K_limit < 0.0 || alpha < 0.0 || beta < 0.0) {
        throw std::invalid_argument("density limits and interaction coefficients must be non-negative");
    }
    if (large_footprint_edge != 2 || small_footprint_voxels != 1) {
        throw std::invalid_argument("unsupported stage footprint configuration");
    }
    if (division_shell_radius != 2) throw std::invalid_argument("division shell radius must be 2 for legacy_like_v1");
    if (initial_radius < 0 || initial_shell_thickness < 0) throw std::invalid_argument("initial geometry must be non-negative");
    if (initial_r_cells > std::numeric_limits<std::uint64_t>::max() - initial_K_cells ||
        initial_r_cells + initial_K_cells >= static_cast<std::uint64_t>(kEmptySlot)) {
        throw std::invalid_argument("initial cell count exceeds stable-slot capacity");
    }
    if (initial_r_growth_rate <= 0.0 || initial_K_growth_rate <= 0.0) throw std::invalid_argument("growth rates must be positive");
    if (initial_r_migration_rate < 0.0 || initial_K_migration_rate < 0.0) throw std::invalid_argument("migration rates must be non-negative");
    if (r_death_delay_hours < 0.0 || K_death_delay_hours < 0.0) throw std::invalid_argument("death delays must be non-negative");
    if (end_time_hours < 0.0 || max_events == 0 || threads <= 0) throw std::invalid_argument("invalid simulation limits");
    if (scheduler_backend != "event_queue_v1") throw std::invalid_argument("unsupported scheduler.backend");
    if (conflict_bucket_hours < 0.0) throw std::invalid_argument("scheduler conflict bucket must be non-negative");
    if (preview_mode != "stable_uid_hash_v1" || full_format != "vtkhdf_points_v1" ||
        checkpoint_format != "hdf5_v1") {
        throw std::invalid_argument("unsupported output strategy");
    }
    if (output_enabled && output_directory.empty()) throw std::invalid_argument("output directory must not be empty");
    if (preview_every_hours < 0.0 || full_every_hours < 0.0 || checkpoint_every_hours < 0.0) throw std::invalid_argument("output intervals must be non-negative");
    if (preview_max_cells == 0) throw std::invalid_argument("preview_max_cells must be positive");
    const double finite_values[] = {
        continue_probability, turn_half_angle_degrees, direction_density_half_angle_degrees,
        direction_density_threshold, distance_weight_exponent, carrying_capacity_scale_2d_to_3d,
        r_limit, K_limit, carrying_capacity_r, carrying_capacity_K, alpha, beta,
        initial_r_growth_rate, initial_K_growth_rate, initial_r_migration_rate,
        initial_K_migration_rate, r_death_delay_hours, K_death_delay_hours,
        end_time_hours, conflict_bucket_hours, preview_every_hours, full_every_hours,
        checkpoint_every_hours, migration_activation_threshold,
    };
    for (const double value : finite_values) {
        if (!std::isfinite(value)) throw std::invalid_argument("configuration contains NaN or infinity");
    }
}

std::string Model3DConfig::to_json() const {
    std::ostringstream out;
    out << std::setprecision(17);
    out << "{\n"
        << "  \"schema_version\": " << schema_version << ",\n"
        << "  \"profile\": \"" << json_escape(profile) << "\",\n"
        << "  \"seed\": " << seed << ",\n"
        << "  \"domain_policy\": \"" << json_escape(domain_policy) << "\",\n"
        << "  \"bounded_domain\": " << (bounded_domain ? "true" : "false") << ",\n"
        << "  \"domain_min\": [" << domain_min.x << ',' << domain_min.y << ',' << domain_min.z << "],\n"
        << "  \"domain_max\": [" << domain_max.x << ',' << domain_max.y << ',' << domain_max.z << "],\n"
        << "  \"chunk_edge\": " << chunk_edge << ",\n"
        << "  \"thin_layer\": " << (thin_layer ? "true" : "false") << ",\n"
        << "  \"direction_set\": \"" << json_escape(direction_set) << "\",\n"
        << "  \"continue_probability\": " << continue_probability << ",\n"
        << "  \"turn_half_angle_degrees\": " << turn_half_angle_degrees << ",\n"
        << "  \"direction_density_radius\": " << direction_density_radius << ",\n"
        << "  \"direction_density_half_angle_degrees\": " << direction_density_half_angle_degrees << ",\n"
        << "  \"direction_density_threshold\": " << direction_density_threshold << ",\n"
        << "  \"persistence_uses_density\": " << (persistence_uses_density ? "true" : "false") << ",\n"
        << "  \"distance_weight_exponent\": " << distance_weight_exponent << ",\n"
        << "  \"migration_activation_enabled\": " << (migration_activation_enabled ? "true" : "false") << ",\n"
        << "  \"migration_activation_window_edge\": " << migration_activation_window_edge << ",\n"
        << "  \"migration_activation_block_edge\": " << migration_activation_block_edge << ",\n"
        << "  \"migration_activation_threshold\": " << migration_activation_threshold << ",\n"
        << "  \"density_backend\": \"" << json_escape(density_backend) << "\",\n"
        << "  \"density_block_edge\": " << density_block_edge << ",\n"
        << "  \"growth_density_window_edge\": " << growth_density_window_edge << ",\n"
        << "  \"carrying_capacity_scale_2d_to_3d\": " << carrying_capacity_scale_2d_to_3d << ",\n"
        << "  \"r_limit\": " << r_limit << ",\n"
        << "  \"K_limit\": " << K_limit << ",\n"
        << "  \"carrying_capacity_r\": " << carrying_capacity_r << ",\n"
        << "  \"carrying_capacity_K\": " << carrying_capacity_K << ",\n"
        << "  \"alpha\": " << alpha << ",\n"
        << "  \"beta\": " << beta << ",\n"
        << "  \"large_footprint_edge\": " << large_footprint_edge << ",\n"
        << "  \"small_footprint_voxels\": " << small_footprint_voxels << ",\n"
        << "  \"division_shell_radius\": " << division_shell_radius << ",\n"
        << "  \"allow_shape_reduction\": " << (allow_shape_reduction ? "true" : "false") << ",\n"
        << "  \"ultrasmall_enabled\": " << (ultrasmall_enabled ? "true" : "false") << ",\n"
        << "  \"initial_r_cells\": " << initial_r_cells << ",\n"
        << "  \"initial_K_cells\": " << initial_K_cells << ",\n"
        << "  \"initial_radius\": " << initial_radius << ",\n"
        << "  \"initial_shell_thickness\": " << initial_shell_thickness << ",\n"
        << "  \"initial_r_growth_rate\": " << initial_r_growth_rate << ",\n"
        << "  \"initial_K_growth_rate\": " << initial_K_growth_rate << ",\n"
        << "  \"initial_r_migration_rate\": " << initial_r_migration_rate << ",\n"
        << "  \"initial_K_migration_rate\": " << initial_K_migration_rate << ",\n"
        << "  \"r_death_delay_hours\": " << r_death_delay_hours << ",\n"
        << "  \"K_death_delay_hours\": " << K_death_delay_hours << ",\n"
        << "  \"end_time_hours\": " << end_time_hours << ",\n"
        << "  \"max_events\": " << max_events << ",\n"
        << "  \"threads\": " << threads << ",\n"
        << "  \"scheduler_backend\": \"" << json_escape(scheduler_backend) << "\",\n"
        << "  \"conflict_bucket_hours\": " << conflict_bucket_hours << ",\n"
        << "  \"output_enabled\": " << (output_enabled ? "true" : "false") << ",\n"
        << "  \"preview_mode\": \"" << json_escape(preview_mode) << "\",\n"
        << "  \"full_format\": \"" << json_escape(full_format) << "\",\n"
        << "  \"checkpoint_format\": \"" << json_escape(checkpoint_format) << "\",\n"
        << "  \"output_directory\": \"" << json_escape(output_directory.string()) << "\",\n"
        << "  \"preview_every_hours\": " << preview_every_hours << ",\n"
        << "  \"full_every_hours\": " << full_every_hours << ",\n"
        << "  \"checkpoint_every_hours\": " << checkpoint_every_hours << ",\n"
        << "  \"preview_max_cells\": " << preview_max_cells << ",\n"
        << "  \"preview_seed\": " << preview_seed << "\n"
        << "}\n";
    return out.str();
}

std::string Model3DConfig::dynamics_json() const {
    Model3DConfig normalized = *this;
    normalized.end_time_hours = 0.0;
    normalized.max_events = 1;
    normalized.threads = 1;
    normalized.output_enabled = false;
    normalized.output_directory.clear();
    normalized.preview_every_hours = 0.0;
    normalized.full_every_hours = 0.0;
    normalized.checkpoint_every_hours = 0.0;
    normalized.preview_max_cells = 1;
    normalized.preview_seed = 0;
    return normalized.to_json();
}

}  // namespace atcg3d
