#include "io/run_manifest.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unordered_set>

#include <yaml-cpp/yaml.h>

namespace atcg3d {
namespace {

std::string escape_json(const std::string& value) {
    std::string result;
    result.reserve(value.size());
    for (const unsigned char character : value) {
        switch (character) {
            case '\\': result += "\\\\"; break;
            case '"': result += "\\\""; break;
            case '\b': result += "\\b"; break;
            case '\f': result += "\\f"; break;
            case '\n': result += "\\n"; break;
            case '\r': result += "\\r"; break;
            case '\t': result += "\\t"; break;
            default:
                if (character < 0x20) {
                    std::ostringstream escaped;
                    escaped << "\\u" << std::hex << std::setw(4)
                            << std::setfill('0') << static_cast<unsigned int>(character);
                    result += escaped.str();
                } else {
                    result.push_back(static_cast<char>(character));
                }
                break;
        }
    }
    return result;
}

[[noreturn]] void malformed(const std::filesystem::path& path,
                            const std::string& detail) {
    throw std::runtime_error("invalid resume metadata " + path.string() + ": " + detail);
}

void require_map(const YAML::Node& node,
                 const std::filesystem::path& path,
                 std::string_view context) {
    if (!node || !node.IsMap()) malformed(path, std::string(context) + " must be an object");
}

void require_sequence(const YAML::Node& node,
                      const std::filesystem::path& path,
                      std::string_view context) {
    if (!node || !node.IsSequence()) {
        malformed(path, std::string(context) + " must be an array");
    }
}

void require_string_sequence(
    const YAML::Node& node,
    const std::filesystem::path& path,
    std::string_view context,
    std::initializer_list<std::string_view> expected) {
    require_sequence(node, path, context);
    if (node.size() != expected.size()) {
        malformed(path, std::string(context) + " has the wrong length");
    }
    std::size_t index = 0;
    for (const std::string_view value : expected) {
        if (!node[index].IsScalar() || node[index].Scalar() != value) {
            malformed(path, std::string(context) +
                                " does not match the schema-v3 array catalog");
        }
        ++index;
    }
}

void require_exact_keys(const YAML::Node& node,
                        const std::filesystem::path& path,
                        std::string_view context,
                        std::initializer_list<std::string_view> expected) {
    require_map(node, path, context);
    std::unordered_set<std::string> allowed;
    for (const std::string_view key : expected) allowed.emplace(key);
    std::unordered_set<std::string> seen;
    for (const auto& item : node) {
        if (!item.first.IsScalar()) malformed(path, std::string(context) + " has a non-string key");
        const std::string key = item.first.Scalar();
        if (!allowed.contains(key)) {
            malformed(path, std::string(context) + " has unknown key " + key);
        }
        if (!seen.insert(key).second) {
            malformed(path, std::string(context) + " has duplicate key " + key);
        }
    }
    for (const std::string_view key : expected) {
        if (!seen.contains(std::string(key))) {
            malformed(path, std::string(context) + " is missing key " + std::string(key));
        }
    }
}

template <class T>
T scalar_as(const YAML::Node& node,
            const std::filesystem::path& path,
            std::string_view context) {
    if (!node || !node.IsScalar()) {
        malformed(path, std::string(context) + " must be a scalar");
    }
    try {
        return node.as<T>();
    } catch (const YAML::Exception&) {
        malformed(path, std::string(context) + " has an invalid value");
    }
}

std::string expected_frame_name(std::string_view relative_directory,
                                std::size_t index) {
    std::ostringstream value;
    value << relative_directory << "/frame_" << std::setw(8) << std::setfill('0')
          << index << ".vtkhdf";
    return value.str();
}

std::vector<SeriesEntry3D> read_series_strict(
    const std::filesystem::path& run_directory,
    const std::filesystem::path& catalog_path,
    std::string_view relative_directory) {
    YAML::Node root;
    try {
        root = YAML::LoadFile(catalog_path.string());
    } catch (const YAML::Exception& error) {
        malformed(catalog_path, error.what());
    }
    require_exact_keys(root, catalog_path, "series",
                       {"file-series-version", "files"});
    if (scalar_as<std::string>(root["file-series-version"], catalog_path,
                               "file-series-version") != "1.0") {
        malformed(catalog_path, "unsupported file-series-version");
    }
    const YAML::Node files = root["files"];
    require_sequence(files, catalog_path, "files");

    std::vector<SeriesEntry3D> result;
    result.reserve(files.size());
    double previous_time = -1.0;
    for (std::size_t index = 0; index < files.size(); ++index) {
        const YAML::Node entry = files[index];
        require_exact_keys(entry, catalog_path, "series entry", {"name", "time"});
        const std::string name = scalar_as<std::string>(entry["name"], catalog_path,
                                                        "series entry name");
        const double time = scalar_as<double>(entry["time"], catalog_path,
                                              "series entry time");
        if (!std::isfinite(time) || time < 0.0) {
            malformed(catalog_path, "series time must be finite and nonnegative");
        }
        if (index != 0 && !(time > previous_time)) {
            malformed(catalog_path, "series times must be strictly increasing");
        }
        const std::string expected = expected_frame_name(relative_directory, index);
        if (name != expected || name.find(".tmp") != std::string::npos) {
            malformed(catalog_path, "non-canonical or temporary frame path " + name);
        }
        const std::filesystem::path frame_path = run_directory / name;
        if (!std::filesystem::is_regular_file(frame_path)) {
            malformed(catalog_path, "referenced frame does not exist: " + name);
        }
        result.push_back({name, time});
        previous_time = time;
    }
    return result;
}

bool same_time(double lhs, double rhs) {
    return std::abs(lhs - rhs) <=
           1e-10 * std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

YAML::Node load_manifest(const std::filesystem::path& path) {
    YAML::Node root;
    try {
        root = YAML::LoadFile(path.string());
    } catch (const YAML::Exception& error) {
        malformed(path, error.what());
    }
    require_exact_keys(
        root, path, "run manifest",
        {"schema", "schema_version", "dimension", "coordinate_system",
         "coordinate_units", "time_units", "preview_series", "full_series",
         "vessel_series", "preview_frames", "full_frames", "vessel_frames",
         "checkpoint_enabled", "cell_topology", "cell_point_arrays",
         "cell_field_arrays", "vessel_topology", "vessel_point_arrays",
         "dynamics_config_json", "effective_config"});
    return root;
}

template <class Writer>
void write_atomic(const std::filesystem::path& path, Writer writer) {
    std::filesystem::create_directories(path.parent_path());
    const std::filesystem::path temporary = path.string() + ".tmp";
    std::filesystem::remove(temporary);
    {
        std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
        if (!stream) throw std::runtime_error("unable to create " + temporary.string());
        writer(stream);
        stream.flush();
        if (!stream) throw std::runtime_error("unable to finish " + temporary.string());
    }
    std::filesystem::rename(temporary, path);
}

}  // namespace

void write_series_atomic(const std::filesystem::path& path,
                         const std::vector<SeriesEntry3D>& entries) {
    double previous_time = -1.0;
    for (std::size_t index = 0; index < entries.size(); ++index) {
        if (!std::isfinite(entries[index].time_hours) || entries[index].time_hours < 0.0 ||
            (index != 0 && !(entries[index].time_hours > previous_time))) {
            throw std::invalid_argument("series times must be finite, nonnegative, and strictly increasing");
        }
        if (entries[index].name.empty() ||
            entries[index].name.find(".tmp") != std::string::npos) {
            throw std::invalid_argument("series entries may not reference empty or temporary paths");
        }
        previous_time = entries[index].time_hours;
    }
    write_atomic(path, [&entries](std::ostream& out) {
        out << std::setprecision(17);
        out << "{\n  \"file-series-version\": \"1.0\",\n  \"files\": [\n";
        for (std::size_t index = 0; index < entries.size(); ++index) {
            out << "    {\"name\": \"" << escape_json(entries[index].name)
                << "\", \"time\": " << entries[index].time_hours << '}';
            if (index + 1 != entries.size()) out << ',';
            out << '\n';
        }
        out << "  ]\n}\n";
    });
}

void write_run_manifest_atomic(const std::filesystem::path& run_directory,
                               const Model3DConfig& config,
                               const std::vector<SeriesEntry3D>& preview,
                               const std::vector<SeriesEntry3D>& full,
                               const std::vector<SeriesEntry3D>& vessels,
                               bool checkpoint_enabled) {
    write_atomic(run_directory / "run.json", [&](std::ostream& out) {
        out << "{\n"
            << "  \"schema\": \"atcg.viz\",\n"
            << "  \"schema_version\": 3,\n"
            << "  \"dimension\": 3,\n"
            << "  \"coordinate_system\": \"cartesian\",\n"
            << "  \"coordinate_units\": \"lattice_voxel\",\n"
            << "  \"time_units\": \"hours\",\n"
            << "  \"preview_series\": \"preview.vtkhdf.series\",\n"
            << "  \"full_series\": \"full.vtkhdf.series\",\n"
            << "  \"vessel_series\": \"vessels.vtkhdf.series\",\n"
            << "  \"preview_frames\": " << preview.size() << ",\n"
            << "  \"full_frames\": " << full.size() << ",\n"
            << "  \"vessel_frames\": " << vessels.size() << ",\n"
            << "  \"checkpoint_enabled\": " << (checkpoint_enabled ? "true" : "false") << ",\n"
            << "  \"cell_topology\": \"vtkPolyData.points_only\",\n"
            << "  \"cell_point_arrays\": [\"cell_id\",\"lesion_id\",\"clone_id\",\"cell_type\",\"stage\",\"viability\",\"display_radius\"],\n"
            << "  \"cell_field_arrays\": [\"total_cell_count\"],\n"
            << "  \"vessel_topology\": \"vtkPolyData.lines\",\n"
            << "  \"vessel_point_arrays\": [\"node_id\",\"vessel_id\",\"source_lesion_id\",\"branch_role\",\"perfused\",\"diameter_voxels\",\"radius_voxels\"],\n"
            << "  \"dynamics_config_json\": \""
            << escape_json(config.dynamics_json()) << "\",\n"
            << "  \"effective_config\": " << config.to_json()
            << "}\n";
    });
}

ExistingRunOutput3D load_existing_run_output(
    const std::filesystem::path& run_directory,
    const Model3DConfig& resume_config) {
    const std::filesystem::path manifest_path = run_directory / "run.json";
    if (!std::filesystem::is_regular_file(manifest_path)) {
        throw std::runtime_error("resume output directory has no run.json: " +
                                 run_directory.string());
    }
    const YAML::Node manifest = load_manifest(manifest_path);
    if (scalar_as<std::string>(manifest["schema"], manifest_path, "schema") !=
            "atcg.viz" ||
        scalar_as<unsigned int>(manifest["schema_version"], manifest_path,
                                "schema_version") != 3U ||
        scalar_as<unsigned int>(manifest["dimension"], manifest_path, "dimension") != 3U) {
        malformed(manifest_path, "unsupported visualization schema");
    }
    if (scalar_as<std::string>(manifest["preview_series"], manifest_path,
                               "preview_series") != "preview.vtkhdf.series" ||
        scalar_as<std::string>(manifest["full_series"], manifest_path,
                               "full_series") != "full.vtkhdf.series" ||
        scalar_as<std::string>(manifest["vessel_series"], manifest_path,
                               "vessel_series") != "vessels.vtkhdf.series") {
        malformed(manifest_path, "non-canonical series path");
    }
    if (scalar_as<std::string>(manifest["coordinate_system"], manifest_path,
                               "coordinate_system") != "cartesian" ||
        scalar_as<std::string>(manifest["coordinate_units"], manifest_path,
                               "coordinate_units") != "lattice_voxel" ||
        scalar_as<std::string>(manifest["time_units"], manifest_path,
                               "time_units") != "hours" ||
        scalar_as<std::string>(manifest["cell_topology"], manifest_path,
                               "cell_topology") != "vtkPolyData.points_only" ||
        scalar_as<std::string>(manifest["vessel_topology"], manifest_path,
                               "vessel_topology") != "vtkPolyData.lines") {
        malformed(manifest_path, "coordinate units or topology do not match schema v3");
    }
    require_string_sequence(
        manifest["cell_point_arrays"], manifest_path, "cell_point_arrays",
        {"cell_id", "lesion_id", "clone_id", "cell_type", "stage",
         "viability", "display_radius"});
    require_string_sequence(manifest["cell_field_arrays"], manifest_path,
                            "cell_field_arrays", {"total_cell_count"});
    require_string_sequence(
        manifest["vessel_point_arrays"], manifest_path, "vessel_point_arrays",
        {"node_id", "vessel_id", "source_lesion_id", "branch_role",
         "perfused", "diameter_voxels", "radius_voxels"});
    (void)scalar_as<bool>(manifest["checkpoint_enabled"], manifest_path,
                          "checkpoint_enabled");
    require_map(manifest["effective_config"], manifest_path,
                "effective_config");
    if (scalar_as<std::string>(manifest["dynamics_config_json"], manifest_path,
                               "dynamics_config_json") != resume_config.dynamics_json()) {
        malformed(manifest_path, "dynamics configuration does not match resume checkpoint");
    }

    ExistingRunOutput3D result;
    result.preview = read_series_strict(run_directory,
                                        run_directory / "preview.vtkhdf.series",
                                        "viz/preview");
    result.full = read_series_strict(run_directory,
                                     run_directory / "full.vtkhdf.series",
                                     "viz/full");
    result.vessels = read_series_strict(run_directory,
                                        run_directory / "vessels.vtkhdf.series",
                                        "viz/vessels");
    const std::size_t manifest_preview = scalar_as<std::size_t>(
        manifest["preview_frames"], manifest_path, "preview_frames");
    const std::size_t manifest_full = scalar_as<std::size_t>(
        manifest["full_frames"], manifest_path, "full_frames");
    const std::size_t manifest_vessels = scalar_as<std::size_t>(
        manifest["vessel_frames"], manifest_path, "vessel_frames");
    if (manifest_preview != result.preview.size() || manifest_full != result.full.size() ||
        manifest_vessels != result.vessels.size()) {
        malformed(manifest_path, "frame counts do not match series catalogs");
    }
    if (result.vessels.size() != result.preview.size()) {
        malformed(manifest_path, "vessel and preview series lengths differ");
    }
    for (std::size_t index = 0; index < result.preview.size(); ++index) {
        if (!same_time(result.preview[index].time_hours,
                       result.vessels[index].time_hours)) {
            malformed(manifest_path, "vessel and preview times differ");
        }
    }
    for (const SeriesEntry3D& full : result.full) {
        const auto preview = std::lower_bound(
            result.preview.begin(), result.preview.end(), full.time_hours,
            [](const SeriesEntry3D& entry, double time) {
                return entry.time_hours < time && !same_time(entry.time_hours, time);
            });
        if (preview == result.preview.end() ||
            !same_time(preview->time_hours, full.time_hours)) {
            malformed(manifest_path, "full frame has no exact-time preview frame");
        }
    }
    return result;
}

}  // namespace atcg3d
