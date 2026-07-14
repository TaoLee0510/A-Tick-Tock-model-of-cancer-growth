#include "io/run_manifest.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace atcg3d {
namespace {

std::string escape_json(const std::string& value) {
    std::string result;
    for (const char character : value) {
        if (character == '\\' || character == '"') result.push_back('\\');
        result.push_back(character);
    }
    return result;
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
                               bool checkpoint_enabled) {
    write_atomic(run_directory / "run.json", [&](std::ostream& out) {
        out << "{\n"
            << "  \"schema\": \"atcg.viz\",\n"
            << "  \"schema_version\": 1,\n"
            << "  \"dimension\": 3,\n"
            << "  \"coordinate_system\": \"cartesian\",\n"
            << "  \"coordinate_units\": \"lattice_voxel\",\n"
            << "  \"time_units\": \"hours\",\n"
            << "  \"preview_series\": \"preview.vtkhdf.series\",\n"
            << "  \"full_series\": \"full.vtkhdf.series\",\n"
            << "  \"preview_frames\": " << preview.size() << ",\n"
            << "  \"full_frames\": " << full.size() << ",\n"
            << "  \"checkpoint_enabled\": " << (checkpoint_enabled ? "true" : "false") << ",\n"
            << "  \"point_arrays\": [\"cell_id\",\"clone_id\",\"cell_type\",\"stage\",\"viability\",\"display_radius\"],\n"
            << "  \"effective_config\": " << config.to_json()
            << "}\n";
    });
}

}  // namespace atcg3d
