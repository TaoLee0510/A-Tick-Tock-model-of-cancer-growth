#include <algorithm>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "config/continuum_config.hpp"
#include "config/nutrient_config.hpp"
#include "config/output_paths.hpp"
#include "config/structured_config.hpp"
#include "engine/simulation.hpp"
#include "io/nutrient_output.hpp"

namespace {

std::string read_text(const std::filesystem::path& path) {
    std::ifstream stream(path, std::ios::binary);
    assert(stream);
    return {std::istreambuf_iterator<char>(stream), {}};
}

std::vector<std::string> columns(const std::string& line) {
    std::istringstream stream(line);
    std::vector<std::string> result;
    for (std::string value; std::getline(stream, value, ',');) {
        result.push_back(value);
    }
    return result;
}

void check_metrics(atcg3d::nutrient::NutrientModelConfig3D config,
                   const std::filesystem::path& directory,
                   bool legacy) {
    using namespace atcg3d;
    using namespace atcg3d::nutrient;
    std::filesystem::remove_all(directory);
    std::filesystem::create_directories(directory / "nutrient");
    const auto path = directory / "nutrient/metrics.csv";
    const std::string legacy_header =
        "time_hours,completed_events,refresh_count,blocks,active_voxels,"
        "source_voxels,consuming_voxels,minimum,maximum,mean,"
        "last_max_update,field_checksum,base_checksum";
    const std::string header = legacy ? legacy_header :
        "time_hours,completed_events,refresh_count,blocks,active_voxels,"
        "source_voxels,consuming_voxels,assembled_r_consumption_per_hour,"
        "assembled_K_consumption_per_hour,minimum,maximum,mean,"
        "last_max_update,field_checksum,base_checksum";
    const std::string original = header + "\r\n" +
        (legacy ? "0,0,0,0,0,0,0,0,0,0,0,0,0\n" :
                  "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0\n");
    std::ofstream(path, std::ios::binary) << original;
    config.base.output_directory = directory;
    config.base.output_enabled = true;
    config.nutrient.metrics_every_hours = 0.001;
    config.nutrient.field_snapshot_every_hours = 0.0;
    auto base = config.base;
    base.output_enabled = false;
    base.control_enabled = false;
    base.angiogenesis.enabled = false;
    base.end_time_hours = 1.0;
    base.max_events = 16;
    Simulation3D simulation(base);
    simulation.initialize();
    simulation.run();
    assert(simulation.clock().time_hours > 0.001);
    NutrientEnvironment3D environment(config.nutrient, base.thin_layer);
    {
        NutrientOutput3D output(config, environment, 0.0);
        output.observe(simulation);
    }
    const std::string text = read_text(path);
    assert(text.starts_with(original));
    const auto values = columns(text.substr(original.size()));
    assert(values.size() == (legacy ? 13 : 15));
    assert(std::stod(values[0]) == simulation.clock().time_hours);
    assert(std::stoull(values[values.size() - 2]) == environment.field_checksum());
    assert(std::stoull(values.back()) == simulation.state_checksum());

    // Reject unknown layouts before opening the stream for append.
    const std::string unknown = "unexpected,columns\n1,2\n";
    std::ofstream(path, std::ios::binary | std::ios::trunc) << unknown;
    bool rejected = false;
    try {
        NutrientOutput3D output(config, environment, 0.0);
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    assert(rejected);
    assert(read_text(path) == unknown);
    std::filesystem::remove_all(directory);
}

}  // namespace

int main() {
    using namespace atcg3d;
    const std::filesystem::path source = ATCG_SOURCE_DIR;
    const auto base = Model3DConfig::load(
        source / "configs/atcg3d_smoke_test_v3.yaml");
    const auto nutrient = nutrient::NutrientModelConfig3D::load(
        source / "ATCG3D_Nutrient/config/nutrient_smoke_v1.yaml");
    const auto continuum = continuum::ContinuumModelConfig3D::load(
        source / "ATCG3D_Continuum/config/continuum_smoke_v1.yaml");
    assert(base.dynamics_json() == read_text(
        source / "tests/3d/fixtures/legacy_smoke_dynamics.json"));
    std::ifstream fingerprints(source / "tests/3d/fixtures/legacy_smoke_fingerprints.txt");
    std::uint64_t nutrient_fingerprint{}, continuum_fingerprint{};
    fingerprints >> nutrient_fingerprint >> continuum_fingerprint;
    assert(fingerprints);
    assert(nutrient.nutrient.fingerprint() == nutrient_fingerprint);
    assert(continuum.dynamics_fingerprint() == continuum_fingerprint);

    auto changed = base;
    changed.direction_guidance_model = "low_density_high_resource_v1";
    assert(changed.dynamics_json() != base.dynamics_json());
    changed = base;
    changed.direction_density_guidance_exponent = 2.0;
    assert(changed.dynamics_json() != base.dynamics_json());
    changed = base;
    changed.direction_resource_guidance_exponent = 2.0;
    assert(changed.dynamics_json() != base.dynamics_json());
    changed = base;
    changed.direction_minimum_guidance_weight = 0.01;
    assert(changed.dynamics_json() != base.dynamics_json());
    changed = base;
    changed.activated_r_normal_multiplier = 2.0;
    assert(changed.dynamics_json() != base.dynamics_json());
    changed.activated_r_migration_rate_model = "normal_multiplier";
    assert(changed.dynamics_json().find("activated_r_normal_multiplier") !=
           std::string::npos);
    changed.activated_r_normal_multiplier = 1.0;
    assert(changed.dynamics_json().find("activated_r_normal_multiplier") !=
           std::string::npos);

    auto changed_nutrient = nutrient.nutrient;
    changed_nutrient.consumption_model = "per_cell_ratio_v2";
    assert(changed_nutrient.fingerprint() != nutrient_fingerprint);
    auto old_continuum = continuum;
    old_continuum.nutrient.initial_value = 0.5;
    old_continuum.nutrient.growth_half_saturation = 0.5;
    old_continuum.nutrient.common_density_limit = 12.0;
    old_continuum.nutrient.common_carrying_capacity = 24.0;
    old_continuum.nutrient.boundary_mode = "ignored_legacy_field";
    old_continuum.nutrient.consumption_model = "ignored_legacy_field";
    old_continuum.nutrient.tumor_front_density_threshold = 0.2;
    assert(old_continuum.dynamics_fingerprint() == continuum_fingerprint);
    auto v2_continuum = continuum;
    v2_continuum.schema_version = 2;
    const auto v2_fingerprint = v2_continuum.dynamics_fingerprint();
    v2_continuum.nutrient.initial_value = 0.5;
    v2_continuum.nutrient.common_density_limit = 12.0;
    v2_continuum.nutrient.boundary_mode = "ignored_legacy_field";
    assert(v2_continuum.dynamics_fingerprint() == v2_fingerprint);
    v2_continuum.nutrient.consumption_model = "per_cell_ratio_v2";
    assert(v2_continuum.dynamics_fingerprint() != v2_fingerprint);
    auto v3_continuum = continuum;
    v3_continuum.schema_version = 3;
    const auto v3_fingerprint = v3_continuum.dynamics_fingerprint();
    v3_continuum.nutrient.initial_value = 0.5;
    assert(v3_continuum.dynamics_fingerprint() != v3_fingerprint);

    // Every shipped configuration has a relative, distinct effective output
    // directory, including wrappers that share the same biological base.
    std::set<std::filesystem::path> directories;
    for (const auto* folder : {"configs", "ATCG3D_Nutrient/config",
             "ATCG3D_Continuum/config", "ATCG3D_StructuredPDE/config",
             "ATCG3D_StructuredPDE_NutrientChemotaxis/config"}) {
        for (const auto& entry : std::filesystem::directory_iterator(source / folder)) {
            if (entry.path().extension() != ".yaml") continue;
            const auto yaml = YAML::LoadFile(entry.path().string());
            const std::string schema = yaml["schema"]["name"].as<std::string>();
            std::filesystem::path directory;
            if (schema == "atcg3d.model_config") {
                directory = Model3DConfig::load(entry.path()).output_directory;
            } else if (schema == "atcg3d.nutrient_model_config") {
                directory = nutrient::NutrientModelConfig3D::load(entry.path())
                    .base.output_directory;
            } else if (schema == "atcg3d.continuum_model_config") {
                directory = continuum::ContinuumModelConfig3D::load(entry.path())
                    .output.directory;
            } else {
                assert(schema == "atcg3d.structured_pde_config");
                directory = structured_pde::StructuredPdeConfig3D::load(entry.path())
                    .continuum.output.directory;
            }
            assert(!directory.empty() && directory.is_relative());
            assert(directories.insert(directory).second);
        }
    }
    const auto root = std::filesystem::temp_directory_path() / "atcg3d_output_root";
    assert(resolve_output_directory("model/run", root) == root / "model/run");
    assert(resolve_output_directory("model/run", {}) == "model/run");
    assert(resolve_output_directory(root / "absolute_run", "other") ==
           root / "absolute_run");
    changed = base;
    changed.output_directory = resolve_output_directory(base.output_directory, root);
    assert(changed.dynamics_json() == base.dynamics_json());
    check_metrics(nutrient, root / "legacy", true);
    check_metrics(nutrient, root / "current", false);
    std::filesystem::remove_all(root);
}
