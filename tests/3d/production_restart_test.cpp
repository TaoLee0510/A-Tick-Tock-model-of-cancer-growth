#include <cassert>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace {
void run(const std::filesystem::path& executable,
         const std::vector<std::string>& arguments,
         const std::filesystem::path& report) {
    const auto child = fork();
    assert(child >= 0);
    if (child == 0) {
        const auto* file = std::freopen(report.c_str(), "w", stdout);
        assert(file != nullptr);
        std::string name = executable.string();
        std::vector<char*> argv{const_cast<char*>(name.c_str())};
        for (const auto& value : arguments) {
            argv.push_back(const_cast<char*>(value.c_str()));
        }
        argv.push_back(nullptr);
        execv(name.c_str(), argv.data());
        _exit(127);
    }
    int status{};
    assert(waitpid(child, &status, 0) == child);
    assert(WIFEXITED(status) && WEXITSTATUS(status) == 0);
}

void write(const std::filesystem::path& path, const YAML::Node& node) {
    std::ofstream output(path);
    output << node << '\n';
    assert(output.good());
}
}  // namespace

int main() {
    const std::filesystem::path source = ATCG_SOURCE_DIR;
    const auto directory = std::filesystem::current_path() / "production-restart-fixture";
    std::filesystem::remove_all(directory);
    std::filesystem::create_directories(directory);
    auto base = YAML::LoadFile((source / "configs/production_2d_2000_r200_v5.yaml").string());
    base["initial"]["explicit_counts"]["r_cells"] = 32;
    base["initial"]["explicit_counts"]["K_cells"] = 32;
    base["initial"]["geometry"]["outer_radius"] = 8;
    base["initial"]["geometry"]["shell_inner_radius"] = 6;
    base["initial"]["geometry"]["inner_small_radius"] = 6;
    base["migration"]["activation_threshold"] = 0.001;
    base["migration"]["activated_r_rate"]["multiplier"] = 20.0;
    auto continuum = YAML::LoadFile((source / "ATCG3D_SharedRules/config/continuum_production_2d_2000_r200_v8.yaml").string());
    continuum["base_config"] = "base.yaml";
    continuum["grid"]["shape"] = std::vector<int>{256, 256, 1};
    continuum["grid"]["origin"] = std::vector<double>{-128, -128, -0.5};
    continuum["migration"]["activated_r_mobility_multiplier"] = 20.0;
    auto rules = YAML::LoadFile((source / "ATCG3D_SharedRules/config/production_pde_2d_2000_r200_v17.yaml").string());
    rules["continuum_config"] = "continuum.yaml";
    rules["structured_migration"]["reactivation_density_threshold"] = 0.0005;
    rules["storage"]["maximum_active_voxels"] = 1000000;
    auto hybrid = YAML::LoadFile((source / "ATCG3D_Hybrid/config/production_2d_2000_r200_v5.yaml").string());
    hybrid["structured_config"] = "rules.yaml";
    write(directory / "base.yaml", base);
    write(directory / "continuum.yaml", continuum);
    write(directory / "rules.yaml", rules);
    write(directory / "hybrid.yaml", hybrid);
    auto legacy_hybrid = YAML::Clone(hybrid);
    legacy_hybrid["schema"]["version"] = 4;
    legacy_hybrid["model"] = "hybrid_invasion_front_v4";
    write(directory / "legacy_hybrid.yaml", legacy_hybrid);
    std::vector<std::string> models{"pde", "hybrid", "legacy_hybrid"};
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
    models.push_back("abm");
#endif
    for (const auto& model : models) {
        const auto executable = std::filesystem::path(ATCG_BINARY_DIR) / "atcg3d_production_benchmark";
        const auto config = directory / (model == "hybrid" ? "hybrid.yaml"
            : model == "legacy_hybrid" ? "legacy_hybrid.yaml" : "rules.yaml");
        const auto checkpoint = directory / (model + ".checkpoint");
        const std::vector<std::string> common{"--config", config.string(), "--model",
            model == "legacy_hybrid" ? "hybrid" : model};
        auto arguments = common;
        arguments.insert(arguments.end(), {"--stop-hours", "0.5", "--threads", "1"});
        run(executable, arguments, directory / "reference.json");
        arguments = common;
        arguments.insert(arguments.end(), {"--stop-hours", "0.25", "--threads", "1",
            "--checkpoint", checkpoint.string()});
        run(executable, arguments, directory / "prefix.json");
        arguments = common;
        arguments.insert(arguments.end(), {"--stop-hours", "0.5", "--threads", "8",
            "--resume", checkpoint.string()});
        run(executable, arguments, directory / "resumed.json");
        const auto reference = YAML::LoadFile((directory / "reference.json").string());
        const auto prefix = YAML::LoadFile((directory / "prefix.json").string());
        const auto resumed = YAML::LoadFile((directory / "resumed.json").string());
        for (const auto* key : {"state_checksum", "field_checksum"}) {
            assert(reference[key].as<std::uint64_t>() == resumed[key].as<std::uint64_t>());
        }
        for (const auto* key : {"mass", "active_mass", "abm_fraction", "radius_99",
                               "vascular_roots", "vascular_length"}) {
            assert(reference[key].as<double>() == resumed[key].as<double>());
        }
        assert(reference["time_hours"].as<double>() == 0.5);
        assert(prefix["time_hours"].as<double>() == 0.25);
        assert(reference["configured_hours"].as<double>() == 2160.0);
        assert(prefix["checkpoint_bytes"].as<std::uint64_t>() > 0);
        assert(reference["active_mass"].as<double>() > 0.0);
        if (model != "legacy_hybrid") {
            assert(reference["vascular_roots"].as<double>() > 0.0);
        }
    }
    std::filesystem::remove_all(directory);
}
