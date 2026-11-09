#include "model/hybrid_model.hpp"
#include "ode/ode_model.hpp"
#include <cassert>
#include <filesystem>
#include <fstream>
#include <set>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace {
int run(const std::filesystem::path &exe,
        const std::vector<std::string> &args) {
    const auto child = fork();
    assert(child >= 0);
    if (child == 0) {
        std::string name = exe.string();
        std::vector<char *> argv{const_cast<char *>(name.c_str())};
        for (const auto &s : args)
            argv.push_back(const_cast<char *>(s.c_str()));
        argv.push_back(nullptr);
        execv(name.c_str(), argv.data());
        _exit(127);
    }
    int status;
    assert(waitpid(child, &status, 0) == child);
    return WIFEXITED(status) ? WEXITSTATUS(status) : 128;
}

void validate_registry(const std::filesystem::path& source,
                       const std::filesystem::path& entry) {
    const auto registry = YAML::LoadFile((source / "docs/configuration_registry.json").string());
    assert(registry["registry_version"].as<int>() == 1);
    std::set<std::filesystem::path> registered;
    const std::set<std::string> statuses{"recommended", "recommended_dependency",
        "benchmark", "legacy", "reproduction"};
    for (const auto& record : registry["configurations"]) {
        const std::filesystem::path path = record["path"].as<std::string>();
        assert(!path.is_absolute());
        assert(std::filesystem::is_regular_file(source / path));
        assert(registered.insert(path).second);
        assert(statuses.contains(record["status"].as<std::string>()));
    }
    std::set<std::filesystem::path> outputs;
    for (const auto* model : {"abm", "pde", "hybrid", "ode"}) {
        const std::filesystem::path path = registry["recommended"][model].as<std::string>();
        assert(registered.contains(path));
        assert(run(entry, {"--model", model, "--config", (source / path).string(),
            "--dry-run"}) == 0);
        const auto output = std::string(model) == "hybrid"
            ? atcg3d::hybrid::HybridConfig3D::load(source / path).output_directory
            : std::string(model) == "ode"
                ? atcg3d::ode::OdeConfig3D::load(source / path).output_directory
                : atcg3d::structured_pde::StructuredPdeConfig3D::load(source / path)
                    .continuum.output.directory;
        assert(!output.is_absolute());
        assert(outputs.insert(output).second);
    }
}
} // namespace
int main() {
    const auto binary = std::filesystem::path(ATCG_BINARY_DIR);
    const auto source = std::filesystem::path(ATCG_SOURCE_DIR);
    const auto entry = binary / "atcg_sim";
    validate_registry(source, entry);
    for (auto [model, path] : std::vector<std::pair<std::string, std::string>>{
             {"abm", "ATCG3D_SharedRules/config/validation_v7.yaml"},
             {"pde", "ATCG3D_StructuredPDE/config/structured_sparse_v9.yaml"},
             {"ode", "ATCG3D_ODE/config/ode_smoke_v1.yaml"},
             {"hybrid", "ATCG3D_Hybrid/config/hybrid_smoke_v1.yaml"},
             {"hybrid", "ATCG3D_Hybrid/config/hybrid_regular_cycle_v2.yaml"},
             {"hybrid", "ATCG3D_Hybrid/config/hybrid_regular_cycle_v3.yaml"}})
        assert(run(entry, {"--model", model, "--config",
                           (source / path).string(), "--dry-run"}) == 0);
    assert(run(entry,
               {"--model", "invalid", "--config",
                (source / "ATCG3D_ODE/config/ode_smoke_v1.yaml").string()}) !=
           0);
    const auto output = std::filesystem::current_path() / "migration-test";
    std::filesystem::remove_all(output);
    assert(
        run(binary / "atcg_config_migrate",
            {"--input",
             (source / "ATCG3D_SharedRules/config/validation_v7.yaml").string(),
             "--output-directory", output.string()}) == 0);
    auto migrated = atcg3d::structured_pde::StructuredPdeConfig3D::load(
        output / "model.yaml");
    assert(migrated.schema_version == 17);
    assert(migrated.continuum.schema_version == 6);
    assert(!migrated.continuum.output.directory.is_absolute());
    assert(!migrated.continuum.base.output_directory.is_absolute());
    assert(
        run(binary / "atcg_config_migrate",
            {"--input",
             (source / "ATCG3D_SharedRules/config/validation_v7.yaml").string(),
             "--output-directory", output.string()}) != 0);
    std::filesystem::remove_all(output);
    for (const auto &input : {"ATCG3D_Nutrient/config/nutrient_smoke_v1.yaml",
                              "ATCG3D_StructuredPDE_NutrientChemotaxis/config/"
                              "structured_smoke_2d_256_v5.yaml",
                              "ATCG3D_ODE/config/ode_smoke_v1.yaml",
                              "ATCG3D_Hybrid/config/hybrid_smoke_v1.yaml",
                              "ATCG3D_Hybrid/config/hybrid_regular_cycle_v2.yaml"}) {
        assert(run(binary / "atcg_config_migrate",
                   {"--input", (source / input).string(), "--output-directory",
                    output.string(), "--grid-edge", "128",
                    "--nutrient-K-per-cell-hour", "0.01"}) == 0);
        std::filesystem::remove_all(output);
    }
    assert(run(entry,
               {"--model", "hybrid", "--config",
                (source / "ATCG3D_Hybrid/config/hybrid_smoke_v1.yaml").string(),
                "--output-root", output.string(), "--report",
                (output / "summary.json").string()}) == 0);
    const auto metrics = output / "atcg3d_hybrid_smoke_v1_run" / "metrics.csv";
    std::ifstream csv(metrics);
    std::string header, initial, final;
    assert(std::getline(csv, header));
    assert(header == "time_hours,total_mass,abm_mass,pde_mass,active_mass,"
                     "mean_nutrient,to_pde,to_abm,state_checksum");
    assert(std::getline(csv, initial) && initial.starts_with("0,"));
    for (std::string row; std::getline(csv, row);)
        final = row;
    assert(final.starts_with("8,"));
    assert(std::filesystem::exists(output / "summary.json"));
    std::filesystem::remove_all(output);

    // Native dispatch and the ensemble adapter must import the same resource-
    // initialized agents, including initial activation and sampled clocks.
    std::filesystem::create_directories(output);
    auto structured = YAML::LoadFile((source / "ATCG3D_SharedRules/config/active_r200_ci_v12.yaml").string());
    auto continuum = YAML::LoadFile((source / "ATCG3D_SharedRules/config/continuum_active_r200_ci_v8.yaml").string());
    auto base = YAML::LoadFile((source / "configs/shared_activation_r200_v1.yaml").string());
    base["initial"]["explicit_counts"]["r_cells"] = 8;
    base["initial"]["explicit_counts"]["K_cells"] = 8;
    base["initial"]["geometry"]["outer_radius"] = 4;
    base["initial"]["geometry"]["shell_inner_radius"] = 3;
    base["initial"]["geometry"]["inner_small_radius"] = 3;
    base["migration"]["activation_threshold"] = 0.001;
    base["simulation"]["threads"] = 1;
    continuum["base_config"] = "base.yaml";
    continuum["grid"]["shape"] = std::vector<int>{12, 12, 1};
    continuum["grid"]["origin"] = std::vector<double>{-6, -6, -0.5};
    continuum["numerics"]["end_time_hours"] = 0.01;
    continuum["numerics"]["time_step_hours"] = 0.01;
    continuum["vascular"]["source_mode"] = "static_voxels";
    continuum["vascular"]["static_sources"] = YAML::Load("[]");
    structured["continuum_config"] = "continuum.yaml";
    structured["structured_migration"]["reactivation_density_threshold"] = 0.0005;
    structured["output"]["directory"] = "native_shared_run";
    for (auto [name, node] : {std::pair{"base.yaml", base}, std::pair{"continuum.yaml", continuum},
                              std::pair{"structured.yaml", structured}}) {
        std::ofstream yaml(output / name);
        yaml << node;
    }
    assert(run(entry, {"--model", "pde", "--config", (output / "structured.yaml").string(),
                       "--output-root", output.string()}) == 0);
    assert(run(binary / "atcg3d_shared_abm", {"--model", "pde", "--config",
        (output / "structured.yaml").string(), "--report", (output / "shared.json").string()}) == 0);
    const auto native = YAML::LoadFile((output / "native_shared_run/final.json").string());
    const auto shared = YAML::LoadFile((output / "shared.json").string());
    assert(native["state_checksum"].as<std::uint64_t>() == shared["state_checksum"].as<std::uint64_t>());
    structured["division_clock"]["model"] = "transported_shifted_geometric_v1";
    structured["division_clock"]["work_bin_width"] = 0.5;
    structured["division_clock"]["maximum_work"] = 128.0;
    {
        std::ofstream yaml(output / "structured.yaml");
        yaml << structured;
    }
    auto hybrid = YAML::LoadFile((source / "ATCG3D_Hybrid/config/hybrid_regular_cycle_v2.yaml").string());
    hybrid["structured_config"] = "structured.yaml";
    {
        std::ofstream yaml(output / "hybrid.yaml");
        yaml << hybrid;
    }
    assert(run(entry, {"--model", "hybrid", "--config", (output / "hybrid.yaml").string(),
        "--mode", "all_pde", "--seed", "2", "--step-hours", "0.01", "--exchange-hours", "0.01",
        "--no-output", "--validation-report", "--report", (output / "hybrid.json").string()}) == 0);
    const auto mixed = YAML::LoadFile((output / "hybrid.json").string());
    assert(mixed["r_mass"].as<double>() > 0 && mixed["K_mass"].as<double>() > 0);
    assert(mixed["radial_mass"].IsSequence());
    assert(!std::filesystem::exists("atcg3d_hybrid_regular_cycle_v2_run"));
    std::filesystem::remove_all(output);
}
