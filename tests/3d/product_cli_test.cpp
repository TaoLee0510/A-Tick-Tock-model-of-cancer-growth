#include "model/hybrid_model.hpp"
#include "ode/ode_model.hpp"
#include <cassert>
#include <filesystem>
#include <fstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

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
} // namespace
int main() {
    const auto binary = std::filesystem::path(ATCG_BINARY_DIR);
    const auto source = std::filesystem::path(ATCG_SOURCE_DIR);
    const auto entry = binary / "atcg_sim";
    for (auto [model, path] : std::vector<std::pair<std::string, std::string>>{
             {"abm", "ATCG3D_SharedRules/config/validation_v7.yaml"},
             {"pde", "ATCG3D_StructuredPDE/config/structured_sparse_v9.yaml"},
             {"ode", "ATCG3D_ODE/config/ode_smoke_v1.yaml"},
             {"hybrid", "ATCG3D_Hybrid/config/hybrid_smoke_v1.yaml"}})
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
    assert(migrated.schema_version == 10);
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
                              "ATCG3D_Hybrid/config/hybrid_smoke_v1.yaml"}) {
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
}
