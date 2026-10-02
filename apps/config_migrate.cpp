#include "config/continuum_config.hpp"
#include "config/model_config.hpp"
#include "config/nutrient_config.hpp"
#include "config/structured_config.hpp"
#include "model/hybrid_model.hpp"
#include "ode/ode_model.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <yaml-cpp/yaml.h>

namespace {
struct Migrator {
    std::filesystem::path directory;
    std::set<std::filesystem::path> active;
    int grid_edge{};
    double nutrient_K_per_cell{-1};
    YAML::Node upgrade(YAML::Node y) {
        const auto name = y["schema"]["name"].as<std::string>();
        const int old = y["schema"]["version"].as<int>();
        int latest = old;
        if (name == "atcg3d.model_config")
            latest = 3;
        else if (name == "atcg3d.nutrient_model_config")
            latest = 3;
        else if (name == "atcg3d.continuum_model_config") {
            if (old < 3)
                throw std::invalid_argument(
                    "continuum v1-v2 migration requires explicit "
                    "shared-resource parameter calibration");
            latest = 6;
            if (!y["angiogenesis"])
                y["angiogenesis"]["model"] = "disabled";
        } else if (name == "atcg3d.structured_pde_config") {
            if (old < 5)
                throw std::invalid_argument(
                    "structured v1-v4 migration requires explicit "
                    "shared-resource parameter calibration");
            latest = 12;
            y["structured_migration"]["activation_stop"] =
                "cohort_clock_refractory_hysteresis_v3";
            if (!y["storage"]) {
                y["storage"]["model"] = "dense_v1";
                y["storage"]["maximum_active_voxels"] = 1000000;
            }
        } else if (name != "atcg3d.ode_config" &&
                   name != "atcg3d.hybrid_config")
            throw std::invalid_argument("unsupported configuration schema");
        if (old < 1 || old > latest)
            throw std::invalid_argument("unsupported source schema version");
        y["schema"]["version"] = latest;
        return y;
    }
    std::filesystem::path copy(const std::filesystem::path &input,
                               const std::string &filename) {
        const auto source = std::filesystem::canonical(input);
        if (!active.insert(source).second)
            throw std::invalid_argument("cyclic configuration references");
        auto y = upgrade(YAML::LoadFile(source.string()));
        for (auto key :
             {"base_config", "continuum_config", "structured_config"})
            if (y[key]) {
                const auto child =
                    source.parent_path() / y[key].as<std::string>();
                const auto target =
                    std::filesystem::path(filename).stem().string() + "_" +
                    key + ".yaml";
                copy(child, target);
                y[key] = target;
            }
        if (y["run"] && y["run"]["mode"]) {
            y["run"]["mode"] = "new";
            y["run"]["resume_checkpoint"] = YAML::Node(YAML::NodeType::Null);
        }
        if (y["initialization"] && y["initialization"]["abm_checkpoint"]) {
            y["initialization"]["mode"] = "base_model";
            y["initialization"]["abm_checkpoint"] =
                YAML::Node(YAML::NodeType::Null);
        }
        if (y["schema"]["name"].as<std::string>() ==
                "atcg3d.nutrient_model_config" &&
            y["nutrient"]["consumption"]["r_per_voxel_hour"]) {
            if (nutrient_K_per_cell <= 0)
                throw std::invalid_argument(
                    "nutrient v1 migration requires --nutrient-K-per-cell-hour "
                    "to select the new uptake contract");
            auto old = y["nutrient"]["consumption"];
            const double r = old["r_per_voxel_hour"].as<double>(),
                         K = old["K_per_voxel_hour"].as<double>();
            if (K <= 0 || r <= 0 ||
                old["r_half_saturation"].as<double>() !=
                    old["K_half_saturation"].as<double>())
                throw std::invalid_argument(
                    "nutrient v1 unequal uptake saturations/rates require "
                    "explicit calibration");
            YAML::Node consumption;
            consumption["K_per_cell_hour"] = nutrient_K_per_cell;
            consumption["r_to_K_ratio"] = r / K;
            consumption["half_saturation"] =
                old["K_half_saturation"].as<double>();
            y["nutrient"]["consumption"] = consumption;
        }
        if (y["schema"]["name"].as<std::string>() ==
                "atcg3d.nutrient_model_config" &&
            !y["vascular_geometry"]) {
            if (grid_edge < 2)
                throw std::invalid_argument(
                    "old nutrient migration requires --grid-edge for the new "
                    "finite resource domain");
            auto base = atcg3d::Model3DConfig::load(
                directory / y["base_config"].as<std::string>());
            auto geometry = y["vascular_geometry"];
            geometry["shape"].push_back(grid_edge);
            geometry["shape"].push_back(grid_edge);
            geometry["shape"].push_back(base.thin_layer ? 1 : grid_edge);
            geometry["origin"].push_back(-grid_edge / 2.0);
            geometry["origin"].push_back(-grid_edge / 2.0);
            geometry["origin"].push_back(base.thin_layer ? -0.5
                                                         : -grid_edge / 2.0);
            geometry["spacing_voxels"] = 1;
            geometry["source_mode"] = "static_voxels";
            geometry["synthetic_axis"] = "y";
            for (int i = 0; i < 3; ++i)
                geometry["synthetic_center"].push_back(0.0);
            geometry["synthetic_radius_voxels"] = 1;
            geometry["static_sources"] = YAML::Node(YAML::NodeType::Sequence);
            const int halo = y["nutrient"]["halo_voxels"].as<int>();
            y["nutrient"]["halo_voxels"] =
                std::max(halo, base.direction_density_radius);
        }
        if (y["nutrient"] && y["nutrient"]["output"] &&
            y["nutrient"]["output"]["directory"])
            y["nutrient"]["output"]["directory"] =
                std::filesystem::path(filename).stem().string() +
                "_nutrient_run";
        for (auto key : {"output", "run"})
            if (y[key] && y[key]["directory"])
                y[key]["directory"] =
                    std::filesystem::path(filename).stem().string() + "_run";
        if (y["output"] && y["output"]["directory"])
            y["output"]["directory"] =
                std::filesystem::path(filename).stem().string() + "_run";
        const auto path = directory / filename;
        std::ofstream out(path);
        out << y << '\n';
        if (!out)
            throw std::runtime_error("migration write failed");
        active.erase(source);
        return path;
    }
};
} // namespace
int main(int argc, char **argv) {
    try {
        std::filesystem::path input, directory;
        int grid_edge = 0;
        double nutrient_K = -1;
        for (int i = 1; i < argc; ++i) {
            const std::string key = argv[i];
            if (key == "--help") {
                std::cout << "Usage: atcg_config_migrate --input YAML "
                             "--output-directory NEW_DIR [--grid-edge "
                             "N] [--nutrient-K-per-cell-hour RATE]\n"
                             "Copies referenced configs and upgrades "
                             "supported shared contracts to latest schemas.\n";
                return 0;
            }
            if (++i == argc)
                throw std::invalid_argument("missing option value");
            if (key == "--input")
                input = argv[i];
            else if (key == "--output-directory")
                directory = argv[i];
            else if (key == "--nutrient-K-per-cell-hour")
                nutrient_K = std::stod(argv[i]);
            else if (key == "--grid-edge")
                grid_edge = std::stoi(argv[i]);
            else
                throw std::invalid_argument("unknown migration option");
        }
        if (input.empty() || directory.empty() ||
            std::filesystem::exists(directory))
            throw std::invalid_argument(
                "migration requires input and a new output directory");
        std::filesystem::create_directories(directory);
        const auto path = Migrator{directory, {}, grid_edge, nutrient_K}.copy(
            input, "model.yaml");
        const auto name =
            YAML::LoadFile(path.string())["schema"]["name"].as<std::string>();
        if (name == "atcg3d.model_config")
            (void)atcg3d::Model3DConfig::load(path);
        else if (name == "atcg3d.nutrient_model_config")
            (void)atcg3d::nutrient::NutrientModelConfig3D::load(path);
        else if (name == "atcg3d.continuum_model_config")
            (void)atcg3d::continuum::ContinuumModelConfig3D::load(path);
        else if (name == "atcg3d.structured_pde_config")
            (void)atcg3d::structured_pde::StructuredPdeConfig3D::load(path);
        else if (name == "atcg3d.ode_config")
            (void)atcg3d::ode::OdeConfig3D::load(path);
        else if (name == "atcg3d.hybrid_config")
            (void)atcg3d::hybrid::HybridConfig3D::load(path);
        std::ofstream note(directory / "MIGRATION.md");
        note
            << "# Configuration migration\n\nSource files were left unchanged. "
               "References and output directories are relative.\nNew schemas "
               "use exact ABM growth windows and transported cohort refractory "
               "clocks;\nupgrading v5/v6 changes the old grid-local cooldown "
               "closure. Restart a new run;\nold checkpoint fingerprints are "
               "intentionally not reused by changed PDE contracts.\n"
               "Nutrient v1 upgrades use the explicitly supplied per-cell K "
               "uptake rate\nand preserve the old r/K uptake ratio and common "
               "saturation. This is a new\ncalibration choice, not an exact "
               "conversion of occupied-voxel uptake.\n";
        std::cout << path.string() << '\n';
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "atcg_config_migrate: " << e.what() << '\n';
        return 1;
    }
}
