#include <cerrno>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace {
std::filesystem::path executable_directory(const char *name) {
    std::filesystem::path path(name);
    if (path.has_parent_path())
        return std::filesystem::canonical(path).parent_path();
    const char *search = std::getenv("PATH");
    std::string paths = search ? search : "";
    for (std::size_t begin = 0; begin <= paths.size();) {
        const auto end = paths.find(':', begin);
        auto candidate =
            std::filesystem::path(paths.substr(begin, end - begin)) / path;
        if (std::filesystem::exists(candidate))
            return std::filesystem::canonical(candidate).parent_path();
        if (end == std::string::npos)
            break;
        begin = end + 1;
    }
    throw std::runtime_error("cannot locate simulator binaries");
}
} // namespace
int main(int argc, char **argv) {
    try {
        std::string model, config;
        std::vector<std::string> forwarded;
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "--help" || arg == "-h") {
                std::cout
                    << "Usage: atcg_sim --model abm|ode|pde|hybrid --config "
                       "YAML [model options]\n"
                    << "ABM accepts native, nutrient or shared structured "
                       "YAML; PDE accepts continuum or structured YAML.\n"
                    << "ODE and hybrid accept their wrapper YAML. All models "
                       "support --dry-run.\n"
                    << "Use docs/simulator.md for checkpoint, output, "
                       "migration and field options.\n";
                return 0;
            }
            if (arg == "--model") {
                if (++i == argc)
                    throw std::invalid_argument("missing model");
                model = argv[i];
            } else {
                forwarded.push_back(arg);
                if (arg == "--config") {
                    if (++i == argc)
                        throw std::invalid_argument("missing config");
                    config = argv[i];
                    forwarded.push_back(config);
                }
            }
        }
        if (config.empty())
            throw std::invalid_argument("--config is required");
        const auto root = YAML::LoadFile(config);
        const auto schema = root["schema"]["name"].as<std::string>();
        std::string binary;
        if (model == "abm") {
            if (schema == "atcg3d.model_config")
                binary = "atcg3d";
            else if (schema == "atcg3d.nutrient_model_config")
                binary = "atcg3d_nutrient";
            else if (schema == "atcg3d.structured_pde_config") {
                binary = "atcg3d_shared_abm";
                forwarded.insert(forwarded.end(), {"--model", "abm"});
            }
        } else if (model == "pde") {
            if (schema == "atcg3d.continuum_model_config")
                binary = "atcg3d_continuum";
            else if (schema == "atcg3d.structured_pde_config")
                binary = "atcg3d_structured_pde";
        } else if (model == "ode" && schema == "atcg3d.ode_config")
            binary = "atcg3d_ode";
        else if (model == "hybrid" && schema == "atcg3d.hybrid_config")
            binary = "atcg3d_hybrid";
        if (binary.empty())
            throw std::invalid_argument("model/configuration schema mismatch");
        const auto path = (executable_directory(argv[0]) / binary).string();
        std::vector<char *> arguments{const_cast<char *>(path.c_str())};
        for (auto &arg : forwarded)
            arguments.push_back(arg.data());
        arguments.push_back(nullptr);
        execv(path.c_str(), arguments.data());
        throw std::runtime_error("cannot execute " + binary + ": " +
                                 std::to_string(errno));
    } catch (const std::exception &e) {
        std::cerr << "atcg_sim: " << e.what() << '\n';
        return 1;
    }
}
