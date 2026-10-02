#include "model/hybrid_model.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

int main(int argc, char **argv) {
    try {
        std::filesystem::path yaml, checkpoint, resume, report;
        std::string mode;
        int threads = 0;
        bool dry = false;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--help") {
                std::cout << "Usage: atcg3d_hybrid --config YAML [--mode "
                             "adaptive|all_abm|all_pde] [--threads N] "
                             "[--checkpoint FILE] "
                             "[--resume-checkpoint FILE] [--report JSON] "
                             "[--dry-run]\n";
                return 0;
            }
            if (arg == "--dry-run") {
                dry = true;
                continue;
            }
            if (i + 1 == argc)
                throw std::invalid_argument("missing argument value");
            const std::string value = argv[++i];
            if (arg == "--config")
                yaml = value;
            else if (arg == "--mode")
                mode = value;
            else if (arg == "--threads")
                threads = std::stoi(value);
            else if (arg == "--checkpoint")
                checkpoint = value;
            else if (arg == "--resume-checkpoint")
                resume = value;
            else if (arg == "--report")
                report = value;
            else
                throw std::invalid_argument("unknown argument: " + arg);
        }
        auto config = atcg3d::hybrid::HybridConfig3D::load(yaml);
        if (!mode.empty())
            config.mode = mode;
        if (threads)
            config.rules.continuum.base.threads = threads;
        config.validate();
        if (dry) {
            std::cout << "{\"model\":\"" << config.model
                      << "\",\"fingerprint\":" << config.fingerprint() << "}\n";
            return 0;
        }
        atcg3d::hybrid::HybridModel3D simulation(config);
        if (!resume.empty())
            simulation.load_checkpoint(resume);
        while (simulation.step()) {
        }
        if (!checkpoint.empty())
            simulation.save_checkpoint(checkpoint);
        const auto d = simulation.diagnostics();
        std::ofstream file;
        if (!report.empty()) {
            if (!report.parent_path().empty())
                std::filesystem::create_directories(report.parent_path());
            file.open(report);
            if (!file)
                throw std::runtime_error("cannot write hybrid report");
        }
        auto &out = report.empty() ? std::cout : file;
        out << std::setprecision(17)
            << "{\"time_hours\":" << simulation.time_hours()
            << ",\"total_mass\":" << d.total_mass
            << ",\"abm_mass\":" << d.abm_mass << ",\"pde_mass\":" << d.pde_mass
            << ",\"active_mass\":" << d.active_mass
            << ",\"to_pde\":" << d.to_pde << ",\"to_abm\":" << d.to_abm
            << ",\"exchanges\":" << d.exchanges
            << ",\"state_checksum\":" << simulation.state_checksum() << "}\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "atcg3d_hybrid: " << e.what() << '\n';
        return 1;
    }
}
