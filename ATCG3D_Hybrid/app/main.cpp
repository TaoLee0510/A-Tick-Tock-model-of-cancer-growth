#include "config/output_paths.hpp"
#include "model/hybrid_model.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <numeric>

int main(int argc, char **argv) {
    try {
        std::filesystem::path yaml, checkpoint, resume, report, output_root;
        std::string mode;
        int threads = 0;
        bool dry = false, no_output = false, validation_report = false;
        std::optional<std::uint64_t> seed;
        std::optional<double> step_hours, exchange_hours;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--help") {
                std::cout << "Usage: atcg3d_hybrid --config YAML [--mode "
                             "adaptive|all_abm|all_pde] [--threads N] "
                             "[--checkpoint FILE] "
                             "[--resume-checkpoint FILE] [--report JSON] "
                             "[--output-root PATH] [--dry-run] [--seed N] "
                             "[--step-hours H] [--exchange-hours H] "
                             "[--no-output] [--validation-report]\n";
                return 0;
            }
            if (arg == "--dry-run") {
                dry = true;
                continue;
            }
            if (arg == "--no-output") { no_output = true; continue; }
            if (arg == "--validation-report") { validation_report = true; continue; }
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
            else if (arg == "--output-root")
                output_root = value;
            else if (arg == "--report")
                report = value;
            else if (arg == "--seed") seed = std::stoull(value);
            else if (arg == "--step-hours") step_hours = std::stod(value);
            else if (arg == "--exchange-hours") exchange_hours = std::stod(value);
            else
                throw std::invalid_argument("unknown argument: " + arg);
        }
        auto config = atcg3d::hybrid::HybridConfig3D::load(yaml);
        if (!mode.empty())
            config.mode = mode;
        if (threads)
            config.rules.continuum.base.threads = threads;
        if (seed) config.rules.continuum.base.seed = *seed;
        if (step_hours) {
            config.rules.continuum.time_step_hours = *step_hours;
            config.rules.continuum.nutrient.refresh_every_hours = *step_hours;
        }
        if (exchange_hours) config.exchange_every_hours = *exchange_hours;
        config.validate();
        if (dry) {
            std::cout << "{\"model\":\"" << config.model
                      << "\",\"fingerprint\":" << config.fingerprint() << "}\n";
            return 0;
        }
        const auto directory = atcg3d::resolve_output_directory(
            config.output_directory, output_root);
        if (!no_output && resume.empty() && std::filesystem::exists(directory) &&
            !std::filesystem::is_empty(directory))
            throw std::runtime_error(
                "refusing to overwrite hybrid output directory");
        std::ofstream metrics;
        if (!no_output) {
            std::filesystem::create_directories(directory);
            metrics.open(directory / (resume.empty() ? "metrics.csv" : "metrics_resumed.csv"));
        }
        if (!no_output && !metrics)
            throw std::runtime_error("cannot write hybrid metrics");
        metrics << "time_hours,total_mass,abm_mass,pde_mass,active_mass,mean_"
                   "nutrient,to_pde,to_abm,state_checksum\n"
                << std::setprecision(17);
        atcg3d::hybrid::HybridModel3D simulation(config);
        if (!resume.empty())
            simulation.load_checkpoint(resume);
        const auto observe = [&] {
            if (no_output) return;
            const auto d = simulation.diagnostics();
            metrics << simulation.time_hours() << ',' << d.total_mass << ','
                    << d.abm_mass << ',' << d.pde_mass << ',' << d.active_mass
                    << ',' << d.mean_nutrient << ',' << d.to_pde << ','
                    << d.to_abm << ',' << simulation.state_checksum() << '\n';
        };
        if (resume.empty())
            simulation.initialize();
        observe();
        while (simulation.step())
            observe();
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
            << ",\"state_checksum\":" << simulation.state_checksum();
        if (validation_report || config.schema_version >= 2) {
            const auto radial = simulation.radial_mass();
            const double total = std::accumulate(radial.begin(), radial.end(), 0.0);
            const auto quantile = [&](double fraction) {
                if (!(total > 0.0)) return 0.0;
                double sum = 0;
                for (std::size_t i = 0; i < radial.size(); ++i) {
                    sum += radial[i];
                    if (sum >= total * fraction) return double(i + 1);
                }
                return 0.0;
            };
            out << ",\"r_mass\":" << d.r_mass << ",\"K_mass\":" << d.K_mass
                << ",\"r_K_ratio\":" << (d.K_mass > 0 ? d.r_mass / d.K_mass : 0.0)
                << ",\"active_fraction\":" << (d.r_mass > 0 ? d.active_mass / d.r_mass : 0.0)
                << ",\"r50\":" << quantile(0.5) << ",\"r90\":" << quantile(0.9)
                << ",\"r99\":" << quantile(0.99) << ",\"radial_mass\":[";
            for (std::size_t i = 0; i < radial.size(); ++i) { if (i) out << ','; out << radial[i]; }
            out << ']';
        }
        out << "}\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "atcg3d_hybrid: " << e.what() << '\n';
        return 1;
    }
}
