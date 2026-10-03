#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>

#include "model/shared_angiogenesis.hpp"

namespace {
struct Arguments {
    int edge{2000}, threads{1}, source_edge{128};
    double hours{24.0}, step{0.25};
};

Arguments arguments(int argc, char** argv) {
    Arguments result;
    for (int index = 1; index < argc; ++index) {
        const std::string option = argv[index];
        if (index + 1 >= argc) throw std::invalid_argument("benchmark option requires a value");
        const std::string value = argv[++index];
        if (option == "--edge") result.edge = std::stoi(value);
        else if (option == "--threads") result.threads = std::stoi(value);
        else if (option == "--source-edge") result.source_edge = std::stoi(value);
        else if (option == "--hours") result.hours = std::stod(value);
        else if (option == "--step-hours") result.step = std::stod(value);
        else throw std::invalid_argument("unknown vascular benchmark option");
    }
    if (result.edge < 8 || result.source_edge < 1 || result.source_edge + 4 > result.edge ||
        result.threads < 1 || !std::isfinite(result.hours) || result.hours <= 0.0 ||
        !std::isfinite(result.step) || result.step <= 0.0) {
        throw std::invalid_argument("invalid vascular benchmark dimensions or time");
    }
    return result;
}

std::uint64_t peak_resident_bytes() {
    struct rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) throw std::runtime_error("unable to read peak memory");
#ifdef __APPLE__
    return static_cast<std::uint64_t>(usage.ru_maxrss);
#else
    return static_cast<std::uint64_t>(usage.ru_maxrss) * 1024ULL;
#endif
}
}  // namespace

int main(int argc, char** argv) {
    try {
        const auto run = arguments(argc, argv);
        atcg3d::continuum::AngiogenesisFieldConfig3D config;
        config.model = "shared_vegf_lattice_v2";
        config.seed_tips_per_hour = 5.0;
        config.tip_branching_per_hour = 0.05;
        config.tip_anastomosis_per_hour = 0.05;
        config.vessel_radius_voxels = 0.5;
        config.validate();
        const std::size_t size = static_cast<std::size_t>(run.edge) * run.edge;
        std::vector<double> cells(size, 0.0), nutrient(size, 0.05), vessels(size, 0.0);
        const int lower = (run.edge - run.source_edge) / 2;
        const int upper = lower + run.source_edge - 1;
        for (int y = lower; y <= upper; ++y) {
            for (int x = lower; x <= upper; ++x) {
                cells[static_cast<std::size_t>(y) * run.edge + x] = 0.5;
            }
        }
        const atcg3d::continuum::VascularConsumerBounds3D range{
            {lower, lower, 0}, {upper, upper, 0}, true};
        atcg3d::continuum::SharedAngiogenesis3D field(
            config, {run.edge, run.edge, 1}, 1.0, true, run.threads);
        field.initialize(std::move(vessels));
        const auto started = std::chrono::steady_clock::now();
        double elapsed = 0.0;
        std::uint64_t updates = 0;
        while (elapsed < run.hours) {
            const double dt = std::min(run.step, run.hours - elapsed);
            field.advance(dt, cells, nutrient, 1.0, true, &range);
            elapsed += dt;
            ++updates;
        }
        const auto seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
        const auto& diagnostic = field.diagnostics();
        std::cout << std::setprecision(17)
                  << "{\"grid_shape\":[" << run.edge << ',' << run.edge << ",1],\"threads\":" << run.threads
                  << ",\"hours\":" << elapsed << ",\"step_hours\":" << run.step << ",\"updates\":" << updates
                  << ",\"source_edge\":" << run.source_edge << ",\"wall_seconds\":" << seconds
                  << ",\"peak_resident_bytes\":" << peak_resident_bytes()
                  << ",\"allocated_field_bytes\":" << field.allocated_bytes()
                  << ",\"active_voxels\":" << diagnostic.active_voxels
                  << ",\"vascular_length\":" << diagnostic.centerline_growth
                  << ",\"vascular_branches\":" << diagnostic.branches
                  << ",\"vascular_anastomoses\":" << diagnostic.anastomoses
                  << ",\"state_checksum\":" << field.checksum() << "}\n";
    } catch (const std::exception& error) {
        std::cerr << "vascular benchmark: " << error.what() << '\n';
        return 1;
    }
}
