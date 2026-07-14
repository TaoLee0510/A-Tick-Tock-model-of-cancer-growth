#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include "config/model_config.hpp"
#include "engine/simulation.hpp"

int main(int argc, char** argv) {
    try {
        std::uint64_t cells = 1000;
        std::uint64_t max_events = 100000;
        int threads = 1;
        for (int index = 1; index < argc; ++index) {
            const std::string argument = argv[index];
            if (argument == "--cells" && index + 1 < argc) cells = std::stoull(argv[++index]);
            else if (argument == "--events" && index + 1 < argc) max_events = std::stoull(argv[++index]);
            else if (argument == "--threads" && index + 1 < argc) threads = std::stoi(argv[++index]);
            else throw std::invalid_argument("unknown or incomplete argument: " + argument);
        }
        if (cells < 2 || max_events == 0 || threads <= 0) {
            throw std::invalid_argument("cells/events/threads are outside valid ranges");
        }
        atcg3d::Model3DConfig config;
        config.output_enabled = false;
        config.initial_r_cells = cells / 2;
        config.initial_K_cells = cells - config.initial_r_cells;
        config.initial_radius = static_cast<int>(std::ceil(std::cbrt(cells * 12.0)));
        config.initial_shell_thickness = std::max(2, config.initial_radius / 4);
        config.end_time_hours = 1000000.0;
        config.max_events = max_events;
        config.threads = threads;
        const auto start = std::chrono::steady_clock::now();
        atcg3d::Simulation3D simulation(config);
        simulation.run();
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        const double rate = seconds > 0.0 ? simulation.clock().completed_events / seconds : 0.0;
        std::cout << std::setprecision(17)
                  << "{\n"
                  << "  \"schema\": \"atcg3d.events.v1\",\n"
                  << "  \"initial_cells\": " << cells << ",\n"
                  << "  \"final_cells\": " << simulation.cells().alive_count() << ",\n"
                  << "  \"completed_events\": " << simulation.clock().completed_events << ",\n"
                  << "  \"simulated_hours\": " << simulation.clock().time_hours << ",\n"
                  << "  \"threads\": " << threads << ",\n"
                  << "  \"wall_seconds\": " << seconds << ",\n"
                  << "  \"events_per_second\": " << rate << ",\n"
                  << "  \"checksum\": " << simulation.state_checksum() << "\n"
                  << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "atcg3d_event_benchmark: " << error.what() << '\n';
        return 1;
    }
}
