#include <cassert>

#include "config/model_config.hpp"
#include "engine/simulation.hpp"

int main() {
    using namespace atcg3d;
    Model3DConfig config;
    config.output_enabled = false;
    config.initial_r_cells = 6;
    config.initial_K_cells = 6;
    config.initial_radius = 8;
    config.end_time_hours = 8.0;
    config.max_events = 100000;
    config.density_block_edge = 2;
    config.threads = 1;

    Simulation3D first(config);
    first.run();
    assert(first.cells().alive_count() > 0);
    const auto checksum = first.state_checksum();

    config.threads = 4;
    Simulation3D second(config);
    second.run();
    assert(second.state_checksum() == checksum);
    assert(second.cells().alive_count() == first.cells().alive_count());
}
