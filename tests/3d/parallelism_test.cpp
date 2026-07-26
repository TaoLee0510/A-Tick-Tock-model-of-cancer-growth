#include <cassert>

#include "config/model_config.hpp"
#include "engine/parallelism.hpp"

int main() {
    using namespace atcg3d;

    Model3DConfig config;
    config.threads = 10;
    config.parallel_mode = "adaptive_cells_and_events_v1";
    config.parallel_min_threads = 1;
    config.parallel_min_events_per_thread = 100;
    config.validate();

    assert(select_worker_count(0, 10000, config, 12) == 5);
    assert(select_worker_count(4999, 10000, config, 12) == 5);
    assert(select_worker_count(5000, 10000, config, 12) == 6);
    assert(select_worker_count(10000, 10000, config, 12) == 7);
    assert(select_worker_count(15000, 10000, config, 12) == 8);
    assert(select_worker_count(20000, 10000, config, 12) == 9);
    assert(select_worker_count(25000, 10000, config, 12) == 10);

    // Event batches independently cap workers, preventing OpenMP startup for
    // tiny batches even when the population is large.
    assert(select_worker_count(25000, 1, config, 12) == 1);
    assert(select_worker_count(25000, 100, config, 12) == 1);
    assert(select_worker_count(25000, 101, config, 12) == 2);
    assert(select_worker_count(25000, 250, config, 12) == 3);
    assert(select_worker_count(25000, 10000, config, 4) == 4);

    config.parallel_mode = "fixed";
    assert(select_worker_count(0, 10000, config, 12) == 10);
    assert(select_worker_count(0, 3, config, 12) == 3);

    Model3DConfig adaptive = config;
    adaptive.parallel_mode = "adaptive_cells_and_events_v1";
    adaptive.parallel_min_events_per_thread = 37;
    assert(adaptive.dynamics_json() == config.dynamics_json());
    adaptive.proposal_window_hours = 7.5;
    adaptive.proposal_window_max_events = 9999;
    adaptive.proposal_dependency_block_edge = 3;
    adaptive.output_async_enabled = true;
    adaptive.output_async_queue_depth = 4;
    assert(adaptive.dynamics_json() == config.dynamics_json());
    return 0;
}
