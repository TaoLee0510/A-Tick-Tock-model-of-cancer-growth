#include "engine/parallelism.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace atcg3d {

int select_worker_count(std::size_t live_cells,
                        std::size_t work_items,
                        const Model3DConfig& config,
                        int available_threads) {
    if (work_items <= 1 || available_threads <= 1 || config.threads <= 1) return 1;

    const std::size_t int_max =
        static_cast<std::size_t>(std::numeric_limits<int>::max());
    const int work_limit = static_cast<int>(std::min(work_items, int_max));
    const int upper = std::max(
        1, std::min({config.threads, available_threads, work_limit}));
    if (config.parallel_mode == "fixed") return upper;

    double fraction = config.parallel_thread_thresholds.front().max_thread_fraction;
    for (const ParallelThreadThresholdConfig& threshold :
         config.parallel_thread_thresholds) {
        if (live_cells < threshold.minimum_cells) break;
        fraction = threshold.max_thread_fraction;
    }
    const int population_limit = std::max(
        config.parallel_min_threads,
        static_cast<int>(std::floor(static_cast<double>(config.threads) * fraction)));

    const std::size_t minimum_batch =
        static_cast<std::size_t>(config.parallel_min_events_per_thread);
    const std::size_t event_workers = 1 + (work_items - 1) / minimum_batch;
    const int event_limit = static_cast<int>(std::min(event_workers, int_max));
    const int selected = std::max(
        config.parallel_min_threads, std::min(population_limit, event_limit));
    return std::clamp(selected, 1, upper);
}

int available_worker_threads() noexcept {
#ifdef _OPENMP
    return std::max(1, omp_get_max_threads());
#else
    return 1;
#endif
}

}  // namespace atcg3d
