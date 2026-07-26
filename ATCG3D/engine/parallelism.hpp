#pragma once

#include <cstddef>

#include "config/model_config.hpp"

namespace atcg3d {

// Returns the number of workers for a read-only proposal batch. Selection is
// deterministic and depends only on effective configuration and workload;
// it never consumes simulation RNG.
int select_worker_count(std::size_t live_cells,
                        std::size_t work_items,
                        const Model3DConfig& config,
                        int available_threads,
                        std::uint64_t minimum_items_per_thread = 0);

int available_worker_threads() noexcept;

// Avoid entering the OpenMP runtime at all for the overwhelmingly common
// one-worker event/refresh path. An OpenMP region with num_threads(1) still
// allocates a serialized team and showed up prominently in macOS sampling.
template <class Function>
void deterministic_parallel_for(std::size_t count,
                                int workers,
                                Function&& function) {
    if (workers <= 1 || count <= 1) {
        for (std::size_t index = 0; index < count; ++index) function(index);
        return;
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(workers)
    for (std::int64_t index = 0;
         index < static_cast<std::int64_t>(count); ++index) {
        function(static_cast<std::size_t>(index));
    }
#else
    for (std::size_t index = 0; index < count; ++index) function(index);
#endif
}

}  // namespace atcg3d
