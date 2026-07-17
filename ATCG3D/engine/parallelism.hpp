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
                        int available_threads);

int available_worker_threads() noexcept;

}  // namespace atcg3d
