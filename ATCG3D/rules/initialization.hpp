#pragma once

#include <vector>

#include "config/model_config.hpp"
#include "core/cell_store.hpp"
#include "rules/lifecycle.hpp"
#include "space/chunk_grid.hpp"
#include "space/density_index.hpp"

namespace atcg3d {

struct InitializationResult {
    CellUid next_uid{1};
    std::vector<LineageEdge> lineage;
};

InitializationResult initialize_sphere_and_shell(CellStore3D& cells,
                                                  SparseChunkGrid3D& grid,
                                                  BlockDensityIndex3D& density,
                                                  const Model3DConfig& config);

}  // namespace atcg3d
