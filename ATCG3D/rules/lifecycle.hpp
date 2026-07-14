#pragma once

#include <cstdint>
#include <vector>

#include "config/model_config.hpp"
#include "core/cell_store.hpp"
#include "space/chunk_grid.hpp"
#include "space/density_index.hpp"

namespace atcg3d {

struct LineageEdge {
    double birth_time{};
    CellUid child_uid{};
    CellUid parent_uid{};
    std::uint32_t clone_id{};
    CellType type{CellType::r};
};

struct DivisionResult {
    bool changed{};
    bool mother_removed{};
    Slot daughter{kEmptySlot};
    std::vector<Vec3i> changed_sites;
};

double sample_division_delay(double density_growth_rate,
                             std::uint64_t seed,
                             CellUid uid,
                             std::uint64_t event_sequence);

void refresh_growth_state(Slot slot,
                          double now,
                          CellStore3D& cells,
                          const BlockDensityIndex3D& density,
                          const Model3DConfig& config);

bool remove_cell(Slot slot,
                 CellStore3D& cells,
                 SparseChunkGrid3D& grid,
                 BlockDensityIndex3D& density);

bool try_stage_recovery(Slot slot,
                        CellStore3D& cells,
                        SparseChunkGrid3D& grid,
                        const Model3DConfig& config,
                        std::uint64_t event_sequence);

DivisionResult divide_cell(Slot mother,
                           double now,
                           CellUid& next_uid,
                           CellStore3D& cells,
                           SparseChunkGrid3D& grid,
                           BlockDensityIndex3D& density,
                           const Model3DConfig& config,
                           std::vector<LineageEdge>& lineage);

}  // namespace atcg3d
