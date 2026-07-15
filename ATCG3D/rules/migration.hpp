#pragma once

#include <cstdint>
#include <vector>

#include "config/model_config.hpp"
#include "core/cell_store.hpp"
#include "space/chunk_grid.hpp"
#include "space/density_index.hpp"

namespace atcg3d {

enum class RngEventKind : std::uint64_t {
    migration_direction = 1,
    division_location = 2,
    division_timing = 3,
    conflict_priority = 4,
    initialization = 5,
    death_timing = 6,
    division_type_conversion = 7,
    division_conflict_priority = 8,
    stage_recovery_conflict_priority = 9,
};

std::vector<DirectionId> feasible_directions(Slot slot,
                                             const CellStore3D& cells,
                                             const SparseChunkGrid3D& grid,
                                             bool thin_layer);

DirectionId select_migration_direction(Slot slot,
                                       const CellStore3D& cells,
                                       const SparseChunkGrid3D& grid,
                                       const BlockDensityIndex3D& density,
                                       const Model3DConfig& config,
                                       std::uint64_t event_sequence);

struct MoveProposal {
    Slot slot{kEmptySlot};
    CellUid uid{};
    DirectionId direction{};
    Vec3i from{};
    Vec3i to{};
    std::vector<Vec3i> reserved_sites;
    std::uint64_t priority{};
};

MoveProposal make_move_proposal(Slot slot,
                                const CellStore3D& cells,
                                const SparseChunkGrid3D& grid,
                                const BlockDensityIndex3D& density,
                                const Model3DConfig& config,
                                std::uint64_t event_sequence,
                                std::uint64_t time_bucket);

bool commit_move(const MoveProposal& proposal,
                 CellStore3D& cells,
                 SparseChunkGrid3D& grid,
                 BlockDensityIndex3D& density);

}  // namespace atcg3d
