#pragma once

#include <cstdint>
#include <vector>

#include "config/model_config.hpp"
#include "core/cell_store.hpp"
#include "space/chunk_grid.hpp"
#include "space/density_index.hpp"

namespace atcg3d {

class VascularInfluenceField3D;

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
    bool stage_recovery{};
    Slot daughter{kEmptySlot};
    std::vector<Vec3i> changed_sites;
};

struct GrowthRefreshResult {
    bool division_time_changed{};
    bool death_time_changed{};
    bool migration_activation_changed{};
};

// Deterministic same-time conflict ordering. The hash domain encodes the
// biological event kind; the draw dimension encodes either the configured
// time bucket or the exact IEEE event time when bucketing is disabled.
std::uint64_t division_conflict_priority(const Model3DConfig& config,
                                         double event_time,
                                         CellUid uid);
std::uint64_t stage_recovery_conflict_priority(const Model3DConfig& config,
                                               double event_time,
                                               CellUid uid);

enum class StageRecoveryAction : std::uint8_t {
    none = 0,
    separate_colocated = 1,
    restore_large = 2,
};

struct StageRecoveryProposal {
    Slot slot{kEmptySlot};
    CellUid uid{};
    std::uint64_t event_sequence{};
    Vec3i from{};
    Vec3i target{};
    StageRecoveryAction action{StageRecoveryAction::none};
    std::vector<Vec3i> reserved_sites;
    bool locks_colocation_group{};
};

enum class DivisionAction : std::uint8_t {
    none = 0,
    large_daughter = 1,
    shape_reduction = 2,
    stage_recovery = 3,
    small_daughter = 4,
    remove_r_mother = 5,
    colocated_daughter = 6,
};

struct DivisionProposal {
    Slot mother{kEmptySlot};
    CellUid mother_uid{};
    std::uint64_t event_sequence{};
    Vec3i mother_anchor{};
    CellStage mother_stage{CellStage::small};
    DivisionAction action{DivisionAction::none};
    Vec3i mother_target{};
    Vec3i daughter_target{};
    std::vector<Vec3i> reserved_sites;
    StageRecoveryProposal stage_recovery;
    bool locks_colocation_group{};
};

double sample_division_delay(double density_growth_rate,
                             std::uint64_t seed,
                             CellUid uid,
                             std::uint64_t event_sequence);
double sample_death_delay(double mean_hours,
                          std::uint64_t seed,
                          CellUid uid,
                          std::uint64_t event_sequence);
double sample_division_delay(double density_growth_rate,
                             const DivisionTimingConfig& timing,
                             std::uint64_t seed,
                             CellUid uid,
                             std::uint64_t event_sequence);

GrowthRefreshResult refresh_growth_state(
    Slot slot,
    double now,
    CellStore3D& cells,
    const BlockDensityIndex3D& density,
    const Model3DConfig& config,
    const VascularInfluenceField3D* vascular_influence = nullptr,
    bool refresh_migration_activation = true);

// Draw exactly one new cell-cycle work requirement. Call only for an initial
// cell or for mother/daughter cells after a successful biological division.
void initialize_division_cycle(Slot slot,
                               double now,
                               CellStore3D& cells,
                               const Model3DConfig& config);

// Re-evaluate the configured migration-density gate without touching growth
// work. Returns true only when the active flag crosses the threshold.
bool refresh_migration_activation_state(Slot slot,
                                        double now,
                                        CellStore3D& cells,
                                        const BlockDensityIndex3D& density,
                                        const Model3DConfig& config);

// Activate a cell after an exact external density index has already established
// that its stage-specific threshold is met. This preserves the same eligibility,
// RNG, and timing rules without repeating the 70^3 density query.
bool activate_migration_state_if_density_high(Slot slot,
                                              double now,
                                              CellStore3D& cells,
                                              const Model3DConfig& config);

// Transition an active cell back to ordinary migration at its exact scheduled
// end time. A fresh event-keyed r-cell normal rate is drawn on every expiry.
bool expire_migration_activation_state(Slot slot,
                                       double now,
                                       CellStore3D& cells,
                                       const Model3DConfig& config);

double effective_migration_rate(Slot slot,
                                const CellStore3D& cells,
                                const Model3DConfig& config);

float sample_normal_migration_rate(CellType type,
                                   float inherent_rate,
                                   const Model3DConfig& config,
                                   CellUid uid,
                                   std::uint64_t event_sequence);

bool migration_allowed_for_cell(Slot slot,
                                const CellStore3D& cells,
                                const Model3DConfig& config);

// 3D mapping of the legacy stage-aware 70x70 r-to-K density query. The
// production path uses the incremental block index and never scans CellStore.
double r_to_K_division_density(const BlockDensityIndex3D& density,
                               Vec3i anchor,
                               CellStage stage,
                               const RToKConversionConfig& config,
                               bool thin_layer = false);

// Pure deterministic conversion rule. It does not consume mutable simulation
// RNG state; the draw is keyed by seed, mother UID, and division event sequence.
bool should_convert_r_daughter(CellType mother_type,
                               double local_density,
                               const RToKConversionConfig& config,
                               std::uint64_t seed,
                               CellUid mother_uid,
                               std::uint64_t division_event_sequence);

bool remove_cell(Slot slot,
                 CellStore3D& cells,
                 SparseChunkGrid3D& grid,
                 BlockDensityIndex3D& density);

bool try_stage_recovery(Slot slot,
                        CellStore3D& cells,
                        SparseChunkGrid3D& grid,
                        const Model3DConfig& config,
                        std::uint64_t event_sequence);

StageRecoveryProposal make_stage_recovery_proposal(
    Slot slot,
    const CellStore3D& cells,
    const SparseChunkGrid3D& grid,
    const Model3DConfig& config,
    std::uint64_t event_sequence);
bool commit_stage_recovery_proposal(const StageRecoveryProposal& proposal,
                                    CellStore3D& cells,
                                    SparseChunkGrid3D& grid);

DivisionProposal make_division_proposal(Slot mother,
                                        const CellStore3D& cells,
                                        const SparseChunkGrid3D& grid,
                                        const Model3DConfig& config);
DivisionResult commit_division_proposal(
    const DivisionProposal& proposal,
    double now,
    CellUid& next_uid,
    CellStore3D& cells,
    SparseChunkGrid3D& grid,
    BlockDensityIndex3D& density,
    const Model3DConfig& config,
    std::vector<LineageEdge>& lineage,
    const VascularInfluenceField3D* vascular_influence = nullptr,
    bool refresh_migration_activation = true);

DivisionResult divide_cell(Slot mother,
                           double now,
                           CellUid& next_uid,
                           CellStore3D& cells,
                           SparseChunkGrid3D& grid,
                           BlockDensityIndex3D& density,
                           const Model3DConfig& config,
                           std::vector<LineageEdge>& lineage,
                           const VascularInfluenceField3D* vascular_influence = nullptr);

}  // namespace atcg3d
