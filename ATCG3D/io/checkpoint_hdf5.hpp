#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <span>
#include <vector>

#include "config/model_config.hpp"
#include "engine/simulation.hpp"

namespace atcg3d {

inline constexpr std::uint32_t kCheckpointBaseSchemaVersion3D = 6;
inline constexpr std::uint32_t kCheckpointFieldDeltaSchemaVersion3D = 7;
inline constexpr std::uint32_t kCheckpointJournalDeltaSchemaVersion3D = 8;

enum class CheckpointCellField3D : std::uint64_t {
    anchor = 1ULL << 0U,
    parent_uid = 1ULL << 1U,
    clone_id = 1ULL << 2U,
    type = 1ULL << 3U,
    stage = 1ULL << 4U,
    viability = 1ULL << 5U,
    flags = 1ULL << 6U,
    last_direction = 1ULL << 7U,
    inherent_growth_rate = 1ULL << 8U,
    density_growth_rate = 1ULL << 9U,
    migration_rate = 1ULL << 10U,
    normal_migration_rate = 1ULL << 11U,
    migration_activation_end_time = 1ULL << 12U,
    division_work_remaining = 1ULL << 13U,
    next_migration_time = 1ULL << 14U,
    next_division_time = 1ULL << 15U,
    death_deadline = 1ULL << 16U,
    last_update_time = 1ULL << 17U,
    swap_ready_time = 1ULL << 18U,
    swap_wait_state = 1ULL << 19U,
    pending_swap_direction = 1ULL << 20U,
    event_sequence = 1ULL << 21U,
    migration_schedule_generation = 1ULL << 22U,
    division_schedule_generation = 1ULL << 23U,
    death_schedule_generation = 1ULL << 24U,
};

constexpr std::uint64_t checkpoint_cell_field_mask(
    CheckpointCellField3D field) noexcept {
    return static_cast<std::uint64_t>(field);
}

struct CheckpointData3D {
    std::vector<CellInit> cells;
    std::vector<Slot> cell_slots;
    std::size_t cell_slot_count{};
    std::vector<Slot> cell_free_slots;
    CellUid next_uid{};
    SimulationClock3D clock;
    SimulationStats3D stats;
    std::vector<LineageEdge> lineage;
    VasculatureState3D vasculature;
    std::uint64_t state_checksum{};
};

struct CheckpointSnapshotView3D {
    std::span<const CellInit> cells;
    std::span<const Slot> cell_slots;
    std::size_t cell_slot_count{};
    std::span<const Slot> cell_free_slots;
    CellUid next_uid{};
    SimulationClock3D clock;
    SimulationStats3D stats;
    std::span<const LineageEdge> lineage;
    const VasculatureState3D& vasculature;
    std::uint64_t state_checksum{};
};

struct CheckpointDeltaSummary3D {
    std::size_t changed_cells{};
    std::size_t removed_cells{};
    std::size_t appended_lineage_edges{};
    double changed_fraction{};
};

struct CheckpointJournalSnapshotView3D {
    std::span<const CheckpointCellMutation3D> cell_mutations;
    std::span<const FreeListMutation3D> free_list_mutations;
    std::size_t cell_slot_count{};
    CellUid next_uid{};
    SimulationClock3D clock;
    SimulationStats3D stats;
    std::span<const LineageEdge> lineage_tail;
    std::size_t lineage_prefix_count{};
    const VasculatureState3D& vasculature;
    std::uint64_t state_checksum{};
};

CheckpointDeltaSummary3D summarize_checkpoint_delta(
    const CheckpointSnapshotView3D& current,
    const CheckpointSnapshotView3D& parent);

void write_hdf5_checkpoint(const std::filesystem::path& path,
                           const Simulation3D& simulation);
void write_hdf5_checkpoint(const std::filesystem::path& path,
                           const CheckpointSnapshotView3D& snapshot,
                           const Model3DConfig& config);
void write_hdf5_delta_checkpoint(
    const std::filesystem::path& path,
    const CheckpointSnapshotView3D& snapshot,
    const CheckpointSnapshotView3D& parent_snapshot,
    const std::filesystem::path& parent_path,
    std::uint64_t chain_length,
    const Model3DConfig& config);
void write_hdf5_journal_delta_checkpoint(
    const std::filesystem::path& path,
    const CheckpointJournalSnapshotView3D& snapshot,
    const std::filesystem::path& parent_path,
    std::uint64_t parent_state_checksum,
    double parent_time_hours,
    std::uint64_t chain_length,
    const Model3DConfig& config);
CheckpointData3D read_hdf5_checkpoint(const std::filesystem::path& path,
                                      const Model3DConfig& expected_config);

}  // namespace atcg3d
