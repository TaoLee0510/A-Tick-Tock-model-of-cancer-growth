#pragma once

#include <cstddef>
#include <compare>
#include <cstdint>
#include <optional>
#include <span>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "core/types.hpp"
#include "vasculature/types.hpp"

namespace atcg3d {

struct ExposedFace3D;
class CellStore3D;
class SparseChunkGrid3D;

enum class LesionConnectivity3D : std::uint8_t {
    face_6 = 6,
    full_26 = 26,
};

// Lesion topology is deliberately independent from Model3DConfig. The YAML
// layer can construct this typed, validated policy once during startup without
// introducing string lookup or virtual dispatch in the simulation hot path.
struct LesionIndexConfig3D {
    int block_edge{8};
    LesionConnectivity3D connectivity{LesionConnectivity3D::full_26};
    double core_activation_occupied_fraction{0.15};
    double core_deactivation_occupied_fraction{0.10};
    std::uint64_t minimum_cells_per_core_block{1};
    double minimum_biological_volume_per_core_block{0.0};

    // Non-core occupied blocks can be attributed to a nearby core component
    // for boundary lookup and lesion statistics. They never participate in
    // connectivity, so a sparse migration trail cannot bridge two lesions.
    int halo_blocks{1};
};

struct Vec3d {
    double x{};
    double y{};
    double z{};

    constexpr auto operator<=>(const Vec3d&) const = default;
};

struct LesionBlockStats3D {
    std::uint64_t cell_count{};
    std::uint64_t occupied_voxel_count{};
    double biological_volume{};
    bool core{};
    LesionId lesion_id{kNoLesionId};
};

struct LesionCoreIdentity3D {
    Vec3i block{};
    LesionId lesion_id{kNoLesionId};

    constexpr auto operator<=>(const LesionCoreIdentity3D&) const = default;
};

// Exact pre-refresh observation for one block that was dirty when a checkpoint
// was taken. `exists == false` records that the old observation map had no such
// block even though the current CellStore/Grid may already contain new data.
// Coordinate sums are fixed-width integers so checkpoint serialization does
// not round biologically relevant centroid state through floating point.
struct LesionDirtyBlockState3D {
    Vec3i block{};
    bool exists{};
    std::uint64_t cell_count{};
    std::uint64_t occupied_voxel_count{};
    double biological_volume{};
    std::int64_t cell_coordinate_sum_x{};
    std::int64_t cell_coordinate_sum_y{};
    std::int64_t cell_coordinate_sum_z{};
    std::int64_t occupied_coordinate_sum_x{};
    std::int64_t occupied_coordinate_sum_y{};
    std::int64_t occupied_coordinate_sum_z{};

    constexpr bool operator==(const LesionDirtyBlockState3D&) const = default;
};

struct LesionMerge3D {
    LesionId result{kNoLesionId};
    // All old lesion IDs with core-block overlap in the merged result,
    // including the retained result ID when applicable.
    std::vector<LesionId> predecessors;
};

struct LesionSplit3D {
    LesionId predecessor{kNoLesionId};
    LesionId retained_child{kNoLesionId};
    // Every resulting child ID, including retained_child.
    std::vector<LesionId> children;
};

struct LesionTopologyDelta3D {
    std::vector<LesionId> created;
    std::vector<LesionId> removed;
    std::vector<LesionMerge3D> merges;
    std::vector<LesionSplit3D> splits;

    bool empty() const noexcept {
        return created.empty() && removed.empty() && merges.empty() && splits.empty();
    }
};

struct LesionBiologicalVolumes3D {
    double stage0_large{8.0};
    double stage1_small{1.0};
    double stage2_ultrasmall{1.0};
};

struct LesionSummary3D {
    LesionId id{kNoLesionId};
    std::uint64_t cell_count{};
    std::uint64_t occupied_voxel_count{};
    double biological_volume{};
    // Geometric centroid of unique occupied lattice sites; falls back to the
    // unique-cell anchor centroid only when no occupancy was supplied.
    Vec3d centroid{};

    // Bounds are conservative lattice-site bounds of all assigned coarse
    // blocks (core plus halo). They do not require a per-voxel hash table.
    Vec3i minimum_site{};
    Vec3i maximum_site{};
    std::size_t assigned_block_count{};
    std::vector<Vec3i> core_blocks;
};

// Sparse coarse-block lesion index.
//
// A caller records each living biological cell exactly once with
// add_cell_anchor(), including colocated stage-2 cells. It records each
// occupied lattice site exactly once with add_occupied_site(), so a stage-0
// cell contributes eight occupied sites but still only one cell. This split
// avoids duplicating CellStore3D or SparseChunkGrid3D inside the index.
//
// begin_rebuild() clears observations while retaining the previous core-block
// identity map. refresh_topology() then applies hysteresis and deterministic
// maximum-overlap ID matching. Local simulation integration may instead use
// the add/remove methods between topology refreshes.
class LesionIndex3D {
public:
    explicit LesionIndex3D(LesionIndexConfig3D config);

    const LesionIndexConfig3D& config() const noexcept { return config_; }
    int block_edge() const noexcept { return config_.block_edge; }

    // Clears observations and identity history, including the ID allocator.
    void reset() noexcept;
    // Clears observations but preserves prior topology for hysteresis and
    // stable lesion-ID continuation across a full sparse rebuild.
    void begin_rebuild();

    void add_cell_anchor(Vec3i anchor, double biological_volume = 1.0);
    void remove_cell_anchor(Vec3i anchor, double biological_volume = 1.0);
    void add_occupied_site(Vec3i site);
    void remove_occupied_site(Vec3i site);

    // Convenience full rebuild used by initialization, checkpoint restore and
    // infrequent dirty-topology refreshes. Every valid CellStore slot
    // contributes one cell; stage-specific biological volume is configurable;
    // grid occupancy is visited once per unique occupied lattice site.
    void rebuild_from(const CellStore3D& cells,
                      const SparseChunkGrid3D& grid,
                      LesionBiologicalVolumes3D biological_volumes = {});

    // Hot-path integration: occupancy-changing events only mark their complete
    // old/new footprint sites. Dirty coordinates are deduplicated at block
    // granularity. A later refresh scans block_edge^3 sites and derives unique
    // biological cells through grid occupants, without a CellStore-wide scan.
    std::size_t mark_dirty_sites(std::span<const Vec3i> changed_sites);
    bool observations_dirty() const noexcept { return !dirty_blocks_.empty(); }
    std::size_t dirty_block_count() const noexcept { return dirty_blocks_.size(); }
    std::size_t refresh_dirty_blocks_from(
        const CellStore3D& cells,
        const SparseChunkGrid3D& grid,
        LesionBiologicalVolumes3D biological_volumes = {});

    // Recomputes sparse core components and halo attribution. Complexity is
    // O(number of observed blocks), not O(domain volume). The engine should
    // call this only at a configured lesion refresh/seed event or when
    // topology_dirty() reports a core-threshold crossing, never per cell event.
    LesionTopologyDelta3D refresh_topology();

    bool topology_dirty() const noexcept { return topology_dirty_; }
    bool statistics_dirty() const noexcept { return statistics_dirty_; }
    bool refresh_needed() const noexcept {
        return observations_dirty() || topology_dirty_ || statistics_dirty_;
    }

    Vec3i block_coordinate(Vec3i site) const noexcept;
    std::optional<LesionBlockStats3D> block_stats(Vec3i coordinate) const;
    bool is_core_block(Vec3i coordinate) const noexcept;

    std::span<const LesionSummary3D> lesions() const noexcept { return lesions_; }
    const LesionSummary3D* find_lesion(LesionId id) const noexcept;
    const LesionSummary3D* lesion(LesionId id) const noexcept { return find_lesion(id); }

    // These are coarse attribution queries. The site/face must already be
    // known occupied/exposed by SparseChunkGrid3D/TumorSurfaceIndex3D.
    std::optional<LesionId> lesion_at_site(Vec3i occupied_site) const noexcept;
    std::optional<LesionId> lesion_for_site(Vec3i occupied_site) const noexcept {
        return lesion_at_site(occupied_site);
    }
    std::optional<LesionId> lesion_for_anchor(Vec3i cell_anchor) const noexcept {
        return lesion_at_site(cell_anchor);
    }
    std::optional<LesionId> lesion_for_face(const ExposedFace3D& face) const noexcept;
    // Resolves a stage-0 surface voxel through its owning cell anchor, avoiding
    // ambiguity when a 2x2x2 footprint crosses coarse lesion blocks.
    std::optional<LesionId> lesion_for_face(
        const ExposedFace3D& face,
        const CellStore3D& cells,
        const SparseChunkGrid3D& grid) const noexcept;
    bool lesion_owns_site(LesionId id, Vec3i occupied_site) const noexcept;
    bool lesion_owns_face(LesionId id, const ExposedFace3D& face) const noexcept;

    LesionId next_lesion_id() const noexcept { return next_lesion_id_; }
    // Checkpoint restoration hook. The value must be greater than every ID
    // currently known to the index and must never be the invalid zero ID.
    void set_next_lesion_id(LesionId value);

    // Versioned-checkpoint hooks. Identity snapshots are lexicographically
    // sorted and contain only core blocks; derived halo attribution is rebuilt.
    std::vector<LesionCoreIdentity3D> snapshot_core_identity() const;
    void restore_core_identity(std::span<const LesionCoreIdentity3D> identity,
                               LesionId next_id);
    std::vector<LesionDirtyBlockState3D> snapshot_dirty_block_state() const;
    // Called after rebuilding current observations from CellStore/Grid. Saved
    // old observations are overlaid only for blocks that were dirty at the
    // checkpoint, the old identity/halo/summary are reconstructed, and those
    // coordinates are then marked dirty again. Thus resume does not observe
    // current cell changes until the originally scheduled lesion refresh.
    void restore_checkpoint_state(
        std::span<const LesionCoreIdentity3D> core_identity,
        LesionId next_id,
        std::span<const LesionDirtyBlockState3D> dirty_states);

    std::size_t observed_block_count() const noexcept { return blocks_.size(); }
    std::size_t core_block_count() const noexcept { return core_lesion_by_block_.size(); }
    std::size_t allocated_bytes() const noexcept;

private:
    struct BlockAggregate {
        std::uint64_t cell_count{};
        std::uint64_t occupied_voxel_count{};
        double biological_volume{};
        std::int64_t cell_sum_x{};
        std::int64_t cell_sum_y{};
        std::int64_t cell_sum_z{};
        std::int64_t occupied_sum_x{};
        std::int64_t occupied_sum_y{};
        std::int64_t occupied_sum_z{};
    };

    using BlockMap = std::unordered_map<Vec3i, BlockAggregate, Vec3iHash>;
    using IdentityMap = std::unordered_map<Vec3i, LesionId, Vec3iHash>;

    static int floor_div(int value, int divisor) noexcept;
    bool qualifies_as_core(Vec3i coordinate, const BlockAggregate& block) const noexcept;
    void erase_block_if_empty(Vec3i coordinate);
    LesionId allocate_lesion_id();

    LesionIndexConfig3D config_;
    std::uint64_t block_voxel_capacity_{};
    BlockMap blocks_;

    // previous_core_identity_ is retained by begin_rebuild() and is the sole
    // source for both threshold hysteresis and overlap-based stable matching.
    IdentityMap previous_core_identity_;
    IdentityMap core_lesion_by_block_;
    IdentityMap assigned_lesion_by_block_;
    std::vector<LesionSummary3D> lesions_;
    std::unordered_set<Vec3i, Vec3iHash> dirty_blocks_;
    LesionId next_lesion_id_{1};
    bool topology_dirty_{};
    bool statistics_dirty_{};
};

}  // namespace atcg3d
