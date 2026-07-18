#pragma once

#include <cstdint>
#include <functional>
#include <queue>
#include <span>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "config/model_config.hpp"
#include "core/cell_store.hpp"
#include "rules/lifecycle.hpp"
#include "space/chunk_grid.hpp"
#include "space/density_index.hpp"
#include "vasculature/angiogenesis_process.hpp"
#include "vasculature/influence_field.hpp"
#include "vasculature/lesion_index.hpp"
#include "vasculature/surface_index.hpp"
#include "vasculature/vessel_grid.hpp"
#include "vasculature/vessel_store.hpp"

namespace atcg3d {

// Numeric order is the deterministic within-time-bucket biological priority.
enum class EventKind : std::uint8_t {
    death = 0,
    lesion_refresh = 1,
    angiogenesis_seed = 2,
    vessel_growth = 3,
    migration_activation_end = 4,
    division = 5,
    migration = 6,
};

struct SimulationClock3D {
    std::uint64_t completed_events{};
    double time_hours{};
};

struct SimulationStats3D {
    std::uint64_t migration_attempts{};
    std::uint64_t migration_commits{};
    std::uint64_t divisions{};
    std::uint64_t deaths{};
    std::uint64_t conflict_rejections{};
    std::uint64_t angiogenesis_seed_attempts{};
    std::uint64_t angiogenesis_roots{};
    std::uint64_t angiogenesis_seed_rejections{};
    std::uint64_t vessel_growth_attempts{};
    std::uint64_t vessel_growth_commits{};
    std::uint64_t vessel_anastomoses{};
    std::uint64_t vascular_displacements{};
};

// Cumulative locality diagnostics for vascular growth refreshes. These
// counters are derived execution data: they neither affect biology nor belong
// in checkpoints/state checksums.
struct VascularRefreshDiagnostics3D {
    std::uint64_t queried_density_blocks{};
    std::uint64_t visited_density_slots{};
    std::uint64_t refreshed_cell_slots{};
};

struct LesionAngiogenesisState3D {
    LesionId lesion_id{kNoLesionId};
    AngiogenesisProcessState3D process;
};

// Immutable vessel provenance remains source_lesion_id.  This sparse table
// maps a historical source to its current topological owner after merges and
// removals.  current_lesion_id == kNoLesionId means the source no longer has a
// living lesion owner.  Identity mappings are implicit and are not stored.
struct LesionSourceOwnership3D {
    LesionId source_lesion_id{kNoLesionId};
    LesionId current_lesion_id{kNoLesionId};

    constexpr bool operator==(const LesionSourceOwnership3D&) const = default;
};

struct LesionSimulationState3D {
    LesionId next_lesion_id{1};
    double last_refresh_time_hours{};
    double next_refresh_time_hours{};
    std::uint32_t refresh_schedule_generation{};
    std::vector<LesionCoreIdentity3D> core_identity;
    std::vector<LesionDirtyBlockState3D> dirty_blocks;
    std::vector<LesionAngiogenesisState3D> processes;
    std::vector<LesionSourceOwnership3D> source_ownership;
};

struct VasculatureState3D {
    // Aggregate compatibility summary. Scheduling uses lesions.processes;
    // counters are the exact sum of every current or retired lesion record.
    AngiogenesisProcessState3D process;
    LesionSimulationState3D lesions;
    std::vector<VesselNodeInit3D> nodes;
    std::vector<VesselTipInit3D> tips;
    std::vector<VesselId> perfused_vessels;
    VesselId next_vessel_id{1};
    VesselNodeUid next_node_uid{1};
    VesselTipUid next_tip_uid{1};
};

// Builds the compatibility summary in lesion-id order.  Active intervals are
// folded through snapshot_time_hours, so accumulated_eligible_hours is the
// exact sum at the snapshot and an active aggregate starts a fresh interval at
// that same time.  Throws on duplicate IDs, invalid time, or counter overflow.
AngiogenesisProcessState3D aggregate_angiogenesis_process_states(
    std::span<const LesionAngiogenesisState3D> processes,
    double snapshot_time_hours);

class Simulation3D {
public:
    explicit Simulation3D(Model3DConfig config);

    void initialize();
    void restore(const std::vector<CellInit>& cells,
                 CellUid next_uid,
                 SimulationClock3D clock,
                 SimulationStats3D stats,
                 std::vector<LineageEdge> lineage,
                 const VasculatureState3D& vasculature = {},
                 std::size_t cell_slot_count = 0,
                 const std::vector<Slot>& cell_slots = {},
                 const std::vector<Slot>& cell_free_slots = {});
    void run(const std::function<void(const Simulation3D&)>& observer = {});
    bool step();

    const Model3DConfig& config() const noexcept { return config_; }
    const CellStore3D& cells() const noexcept { return cells_; }
    CellStore3D& cells() noexcept { return cells_; }
    const SparseChunkGrid3D& grid() const noexcept { return grid_; }
    const BlockDensityIndex3D& density() const noexcept { return density_; }
    const std::vector<LineageEdge>& lineage() const noexcept { return lineage_; }
    const SimulationClock3D& clock() const noexcept { return clock_; }
    const SimulationStats3D& stats() const noexcept { return stats_; }
    CellUid next_uid() const noexcept { return next_uid_; }
    std::size_t pending_event_count() const noexcept { return events_.size(); }
    std::uint64_t event_queue_rebuild_count() const noexcept {
        return event_queue_rebuild_count_;
    }
    std::uint64_t migration_activation_bulk_slot_visits() const noexcept {
        return migration_activation_bulk_slot_visits_;
    }
    std::uint64_t migration_activation_direct_slot_visits() const noexcept {
        return migration_activation_direct_slot_visits_;
    }
    std::uint64_t migration_activation_class_recomputes() const noexcept {
        return migration_activation_class_recomputes_;
    }
    std::size_t migration_activation_class_cache_size() const noexcept {
        return migration_activation_class_cache_.size();
    }
    const VascularRefreshDiagnostics3D& vascular_refresh_diagnostics() const noexcept {
        return vascular_refresh_diagnostics_;
    }

    const VesselNodeStore3D& vessel_nodes() const noexcept { return vessel_nodes_; }
    const VesselTipStore3D& vessel_tips() const noexcept { return vessel_tips_; }
    const SparseVesselGrid3D& vessel_grid() const noexcept { return vessel_grid_; }
    const VascularInfluenceField3D& vascular_influence() const noexcept {
        return vascular_influence_;
    }
    const TumorSurfaceIndex3D& tumor_surface() const noexcept { return tumor_surface_; }
    const LesionIndex3D& lesion_index() const noexcept { return lesion_index_; }
    AngiogenesisProcessState3D angiogenesis_state() const;
    std::size_t active_vessel_tip_count() const;
    std::size_t active_vessel_tip_count(LesionId source_lesion_id) const;
    double biological_tumor_volume() const noexcept;
    VasculatureState3D snapshot_vasculature() const;

    std::uint64_t state_checksum() const;
    std::vector<CellInit> snapshot_cells() const;
    std::vector<Slot> snapshot_cell_slots() const;

private:
    struct Event {
        double time{};
        std::uint32_t slot{kEmptySlot};
        std::uint64_t uid{};
        EventKind kind{EventKind::migration};
        std::uint32_t generation{};
    };

    struct EventLater {
        bool operator()(const Event& lhs, const Event& rhs) const noexcept {
            if (lhs.time != rhs.time) return lhs.time > rhs.time;
            if (lhs.kind != rhs.kind) return lhs.kind > rhs.kind;
            return lhs.uid > rhs.uid;
        }
    };

    struct VesselGrowthProposal {
        Event event;
        Vec3i from{};
        Vec3i to{};
        DirectionId direction{kStayDirection};
        std::vector<Vec3i> capsule;
        std::uint64_t priority{};
        bool anastomosis{};
        bool valid{};
    };

    bool current(const Event& event) const;
    void schedule_cell(Slot slot);
    void restore_cell_events(Slot slot);
    void reschedule_event(EventKind kind, Slot slot);
    void apply_growth_refresh(Slot slot, const GrowthRefreshResult& refresh);
    void synchronize_migration_schedule(Slot slot, bool activation_changed);
    double event_time(EventKind kind, Slot slot) const;
    std::uint32_t event_generation(EventKind kind, Slot slot) const;
    std::uint32_t bump_event_generation(EventKind kind, Slot slot);
    void schedule(EventKind kind, std::uint32_t slot, std::uint64_t uid,
                  double time, std::uint32_t generation);
    void reset_event_queue_rebuild_threshold();
    void maybe_compact_event_queue();
    void process_non_migration(const Event& event);
    void process_deaths(const std::vector<Event>& events);
    void process_divisions(const std::vector<Event>& events);
    void process_migrations(const std::vector<Event>& events);

    void rebuild_tumor_surface();
    void refresh_tumor_surface(std::span<const Vec3i> changed_sites);
    void rebuild_lesion_index();
    void mark_lesion_dirty(std::span<const Vec3i> changed_sites);
    bool refresh_lesion_index(bool force);
    void schedule_lesion_refresh_event();
    void sync_angiogenesis_eligibility(bool force_refresh = false);
    void recompute_aggregate_angiogenesis_state();
    LesionId current_lesion_for_source(LesionId source_lesion_id) const noexcept;
    void update_lesion_source_ownership(
        const std::unordered_map<LesionId, LesionId>& successors,
        const std::unordered_set<LesionId>& current_lesions);
    void schedule_seed_event(LesionId lesion_id);
    bool process_seed_event(const Event& event);
    bool create_vessel_root(LesionId source_lesion_id,
                            const ExposedFace3D& face,
                            Vec3i inward_target,
                            std::uint64_t seed_event_sequence);
    bool root_has_local_support(LesionId source_lesion_id,
                                std::span<const Vec3i> root_capsule) const;

    std::vector<DirectionId> feasible_vessel_directions(VesselTipSlot slot) const;
    void schedule_vessel_tip(VesselTipSlot slot);
    void restore_vessel_tip_event(VesselTipSlot slot);
    VesselGrowthProposal make_vessel_growth_proposal(const Event& event) const;
    void process_vessel_growth(const std::vector<Event>& events);
    bool commit_vessel_growth(const VesselGrowthProposal& proposal);
    void activate_vessel_perfusion(VesselId vessel_id);
    bool vessel_is_perfused(VesselId vessel_id) const;

    std::vector<Vec3i> displace_cells_for_vessel(std::span<const Vec3i> capsule);
    void refresh_growth_near_vessel(std::span<const Vec3i> vessel_sites);
    void refresh_neighborhood(const std::vector<Vec3i>& changed_sites);
    Vec3i migration_activation_query_block(Vec3i anchor) const noexcept;
    std::uint8_t migration_activation_class(Vec3i query_block) const;
    void rebuild_migration_activation_class_cache();
    void apply_migration_activation_class(Slot slot, std::uint8_t classes);
    void refresh_migration_activation_near(
        const std::vector<Vec3i>& changed_sites);
    void recover_neighborhood(std::vector<Vec3i>& changed_sites);
    std::vector<Slot> nearby_slots(const std::vector<Vec3i>& sites, int radius) const;

    Model3DConfig config_;
    CellStore3D cells_;
    DomainPolicy domain_;
    SparseVesselGrid3D vessel_grid_;
    SparseChunkGrid3D grid_;
    BlockDensityIndex3D density_;
    VesselNodeStore3D vessel_nodes_;
    VesselTipStore3D vessel_tips_;
    VascularInfluenceField3D vascular_influence_;
    TumorSurfaceIndex3D tumor_surface_;
    LesionIndex3D lesion_index_;
    std::unordered_map<LesionId, AngiogenesisProcess3D>
        lesion_angiogenesis_processes_;
    std::unordered_map<LesionId, LesionId> lesion_source_ownership_;
    AngiogenesisProcessState3D aggregate_angiogenesis_state_{};
    double last_lesion_refresh_time_hours_{};
    double next_lesion_refresh_time_hours_{};
    std::uint32_t lesion_refresh_schedule_generation_{};
    std::unordered_map<Vec3i, VesselNodeSlot, Vec3iHash> centerline_nodes_;
    std::unordered_set<VesselId> perfused_vessels_;
    std::vector<LineageEdge> lineage_;
    CellUid next_uid_{1};
    VesselId next_vessel_id_{1};
    VesselNodeUid next_vessel_node_uid_{1};
    VesselTipUid next_vessel_tip_uid_{1};
    SimulationClock3D clock_;
    SimulationStats3D stats_;
    std::priority_queue<Event, std::vector<Event>, EventLater> events_;
    // Derived, deterministic state. Bit 0 is the small/ultrasmall activation
    // class and bit 1 is the large-cell class for a quantized query block.
    // It is rebuilt from CellStore+density after initialization/resume and is
    // intentionally absent from checkpoints and biological checksums.
    std::unordered_map<Vec3i, std::uint8_t, Vec3iHash>
        migration_activation_class_cache_;
    std::uint64_t migration_activation_bulk_slot_visits_{};
    std::uint64_t migration_activation_direct_slot_visits_{};
    std::uint64_t migration_activation_class_recomputes_{};
    VascularRefreshDiagnostics3D vascular_refresh_diagnostics_;
    std::size_t event_queue_rebuild_threshold_{256};
    std::uint64_t event_queue_rebuild_count_{};
    bool initialized_{};
};

}  // namespace atcg3d
