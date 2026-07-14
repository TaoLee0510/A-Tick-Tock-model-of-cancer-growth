#pragma once

#include <cstdint>
#include <functional>
#include <queue>
#include <vector>

#include "config/model_config.hpp"
#include "core/cell_store.hpp"
#include "rules/lifecycle.hpp"
#include "space/chunk_grid.hpp"
#include "space/density_index.hpp"

namespace atcg3d {

enum class EventKind : std::uint8_t {
    death = 0,
    division = 1,
    migration = 2,
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
};

class Simulation3D {
public:
    explicit Simulation3D(Model3DConfig config);

    void initialize();
    void restore(const std::vector<CellInit>& cells,
                 CellUid next_uid,
                 SimulationClock3D clock,
                 SimulationStats3D stats,
                 std::vector<LineageEdge> lineage);
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

    std::uint64_t state_checksum() const;
    std::vector<CellInit> snapshot_cells() const;

private:
    struct Event {
        double time{};
        Slot slot{kEmptySlot};
        CellUid uid{};
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

    bool current(const Event& event) const;
    void schedule_cell(Slot slot);
    void schedule(EventKind kind, Slot slot, double time, std::uint32_t generation);
    void process_non_migration(const Event& event);
    void process_migrations(const std::vector<Event>& events);
    void refresh_neighborhood(const std::vector<Vec3i>& changed_sites);
    void recover_neighborhood(const std::vector<Vec3i>& changed_sites);
    std::vector<Slot> nearby_slots(const std::vector<Vec3i>& sites, int radius) const;

    Model3DConfig config_;
    CellStore3D cells_;
    DomainPolicy domain_;
    SparseChunkGrid3D grid_;
    BlockDensityIndex3D density_;
    std::vector<LineageEdge> lineage_;
    CellUid next_uid_{1};
    SimulationClock3D clock_;
    SimulationStats3D stats_;
    std::priority_queue<Event, std::vector<Event>, EventLater> events_;
    bool initialized_{};
};

}  // namespace atcg3d
