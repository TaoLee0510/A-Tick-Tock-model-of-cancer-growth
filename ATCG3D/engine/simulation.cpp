#include "engine/simulation.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <stdexcept>
#include <unordered_set>

#include "core/stateless_rng.hpp"
#include "rules/density.hpp"
#include "rules/initialization.hpp"
#include "rules/migration.hpp"

namespace atcg3d {
namespace {

bool same_time(double lhs, double rhs) {
    return std::abs(lhs - rhs) <= 1e-10 * std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

std::uint64_t hash_combine(std::uint64_t current, std::uint64_t value) {
    return splitmix64(current ^ splitmix64(value + 0x9e3779b97f4a7c15ULL));
}

}  // namespace

Simulation3D::Simulation3D(Model3DConfig config)
    : config_(std::move(config)),
      domain_(config_),
      grid_(config_.chunk_edge, domain_),
      density_(config_.density_block_edge) {
    config_.validate();
}

void Simulation3D::initialize() {
    if (initialized_) {
        throw std::logic_error("simulation is already initialized");
    }
    InitializationResult result = initialize_sphere_and_shell(cells_, grid_, density_, config_);
    next_uid_ = result.next_uid;
    lineage_ = std::move(result.lineage);
    for (const Slot slot : cells_.alive_slots()) {
        schedule_cell(slot);
    }
    initialized_ = true;
}

void Simulation3D::restore(const std::vector<CellInit>& restored_cells,
                           CellUid next_uid,
                           SimulationClock3D clock,
                           SimulationStats3D stats,
                           std::vector<LineageEdge> lineage) {
    if (initialized_ || cells_.slot_count() != 0) {
        throw std::logic_error("restore requires a fresh Simulation3D instance");
    }
    cells_.reserve(restored_cells.size());
    for (const CellInit& cell : restored_cells) {
        const Slot slot = cells_.create(cell);
        bool placed = false;
        if (cell.stage == CellStage::large) {
            placed = grid_.place_large(cell.anchor, slot);
        } else if (grid_.empty(cell.anchor)) {
            placed = grid_.place_single(cell.anchor, slot);
        } else if (cell.stage == CellStage::ultrasmall) {
            placed = grid_.add_colocated(cell.anchor, slot);
        }
        if (!placed) {
            throw std::runtime_error("checkpoint contains conflicting 3D footprints");
        }
        density_.add(cell.anchor, cell.type);
    }
    next_uid_ = next_uid;
    clock_ = clock;
    stats_ = stats;
    lineage_ = std::move(lineage);
    initialized_ = true;
    for (const Slot slot : cells_.alive_slots()) {
        schedule_cell(slot);
    }
}

void Simulation3D::run(const std::function<void(const Simulation3D&)>& observer) {
    if (!initialized_) {
        initialize();
    }
    if (observer) {
        observer(*this);
    }
    double last_observer_time = clock_.time_hours;
    while (clock_.completed_events < config_.max_events && step()) {
        if (observer && (!same_time(clock_.time_hours, last_observer_time) || events_.empty())) {
            observer(*this);
            last_observer_time = clock_.time_hours;
        }
    }
}

bool Simulation3D::step() {
    if (!initialized_) {
        initialize();
    }
    while (!events_.empty() && !current(events_.top())) {
        events_.pop();
    }
    if (events_.empty() || events_.top().time > config_.end_time_hours) {
        clock_.time_hours = std::min(config_.end_time_hours,
                                     events_.empty() ? config_.end_time_hours : events_.top().time);
        return false;
    }

    const double batch_time = events_.top().time;
    std::vector<Event> batch;
    while (!events_.empty() && same_time(events_.top().time, batch_time)) {
        Event event = events_.top();
        events_.pop();
        if (current(event)) {
            batch.push_back(event);
        }
    }
    if (batch.empty()) {
        return !events_.empty();
    }
    clock_.time_hours = batch_time;

    std::vector<Event> migrations;
    for (const Event& event : batch) {
        if (event.kind == EventKind::migration) {
            migrations.push_back(event);
        } else {
            process_non_migration(event);
            ++clock_.completed_events;
        }
    }
    if (!migrations.empty()) {
        process_migrations(migrations);
        clock_.completed_events += migrations.size();
    }
    return true;
}

bool Simulation3D::current(const Event& event) const {
    if (!cells_.valid(event.slot) || cells_.uid(event.slot) != event.uid) {
        return false;
    }
    double expected = 0.0;
    switch (event.kind) {
        case EventKind::migration: expected = cells_.next_migration_time(event.slot); break;
        case EventKind::division: expected = cells_.next_division_time(event.slot); break;
        case EventKind::death: expected = cells_.death_deadline(event.slot); break;
    }
    return event.generation == cells_.schedule_generation(event.slot) &&
           expected > 0.0 && same_time(expected, event.time);
}

void Simulation3D::schedule(EventKind kind, Slot slot, double time, std::uint32_t generation) {
    if (cells_.valid(slot) && time > clock_.time_hours && std::isfinite(time)) {
        events_.push({time, slot, cells_.uid(slot), kind, generation});
    }
}

void Simulation3D::schedule_cell(Slot slot) {
    if (!cells_.valid(slot)) {
        return;
    }
    const std::uint32_t generation = cells_.bump_schedule_generation(slot);
    schedule(EventKind::migration, slot, cells_.next_migration_time(slot), generation);
    schedule(EventKind::division, slot, cells_.next_division_time(slot), generation);
    schedule(EventKind::death, slot, cells_.death_deadline(slot), generation);
}

void Simulation3D::process_non_migration(const Event& event) {
    if (!current(event)) {
        return;
    }
    if (event.kind == EventKind::death) {
        const double rate = density_growth_rate_for_cell(cells_, event.slot, density_, config_);
        if (rate > 0.0) {
            refresh_growth_state(event.slot, clock_.time_hours, cells_, density_, config_);
            schedule_cell(event.slot);
            return;
        }
        const Vec3i site = cells_.anchor(event.slot);
        if (remove_cell(event.slot, cells_, grid_, density_)) {
            ++stats_.deaths;
            recover_neighborhood({site});
            refresh_neighborhood({site});
        }
        return;
    }

    DivisionResult result = divide_cell(event.slot, clock_.time_hours, next_uid_, cells_, grid_,
                                        density_, config_, lineage_);
    if (result.changed) {
        if (result.daughter != kEmptySlot) {
            ++stats_.divisions;
            const double rate = cells_.migration_rate(result.daughter);
            cells_.set_next_migration_time(result.daughter,
                rate > 0.0 ? clock_.time_hours + 1.0 / rate : 0.0);
            schedule_cell(result.daughter);
        } else if (result.mother_removed) {
            ++stats_.deaths;
        }
        recover_neighborhood(result.changed_sites);
        refresh_neighborhood(result.changed_sites);
    } else if (cells_.valid(event.slot)) {
        cells_.set_next_division_time(event.slot, clock_.time_hours + 1.0);
    }
    if (cells_.valid(event.slot)) {
        schedule_cell(event.slot);
    }
}

void Simulation3D::process_migrations(const std::vector<Event>& events) {
    struct Pending {
        Event event;
        std::uint64_t sequence{};
    };
    std::vector<Pending> pending;
    pending.reserve(events.size());
    for (const Event& event : events) {
        if (current(event)) {
            pending.push_back({event, cells_.consume_event_sequence(event.slot)});
        }
    }
    std::vector<MoveProposal> proposals(pending.size());
    const std::uint64_t time_bucket = config_.conflict_bucket_hours > 0.0
        ? static_cast<std::uint64_t>(std::floor(clock_.time_hours / config_.conflict_bucket_hours))
        : std::bit_cast<std::uint64_t>(clock_.time_hours);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(config_.threads)
#endif
    for (std::int64_t index = 0; index < static_cast<std::int64_t>(pending.size()); ++index) {
        proposals[static_cast<std::size_t>(index)] = make_move_proposal(
            pending[static_cast<std::size_t>(index)].event.slot, cells_, grid_, density_, config_,
            pending[static_cast<std::size_t>(index)].sequence, time_bucket);
    }
    std::sort(proposals.begin(), proposals.end(), [](const MoveProposal& lhs, const MoveProposal& rhs) {
        if (lhs.priority != rhs.priority) return lhs.priority > rhs.priority;
        return lhs.uid < rhs.uid;
    });

    std::unordered_set<Vec3i, Vec3iHash> reserved;
    std::vector<Vec3i> changed_sites;
    for (const MoveProposal& proposal : proposals) {
        ++stats_.migration_attempts;
        if (proposal.direction == kStayDirection) {
            if (cells_.valid(proposal.slot)) {
                cells_.set_last_direction(proposal.slot, kStayDirection);
            }
        } else {
            const bool conflict = std::any_of(proposal.reserved_sites.begin(), proposal.reserved_sites.end(),
                                              [&reserved](Vec3i site) { return reserved.contains(site); });
            if (conflict) {
                ++stats_.conflict_rejections;
            } else if (commit_move(proposal, cells_, grid_, density_)) {
                reserved.insert(proposal.reserved_sites.begin(), proposal.reserved_sites.end());
                changed_sites.push_back(proposal.from);
                changed_sites.push_back(proposal.to);
                ++stats_.migration_commits;
            }
        }
        if (cells_.valid(proposal.slot)) {
            const double rate = cells_.migration_rate(proposal.slot);
            cells_.set_next_migration_time(proposal.slot,
                rate > 0.0 ? clock_.time_hours + 1.0 / rate : 0.0);
            schedule_cell(proposal.slot);
        }
    }
    if (!changed_sites.empty()) {
        recover_neighborhood(changed_sites);
        refresh_neighborhood(changed_sites);
    }
}

std::vector<Slot> Simulation3D::nearby_slots(const std::vector<Vec3i>& sites, int radius) const {
    std::unordered_set<Slot> unique;
    for (const Vec3i center : sites) {
        for (int dx = -radius; dx <= radius; ++dx) {
            for (int dy = -radius; dy <= radius; ++dy) {
                for (int dz = config_.thin_layer ? 0 : -radius;
                     dz <= (config_.thin_layer ? 0 : radius); ++dz) {
                    for (const Slot slot : grid_.occupants(center + Vec3i{dx, dy, dz})) {
                        if (cells_.valid(slot)) {
                            unique.insert(slot);
                        }
                    }
                }
            }
        }
    }
    return {unique.begin(), unique.end()};
}

void Simulation3D::recover_neighborhood(const std::vector<Vec3i>& changed_sites) {
    std::vector<Slot> candidates = nearby_slots(changed_sites, 2);
    std::sort(candidates.begin(), candidates.end(), [this](Slot lhs, Slot rhs) {
        return cells_.uid(lhs) < cells_.uid(rhs);
    });
    for (const Slot slot : candidates) {
        if (!cells_.valid(slot) || cells_.stage(slot) == CellStage::large) {
            continue;
        }
        const Vec3i before = cells_.anchor(slot);
        const std::uint64_t sequence = cells_.consume_event_sequence(slot);
        if (try_stage_recovery(slot, cells_, grid_, config_, sequence)) {
            const Vec3i after = cells_.anchor(slot);
            density_.move(before, after, cells_.type(slot));
        }
    }
}

void Simulation3D::refresh_neighborhood(const std::vector<Vec3i>& changed_sites) {
    std::vector<Slot> slots = nearby_slots(changed_sites, config_.growth_density_window_edge / 2 + 1);
    std::sort(slots.begin(), slots.end(), [this](Slot lhs, Slot rhs) {
        return cells_.uid(lhs) < cells_.uid(rhs);
    });
    for (const Slot slot : slots) {
        refresh_growth_state(slot, clock_.time_hours, cells_, density_, config_);
        schedule_cell(slot);
    }
}

std::uint64_t Simulation3D::state_checksum() const {
    std::vector<Slot> slots = cells_.alive_slots();
    std::sort(slots.begin(), slots.end(), [this](Slot lhs, Slot rhs) {
        return cells_.uid(lhs) < cells_.uid(rhs);
    });
    std::uint64_t checksum = splitmix64(slots.size());
    for (const Slot slot : slots) {
        const Vec3i anchor = cells_.anchor(slot);
        checksum = hash_combine(checksum, cells_.uid(slot));
        checksum = hash_combine(checksum, static_cast<std::uint32_t>(anchor.x));
        checksum = hash_combine(checksum, static_cast<std::uint32_t>(anchor.y));
        checksum = hash_combine(checksum, static_cast<std::uint32_t>(anchor.z));
        checksum = hash_combine(checksum, static_cast<std::uint8_t>(cells_.type(slot)));
        checksum = hash_combine(checksum, static_cast<std::uint8_t>(cells_.stage(slot)));
        checksum = hash_combine(checksum, cells_.last_direction(slot));
        checksum = hash_combine(checksum, cells_.event_sequence(slot));
    }
    return checksum;
}

std::vector<CellInit> Simulation3D::snapshot_cells() const {
    std::vector<CellInit> result;
    result.reserve(cells_.alive_count());
    for (const Slot slot : cells_.alive_slots()) {
        result.push_back(cells_.snapshot(slot));
    }
    std::sort(result.begin(), result.end(), [](const CellInit& lhs, const CellInit& rhs) {
        return lhs.uid < rhs.uid;
    });
    return result;
}

}  // namespace atcg3d
