#include "engine/simulation.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

#include "core/stateless_rng.hpp"
#include "engine/parallelism.hpp"
#include "geometry/directions.hpp"
#include "geometry/footprint.hpp"
#include "rules/density.hpp"
#include "rules/initialization.hpp"
#include "rules/migration.hpp"
#include "vasculature/geometry.hpp"
#include "vasculature/growth_rules.hpp"

namespace atcg3d {
namespace {

constexpr std::uint8_t kSmallMigrationActivationClass = 1U << 0U;
constexpr std::uint8_t kLargeMigrationActivationClass = 1U << 1U;

constexpr std::uint64_t kVesselConflictEvent = 0x564553434f4e464cULL;

bool same_time(double lhs, double rhs) {
    return std::abs(lhs - rhs) <= 1e-10 * std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

double next_periodic_output_boundary(const Model3DConfig& config,
                                     double after_time) noexcept {
    if (!config.output_enabled) return std::numeric_limits<double>::infinity();
    double next = std::numeric_limits<double>::infinity();
    const auto consider = [&](double interval) {
        if (!(interval > 0.0)) return;
        const double completed_intervals = std::floor(after_time / interval);
        double candidate = (completed_intervals + 1.0) * interval;
        if (candidate < after_time || same_time(candidate, after_time)) {
            candidate += interval;
        }
        next = std::min(next, candidate);
    };
    consider(config.preview_every_hours);
    consider(config.full_every_hours);
    consider(config.checkpoint_every_hours);
    return next;
}

double float_storage_time_tolerance(double lhs, double rhs) noexcept {
    return 2.0 * static_cast<double>(std::numeric_limits<float>::epsilon()) *
           std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

std::uint64_t hash_combine(std::uint64_t current, std::uint64_t value) {
    return splitmix64(current ^ splitmix64(value + 0x9e3779b97f4a7c15ULL));
}

std::uint64_t float_bits(float value) {
    return std::bit_cast<std::uint32_t>(value);
}

std::uint64_t double_bits(double value) {
    return std::bit_cast<std::uint64_t>(value);
}

std::int64_t squared_distance(Vec3i lhs, Vec3i rhs) {
    const std::int64_t dx = static_cast<std::int64_t>(lhs.x) - rhs.x;
    const std::int64_t dy = static_cast<std::int64_t>(lhs.y) - rhs.y;
    const std::int64_t dz = static_cast<std::int64_t>(lhs.z) - rhs.z;
    return dx * dx + dy * dy + dz * dz;
}

std::int64_t floor_div_positive(std::int64_t value,
                                std::int64_t divisor) noexcept {
    std::int64_t quotient = value / divisor;
    if (value < 0 && value % divisor != 0) --quotient;
    return quotient;
}

bool influenced_by_any_source(Vec3i anchor,
                              const std::vector<Vec3i>& sources,
                              double cutoff_squared) noexcept {
    for (const Vec3i source : sources) {
        const double dx = static_cast<double>(anchor.x) - source.x;
        const double dy = static_cast<double>(anchor.y) - source.y;
        const double dz = static_cast<double>(anchor.z) - source.z;
        if (dx * dx + dy * dy + dz * dz < cutoff_squared) return true;
    }
    return false;
}

std::vector<Vec3i> occupied_sites_for_cell(const CellStore3D& cells, Slot slot) {
    if (!cells.valid(slot)) return {};
    if (cells.stage(slot) == CellStage::large) {
        const auto footprint = large_footprint(cells.anchor(slot));
        return {footprint.begin(), footprint.end()};
    }
    return {cells.anchor(slot)};
}

bool has_specific_generations(const CellInit& cell) {
    return cell.migration_schedule_generation != 0 ||
           cell.division_schedule_generation != 0 ||
           cell.death_schedule_generation != 0;
}

VascularInfluenceProfile3D influence_profile_from_config(
    const std::string& profile) {
    if (profile == "linear_cutoff") {
        return VascularInfluenceProfile3D::linear_cutoff;
    }
    if (profile == "exponential") {
        return VascularInfluenceProfile3D::exponential;
    }
    throw std::invalid_argument("unsupported vascular influence profile: " + profile);
}

LesionIndexConfig3D lesion_index_config(const Model3DConfig& config) {
    LesionIndexConfig3D result;
    result.block_edge = config.angiogenesis.lesion_block_edge;
    result.connectivity = config.angiogenesis.lesion_connectivity == 6
        ? LesionConnectivity3D::face_6
        : LesionConnectivity3D::full_26;
    result.core_activation_occupied_fraction =
        config.angiogenesis.lesion_core_activation_occupied_fraction;
    result.core_deactivation_occupied_fraction =
        config.angiogenesis.lesion_core_deactivation_occupied_fraction;
    result.minimum_cells_per_core_block =
        config.angiogenesis.lesion_minimum_cells_per_core_block;
    result.minimum_biological_volume_per_core_block =
        config.angiogenesis.lesion_minimum_biological_volume_per_core_block;
    result.halo_blocks = config.angiogenesis.lesion_halo_blocks;
    return result;
}

LesionBiologicalVolumes3D lesion_biological_volumes(
    const Model3DConfig& config) noexcept {
    return {
        config.angiogenesis.stage0_biological_volume_voxels3,
        config.angiogenesis.stage1_biological_volume_voxels3,
        config.angiogenesis.stage2_biological_volume_voxels3,
    };
}

void checked_counter_add(std::uint64_t& target,
                         std::uint64_t value,
                         const char* name) {
    if (value > std::numeric_limits<std::uint64_t>::max() - target) {
        throw std::overflow_error(std::string("angiogenesis ") + name +
                                  " overflow during lesion aggregation");
    }
    target += value;
}

AngiogenesisProcessState3D merge_angiogenesis_process_states(
    std::span<const LesionAngiogenesisState3D> states,
    double now_hours) {
    AngiogenesisProcessState3D merged =
        aggregate_angiogenesis_process_states(states, now_hours);
    if (merged.schedule_generation ==
        std::numeric_limits<std::uint32_t>::max()) {
        throw std::overflow_error(
            "angiogenesis schedule generation overflow during lesion merge");
    }
    ++merged.schedule_generation;
    return merged;
}

Vec3i rounded_lesion_centroid(const LesionSummary3D& lesion) noexcept {
    return {static_cast<std::int32_t>(std::llround(lesion.centroid.x)),
            static_cast<std::int32_t>(std::llround(lesion.centroid.y)),
            static_cast<std::int32_t>(std::llround(lesion.centroid.z))};
}

void validate_restored_cell_schedule(const CellInit& cell,
                                     const SimulationClock3D& clock,
                                     const Model3DConfig& config) {
    const double tolerance = float_storage_time_tolerance(
        clock.time_hours, cell.last_update_time);
    if (!std::isfinite(cell.last_update_time) ||
        cell.last_update_time > clock.time_hours + tolerance) {
        throw std::runtime_error(
            "restored cell update time exceeds the simulation clock");
    }
    for (const double time : {cell.next_migration_time,
                              cell.next_division_time,
                              cell.death_deadline}) {
        if (!std::isfinite(time) || time < 0.0 ||
            (time > 0.0 && time < clock.time_hours)) {
            throw std::runtime_error("restored cell contains an invalid event time");
        }
    }

    const bool active =
        (cell.flags & static_cast<std::uint8_t>(kMigrationActive)) != 0;
    if ((!config.migration_activation_enabled && active) ||
        active != (cell.migration_activation_end_time > clock.time_hours) ||
        (!active && cell.migration_activation_end_time != 0.0)) {
        throw std::runtime_error(
            "restored cell migration activation state/end time are inconsistent");
    }
    const double effective_rate = active ? cell.migration_rate
        : (cell.normal_migration_rate > 0.0F
               ? cell.normal_migration_rate : cell.migration_rate);
    const bool expects_migration = effective_rate > 0.0;
    if ((cell.next_migration_time > 0.0) != expects_migration) {
        throw std::runtime_error(
            "restored cell migration activation/rate and event time are inconsistent");
    }

    const bool growth_active =
        cell.density_growth_rate > config.death_growth_rate_threshold;
    if ((cell.next_division_time > 0.0) != growth_active ||
        (cell.death_deadline > 0.0) == growth_active) {
        throw std::runtime_error(
            "restored cell growth, division, and death schedules are inconsistent");
    }
}

}  // namespace

AngiogenesisProcessState3D aggregate_angiogenesis_process_states(
    std::span<const LesionAngiogenesisState3D> processes,
    double snapshot_time_hours) {
    if (!std::isfinite(snapshot_time_hours) || snapshot_time_hours < 0.0) {
        throw std::invalid_argument(
            "angiogenesis aggregate snapshot time must be finite and nonnegative");
    }

    std::vector<const LesionAngiogenesisState3D*> ordered;
    ordered.reserve(processes.size());
    for (const LesionAngiogenesisState3D& entry : processes) {
        ordered.push_back(&entry);
    }
    std::sort(ordered.begin(), ordered.end(), [](const auto* lhs, const auto* rhs) {
        return lhs->lesion_id < rhs->lesion_id;
    });

    AngiogenesisProcessState3D aggregate;
    double next_time = std::numeric_limits<double>::infinity();
    LesionId previous = kNoLesionId;
    bool first = true;
    for (const LesionAngiogenesisState3D* entry : ordered) {
        if (entry->lesion_id == kNoLesionId ||
            (!first && entry->lesion_id == previous)) {
            throw std::invalid_argument(
                "angiogenesis aggregate lesion IDs must be nonzero and unique");
        }
        first = false;
        previous = entry->lesion_id;

        const AngiogenesisProcessState3D& state = entry->process;
        if (!std::isfinite(state.accumulated_eligible_hours) ||
            state.accumulated_eligible_hours < 0.0 ||
            !std::isfinite(state.eligibility_started_hours) ||
            state.eligibility_started_hours < 0.0 ||
            !std::isfinite(state.next_seed_time_hours) ||
            state.next_seed_time_hours < 0.0) {
            throw std::invalid_argument(
                "angiogenesis aggregate contains an invalid process time");
        }

        double elapsed = state.accumulated_eligible_hours;
        if (state.eligible) {
            if (state.eligibility_started_hours > snapshot_time_hours) {
                throw std::invalid_argument(
                    "angiogenesis eligibility starts after aggregate snapshot time");
            }
            elapsed += snapshot_time_hours - state.eligibility_started_hours;
            aggregate.eligible = true;
            if (state.next_seed_time_hours > 0.0) {
                next_time = std::min(next_time, state.next_seed_time_hours);
            }
        }
        aggregate.accumulated_eligible_hours += elapsed;
        if (!std::isfinite(elapsed) ||
            !std::isfinite(aggregate.accumulated_eligible_hours)) {
            throw std::overflow_error(
                "angiogenesis eligible time overflow during lesion aggregation");
        }
        aggregate.event_sequence =
            std::max(aggregate.event_sequence, state.event_sequence);
        aggregate.schedule_generation =
            std::max(aggregate.schedule_generation, state.schedule_generation);
        checked_counter_add(aggregate.attempted_events,
                            state.attempted_events, "attempt counter");
        checked_counter_add(aggregate.committed_roots,
                            state.committed_roots, "root counter");
        checked_counter_add(aggregate.rejected_events,
                            state.rejected_events, "rejection counter");
    }
    aggregate.next_seed_time_hours =
        std::isfinite(next_time) ? next_time : 0.0;
    aggregate.eligibility_started_hours =
        aggregate.eligible ? snapshot_time_hours : 0.0;
    return aggregate;
}

Simulation3D::Simulation3D(Model3DConfig config)
    : config_(std::move(config)),
      domain_(config_),
      vessel_grid_(config_.chunk_edge, domain_),
      grid_(config_.chunk_edge, domain_),
      density_(config_.density_block_edge),
      vascular_influence_(config_.density_block_edge,
                          static_cast<float>(config_.angiogenesis.influence_cutoff_radius_voxels),
                          static_cast<float>(config_.angiogenesis.influence_max_relief_fraction),
                          influence_profile_from_config(
                              config_.angiogenesis.influence_profile),
                          static_cast<float>(
                              config_.angiogenesis.influence_decay_length_voxels)),
      lesion_index_(lesion_index_config(config_)) {
    config_.validate();
    grid_.attach_vessel_grid(&vessel_grid_);
}

void Simulation3D::initialize() {
    if (initialized_) {
        throw std::logic_error("simulation is already initialized");
    }
    InitializationResult result = initialize_sphere_and_shell(cells_, grid_, density_, config_);
    next_uid_ = result.next_uid;
    lineage_ = std::move(result.lineage);
    rebuild_migration_activation_class_cache();
    if (config_.angiogenesis.enabled) {
        rebuild_tumor_surface();
        rebuild_lesion_index();
    }
    for (const Slot slot : cells_.alive_slots()) schedule_cell(slot);
    initialized_ = true;
    sync_angiogenesis_eligibility(true);
    reset_event_queue_rebuild_threshold();
}

void Simulation3D::restore(const std::vector<CellInit>& restored_cells,
                           CellUid next_uid,
                           SimulationClock3D clock,
                           SimulationStats3D stats,
                           std::vector<LineageEdge> lineage,
                           const VasculatureState3D& vasculature,
                           std::size_t cell_slot_count,
                           const std::vector<Slot>& cell_slots,
                           const std::vector<Slot>& cell_free_slots) {
    if (initialized_ || cells_.slot_count() != 0) {
        throw std::logic_error("restore requires a fresh Simulation3D instance");
    }
    const bool has_explicit_layout =
        cell_slot_count != 0 || !cell_slots.empty() || !cell_free_slots.empty();
    if (has_explicit_layout) {
        cells_.restore_layout(cell_slot_count, cell_slots, restored_cells,
                              cell_free_slots);
    } else {
        cells_.reserve(restored_cells.size());
    }
    for (std::size_t index = 0; index < restored_cells.size(); ++index) {
        const CellInit& cell = restored_cells[index];
        validate_restored_cell_schedule(cell, clock, config_);
        const Slot slot = has_explicit_layout
            ? cell_slots[index] : cells_.create(cell);
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
        density_.add(cell.anchor, cell.type, slot);
    }

    vessel_nodes_.reserve(vasculature.nodes.size());
    for (const VesselNodeInit3D& node : vasculature.nodes) {
        const VesselNodeSlot slot = vessel_nodes_.create(node);
        if (node.parent_node_slot != kEmptyVesselNodeSlot &&
            (!vessel_nodes_.valid(node.parent_node_slot) || node.parent_node_slot >= slot)) {
            throw std::runtime_error("checkpoint vessel parent slot is invalid");
        }
        centerline_nodes_.try_emplace(node.position, slot);
    }
    for (const VesselId vessel_id : vasculature.perfused_vessels) {
        perfused_vessels_.insert(vessel_id);
    }
    for (const VesselNodeSlot slot : vessel_nodes_.alive_slots()) {
        const VesselNodeInit3D node = vessel_nodes_.snapshot(slot);
        std::vector<Vec3i> capsule;
        if (node.parent_node_slot == kEmptyVesselNodeSlot) {
            capsule = rasterize_capsule(node.position, node.position, node.diameter_voxels);
        } else {
            capsule = rasterize_capsule(vessel_nodes_.position(node.parent_node_slot),
                                        node.position, node.diameter_voxels);
        }
        if (!vessel_grid_.all_in_domain(capsule)) {
            throw std::runtime_error("checkpoint vessel lies outside the configured domain");
        }
        for (const Vec3i site : capsule) {
            if (grid_.owner(site) != kEmptySlot) {
                throw std::runtime_error("checkpoint vessel overlaps a biological cell");
            }
        }
        const bool perfused = node.perfused || vessel_is_perfused(node.vessel_id);
        if (!vessel_grid_.add_sites(
                capsule, node.role, perfused, node.vessel_id).placed) {
            throw std::runtime_error("checkpoint vessel occupancy cannot be rebuilt");
        }
        if (perfused) {
            perfused_vessels_.insert(node.vessel_id);
            vessel_nodes_.set_perfused(slot, true);
            vascular_influence_.add_sources(capsule);
        }
    }

    vessel_tips_.reserve(vasculature.tips.size());
    for (const VesselTipInit3D& tip : vasculature.tips) {
        if (!vessel_nodes_.valid(tip.current_node_slot) ||
            vessel_nodes_.uid(tip.current_node_slot) != tip.current_node_uid) {
            throw std::runtime_error("checkpoint vessel tip references an invalid node");
        }
        vessel_tips_.create(tip);
    }

    next_uid_ = next_uid;
    next_vessel_id_ = vasculature.next_vessel_id;
    next_vessel_node_uid_ = vasculature.next_node_uid;
    next_vessel_tip_uid_ = vasculature.next_tip_uid;
    clock_ = clock;
    stats_ = stats;
    lineage_ = std::move(lineage);
    rebuild_migration_activation_class_cache();
    if (!config_.angiogenesis.enabled &&
        !vasculature.lesions.source_ownership.empty()) {
        throw std::runtime_error(
            "checkpoint has lesion source ownership while angiogenesis is disabled");
    }
    if (config_.angiogenesis.enabled) {
        rebuild_tumor_surface();
        rebuild_lesion_index();
        if (!vasculature.lesions.core_identity.empty() ||
            !vasculature.lesions.dirty_blocks.empty()) {
            lesion_index_.restore_checkpoint_state(
                vasculature.lesions.core_identity,
                vasculature.lesions.next_lesion_id,
                vasculature.lesions.dirty_blocks);
        } else if (vasculature.lesions.next_lesion_id >
                   lesion_index_.next_lesion_id()) {
            lesion_index_.set_next_lesion_id(
                vasculature.lesions.next_lesion_id);
        }
        if (!std::isfinite(vasculature.lesions.last_refresh_time_hours) ||
            vasculature.lesions.last_refresh_time_hours < 0.0 ||
            vasculature.lesions.last_refresh_time_hours > clock_.time_hours) {
            throw std::runtime_error(
                "checkpoint lesion refresh time is invalid");
        }
        last_lesion_refresh_time_hours_ =
            vasculature.lesions.last_refresh_time_hours;
        if (!std::isfinite(vasculature.lesions.next_refresh_time_hours) ||
            vasculature.lesions.next_refresh_time_hours < 0.0 ||
            (vasculature.lesions.next_refresh_time_hours > 0.0 &&
             vasculature.lesions.next_refresh_time_hours < clock_.time_hours &&
             !same_time(vasculature.lesions.next_refresh_time_hours,
                        clock_.time_hours))) {
            throw std::runtime_error(
                "checkpoint next lesion refresh time is invalid");
        }
        next_lesion_refresh_time_hours_ =
            vasculature.lesions.next_refresh_time_hours;
        lesion_refresh_schedule_generation_ =
            vasculature.lesions.refresh_schedule_generation;
        if (vasculature.lesions.dirty_blocks.empty() !=
            (next_lesion_refresh_time_hours_ == 0.0)) {
            throw std::runtime_error(
                "checkpoint lesion dirty state and refresh event disagree");
        }
        if (next_lesion_refresh_time_hours_ > 0.0) {
            const double expected = last_lesion_refresh_time_hours_ +
                config_.angiogenesis.lesion_refresh_interval_hours;
            if (!same_time(next_lesion_refresh_time_hours_, expected)) {
                throw std::runtime_error(
                    "checkpoint lesion refresh time does not match config interval");
            }
        }

        std::unordered_set<LesionId> current_lesion_ids;
        current_lesion_ids.reserve(lesion_index_.lesions().size());
        for (const LesionSummary3D& lesion : lesion_index_.lesions()) {
            current_lesion_ids.insert(lesion.id);
        }
        LesionId previous_source = kNoLesionId;
        for (const LesionSourceOwnership3D& ownership :
             vasculature.lesions.source_ownership) {
            if (ownership.source_lesion_id == kNoLesionId ||
                ownership.source_lesion_id >= lesion_index_.next_lesion_id() ||
                ownership.source_lesion_id <= previous_source ||
                ownership.current_lesion_id == ownership.source_lesion_id ||
                ownership.current_lesion_id >= lesion_index_.next_lesion_id() ||
                current_lesion_ids.contains(ownership.source_lesion_id) ||
                (ownership.current_lesion_id != kNoLesionId &&
                 !current_lesion_ids.contains(ownership.current_lesion_id))) {
                throw std::runtime_error(
                    "checkpoint lesion source ownership is invalid");
            }
            lesion_source_ownership_.emplace(
                ownership.source_lesion_id,
                ownership.current_lesion_id);
            previous_source = ownership.source_lesion_id;
        }
        const auto known_source = [this, &current_lesion_ids](LesionId source) {
            // source==0 remains available only to old in-memory geometry
            // fixtures that do not claim checkpoint provenance.
            return source == kNoLesionId ||
                   current_lesion_ids.contains(source) ||
                   lesion_source_ownership_.contains(source);
        };
        for (const VesselNodeSlot slot : vessel_nodes_.alive_slots()) {
            if (!known_source(vessel_nodes_.source_lesion_id(slot))) {
                throw std::runtime_error(
                    "vessel node historical source has no ownership record");
            }
        }
        for (const VesselTipSlot slot : vessel_tips_.alive_slots()) {
            if (!known_source(vessel_tips_.source_lesion_id(slot))) {
                throw std::runtime_error(
                    "vessel tip historical source has no ownership record");
            }
        }

        std::unordered_set<LesionId> restored_process_ids;
        for (const LesionAngiogenesisState3D& entry :
             vasculature.lesions.processes) {
            if (entry.lesion_id == kNoLesionId ||
                !restored_process_ids.insert(entry.lesion_id).second) {
                throw std::runtime_error(
                    "checkpoint lesion process ID is invalid or duplicated");
            }
            auto [found, inserted] = lesion_angiogenesis_processes_.try_emplace(
                entry.lesion_id, config_.seed, entry.lesion_id);
            if (!inserted) {
                throw std::runtime_error("duplicate checkpoint lesion process");
            }
            found->second.restore(entry.process);
            if (entry.process.eligible &&
                lesion_index_.find_lesion(entry.lesion_id) == nullptr) {
                throw std::runtime_error(
                    "eligible checkpoint process references a retired lesion");
            }
        }
        // In-memory/manual fixtures from the pre-lesion API may still supply
        // one aggregate process. Map it only when exactly one lesion exists;
        // versioned checkpoint v3 always writes explicit lesion processes.
        if (vasculature.lesions.processes.empty() &&
            (vasculature.process.eligible ||
             vasculature.process.event_sequence != 0 ||
             vasculature.process.attempted_events != 0) &&
            lesion_index_.lesions().size() == 1U) {
            const LesionId id = lesion_index_.lesions().front().id;
            auto [found, inserted] = lesion_angiogenesis_processes_.try_emplace(
                id, config_.seed, id);
            (void)inserted;
            found->second.restore(vasculature.process);
        }
        recompute_aggregate_angiogenesis_state();
    }
    initialized_ = true;

    for (const Slot slot : cells_.alive_slots()) {
        if (has_specific_generations(cells_.snapshot(slot))) restore_cell_events(slot);
        else schedule_cell(slot);
    }
    for (const VesselTipSlot slot : vessel_tips_.alive_slots()) restore_vessel_tip_event(slot);
    const bool initialize_missing_lesion_processes =
        config_.angiogenesis.enabled &&
        lesion_angiogenesis_processes_.empty() &&
        !lesion_index_.lesions().empty();
    if (initialize_missing_lesion_processes) {
        // Direct in-memory restores used by tests and embedding applications
        // may provide cells without serialized lesion clocks.  Start each
        // independently eligible lesion at the restored clock.  Versioned
        // checkpoint files always provide explicit process records and take
        // the branch below instead.
        sync_angiogenesis_eligibility(true);
    } else {
        for (const auto& [lesion_id, process] : lesion_angiogenesis_processes_) {
            if (process.state().eligible) schedule_seed_event(lesion_id);
        }
    }
    if (next_lesion_refresh_time_hours_ > 0.0) {
        schedule(EventKind::lesion_refresh, kEmptySlot, kNoLesionId,
                 next_lesion_refresh_time_hours_,
                 lesion_refresh_schedule_generation_);
    }
    reset_event_queue_rebuild_threshold();
}

void Simulation3D::run(const std::function<void(const Simulation3D&)>& observer) {
    if (!initialized_) initialize();
    if (observer) observer(*this);
    double last_observer_time = clock_.time_hours;
    while (clock_.completed_events < config_.max_events) {
        // Output is part of the event-driven clock, but never of the
        // biological event queue. If no biological event occurs for several
        // hours, advance through each requested output boundary and expose the
        // unchanged lazy state there. This produces exact timeline/checkpoint
        // times without drawing RNG or changing event ordering.
        while (!events_.empty() && !current(events_.top())) events_.pop();
        const double next_event_time = events_.empty()
            ? std::numeric_limits<double>::infinity()
            : events_.top().time;
        const double next_output_time = observer
            ? next_periodic_output_boundary(config_, clock_.time_hours)
            : std::numeric_limits<double>::infinity();
        if (next_output_time <= config_.end_time_hours &&
            next_event_time > next_output_time &&
            !same_time(next_event_time, next_output_time)) {
            clock_.time_hours = next_output_time;
            observer(*this);
            last_observer_time = clock_.time_hours;
            continue;
        }

        if (!step()) break;
        if (observer && (!same_time(clock_.time_hours, last_observer_time) || events_.empty())) {
            observer(*this);
            last_observer_time = clock_.time_hours;
        }
    }
}

bool Simulation3D::step() {
    if (!initialized_) initialize();
    if (clock_.completed_events >= config_.max_events) return false;
    while (!events_.empty() && !current(events_.top())) events_.pop();
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
        if (current(event)) batch.push_back(event);
    }
    if (batch.empty()) return !events_.empty();
    const std::uint64_t remaining_budget =
        config_.max_events - clock_.completed_events;
    if (batch.size() > remaining_budget) {
        // Proposal/conflict resolution is defined over the complete same-time
        // batch. Stop before it instead of executing a partial batch with
        // different conflict semantics or exceeding the hard event cap.
        for (const Event& event : batch) events_.push(event);
        return false;
    }
    clock_.time_hours = batch_time;

    std::vector<Event> death_events;
    std::vector<Event> vessel_events;
    std::vector<Event> division_events;
    std::vector<Event> migration_events;
    for (const EventKind kind : {EventKind::death, EventKind::lesion_refresh,
                                 EventKind::angiogenesis_seed,
                                 EventKind::vessel_growth,
                                 EventKind::migration_activation_end,
                                 EventKind::division, EventKind::migration}) {
        for (const Event& event : batch) {
            if (event.kind != kind) continue;
            if (kind == EventKind::migration_activation_end) {
                process_non_migration(event);
            } else if (kind == EventKind::lesion_refresh) {
                sync_angiogenesis_eligibility(true);
            } else if (kind == EventKind::death) {
                death_events.push_back(event);
            } else if (kind == EventKind::angiogenesis_seed) {
                (void)process_seed_event(event);
            } else if (kind == EventKind::vessel_growth) {
                vessel_events.push_back(event);
            } else if (kind == EventKind::division) {
                division_events.push_back(event);
            } else {
                migration_events.push_back(event);
            }
            if (kind != EventKind::death &&
                kind != EventKind::vessel_growth &&
                kind != EventKind::division &&
                kind != EventKind::migration) {
                ++clock_.completed_events;
            }
        }
        if (kind == EventKind::death && !death_events.empty()) {
            process_deaths(death_events);
            clock_.completed_events += death_events.size();
        }
        if (kind == EventKind::vessel_growth && !vessel_events.empty()) {
            process_vessel_growth(vessel_events);
            clock_.completed_events += vessel_events.size();
        }
        if (kind == EventKind::division && !division_events.empty()) {
            process_divisions(division_events);
            clock_.completed_events += division_events.size();
        }
        if (kind == EventKind::migration && !migration_events.empty()) {
            process_migrations(migration_events);
            clock_.completed_events += migration_events.size();
        }
    }
    sync_angiogenesis_eligibility();
    maybe_compact_event_queue();
    return true;
}

bool Simulation3D::current(const Event& event) const {
    if (event.kind == EventKind::lesion_refresh) {
        return config_.angiogenesis.enabled &&
               next_lesion_refresh_time_hours_ > 0.0 &&
               event.generation == lesion_refresh_schedule_generation_ &&
               same_time(event.time, next_lesion_refresh_time_hours_);
    }
    if (event.kind == EventKind::angiogenesis_seed) {
        if (!config_.angiogenesis.enabled) return false;
        const auto found = lesion_angiogenesis_processes_.find(event.uid);
        return found != lesion_angiogenesis_processes_.end() &&
               found->second.event_current(event.time, event.generation);
    }
    if (event.kind == EventKind::vessel_growth) {
        const VesselTipSlot slot = event.slot;
        return vessel_tips_.valid(slot) && vessel_tips_.uid(slot) == event.uid &&
               vessel_tips_.status(slot) == VesselTipStatus::active &&
               vessel_tips_.schedule_generation(slot) == event.generation &&
               vessel_tips_.pending_direction(slot) != kStayDirection &&
               same_time(vessel_tips_.next_growth_time(slot), event.time);
    }
    const Slot slot = event.slot;
    if (!cells_.valid(slot) || cells_.uid(slot) != event.uid) return false;
    if (event.kind == EventKind::migration_activation_end) {
        return (cells_.flags(slot) & static_cast<std::uint8_t>(kMigrationActive)) != 0 &&
               event.generation == cells_.migration_schedule_generation(slot) &&
               cells_.migration_activation_end_time(slot) > 0.0 &&
               same_time(cells_.migration_activation_end_time(slot), event.time);
    }
    const double expected = event_time(event.kind, slot);
    return event.generation == event_generation(event.kind, slot) &&
           expected > 0.0 && same_time(expected, event.time);
}

double Simulation3D::event_time(EventKind kind, Slot slot) const {
    switch (kind) {
        case EventKind::migration: return cells_.next_migration_time(slot);
        case EventKind::migration_activation_end:
            return cells_.migration_activation_end_time(slot);
        case EventKind::division: return cells_.next_division_time(slot);
        case EventKind::death: return cells_.death_deadline(slot);
        default: break;
    }
    throw std::logic_error("event kind is not a cell event");
}

std::uint32_t Simulation3D::event_generation(EventKind kind, Slot slot) const {
    switch (kind) {
        case EventKind::migration: return cells_.migration_schedule_generation(slot);
        case EventKind::migration_activation_end:
            return cells_.migration_schedule_generation(slot);
        case EventKind::division: return cells_.division_schedule_generation(slot);
        case EventKind::death: return cells_.death_schedule_generation(slot);
        default: break;
    }
    throw std::logic_error("event kind is not a cell event");
}

std::uint32_t Simulation3D::bump_event_generation(EventKind kind, Slot slot) {
    switch (kind) {
        case EventKind::migration: return cells_.bump_migration_schedule_generation(slot);
        case EventKind::division: return cells_.bump_division_schedule_generation(slot);
        case EventKind::death: return cells_.bump_death_schedule_generation(slot);
        default: break;
    }
    throw std::logic_error("event kind is not a cell event");
}

void Simulation3D::schedule(EventKind kind, std::uint32_t slot, std::uint64_t uid,
                            double time, std::uint32_t generation) {
    if (time > 0.0 && time >= clock_.time_hours && std::isfinite(time)) {
        events_.push({time, slot, uid, kind, generation});
    }
}

void Simulation3D::reset_event_queue_rebuild_threshold() {
    constexpr std::size_t kMinimumSlack = 256;
    const std::size_t size = events_.size();
    const std::size_t slack = std::max(kMinimumSlack, size / 4U);
    event_queue_rebuild_threshold_ =
        size > std::numeric_limits<std::size_t>::max() - slack
            ? std::numeric_limits<std::size_t>::max()
            : size + slack;
}

void Simulation3D::maybe_compact_event_queue() {
    if (events_.size() <= event_queue_rebuild_threshold_) return;

    // Generation changes intentionally make old heap entries inert. Rebuild
    // only after the heap has grown materially (25%, with a small fixed slack)
    // so stale entries remain bounded without an O(N) pass on every event.
    // Canonical next times/generations are copied from stores without drawing
    // RNG or bumping a generation, so maintenance is thread-neutral and cannot
    // affect future biology.
    events_ = decltype(events_){};
    for (const Slot slot : cells_.alive_slots()) restore_cell_events(slot);
    for (const VesselTipSlot slot : vessel_tips_.alive_slots()) {
        restore_vessel_tip_event(slot);
    }
    for (const auto& [lesion_id, process] : lesion_angiogenesis_processes_) {
        if (process.state().eligible) schedule_seed_event(lesion_id);
    }
    if (next_lesion_refresh_time_hours_ > 0.0) {
        schedule(EventKind::lesion_refresh, kEmptySlot, kNoLesionId,
                 next_lesion_refresh_time_hours_,
                 lesion_refresh_schedule_generation_);
    }
    ++event_queue_rebuild_count_;
    reset_event_queue_rebuild_threshold();
}

void Simulation3D::schedule_cell(Slot slot) {
    if (!cells_.valid(slot)) return;
    reschedule_event(EventKind::migration, slot);
    reschedule_event(EventKind::division, slot);
    reschedule_event(EventKind::death, slot);
}

void Simulation3D::restore_cell_events(Slot slot) {
    if (!cells_.valid(slot)) return;
    const std::uint32_t migration_generation =
        cells_.migration_schedule_generation(slot);
    schedule(EventKind::migration, slot, cells_.uid(slot),
             cells_.next_migration_time(slot), migration_generation);
    schedule(EventKind::migration_activation_end, slot, cells_.uid(slot),
             cells_.migration_activation_end_time(slot), migration_generation);
    for (const EventKind kind : {EventKind::division, EventKind::death}) {
        schedule(kind, slot, cells_.uid(slot), event_time(kind, slot),
                 event_generation(kind, slot));
    }
}

void Simulation3D::reschedule_event(EventKind kind, Slot slot) {
    if (!cells_.valid(slot)) return;
    if (kind == EventKind::migration_activation_end) {
        throw std::logic_error(
            "migration activation end must be rescheduled with migration");
    }
    const std::uint32_t generation = bump_event_generation(kind, slot);
    schedule(kind, slot, cells_.uid(slot), event_time(kind, slot), generation);
    if (kind == EventKind::migration) {
        schedule(EventKind::migration_activation_end, slot, cells_.uid(slot),
                 cells_.migration_activation_end_time(slot), generation);
    }
}

void Simulation3D::synchronize_migration_schedule(Slot slot,
                                                  bool activation_changed) {
    if (!activation_changed || !cells_.valid(slot)) return;
    const double rate = effective_migration_rate(slot, cells_, config_);
    cells_.set_next_migration_time(
        slot, migration_allowed_for_cell(slot, cells_, config_) && rate > 0.0
                  ? clock_.time_hours + 1.0 / rate
                  : 0.0);
    reschedule_event(EventKind::migration, slot);
}

void Simulation3D::apply_growth_refresh(Slot slot,
                                        const GrowthRefreshResult& refresh) {
    if (!cells_.valid(slot)) return;
    if (refresh.division_time_changed) {
        reschedule_event(EventKind::division, slot);
    }
    if (refresh.death_time_changed) {
        reschedule_event(EventKind::death, slot);
    }
    synchronize_migration_schedule(slot, refresh.migration_activation_changed);
}

void Simulation3D::process_non_migration(const Event& event) {
    if (!current(event)) return;
    if (event.kind != EventKind::migration_activation_end) {
        throw std::logic_error("event kind requires a same-time batch processor");
    }
    if (expire_migration_activation_state(
            event.slot, clock_.time_hours, cells_, config_)) {
        const double rate = effective_migration_rate(
            event.slot, cells_, config_);
        cells_.set_next_migration_time(
            event.slot, rate > 0.0 ? clock_.time_hours + 1.0 / rate : 0.0);
        reschedule_event(EventKind::migration, event.slot);
    }
}

void Simulation3D::process_deaths(const std::vector<Event>& events) {
    struct DeathDecision {
        Event event;
        bool remove{};
    };
    std::vector<Event> pending;
    pending.reserve(events.size());
    for (const Event& event : events) {
        if (current(event)) pending.push_back(event);
    }
    std::vector<DeathDecision> decisions(pending.size());
    std::vector<Event> removals;
    removals.reserve(events.size());
    const VascularInfluenceField3D* influence = config_.angiogenesis.enabled
        ? &vascular_influence_ : nullptr;

    // Decide every same-time death against the same occupancy/density
    // snapshot. No removal is visible while another death is being judged.
    // This read-only phase is safe to parallelize; refresh and commit remain
    // ordered below so thread scheduling cannot affect model evolution.
#ifdef _OPENMP
    const int worker_count = select_worker_count(
        cells_.alive_count(), pending.size(), config_, available_worker_threads());
#pragma omp parallel for schedule(static) num_threads(worker_count)
#endif
    for (std::int64_t index = 0;
         index < static_cast<std::int64_t>(pending.size()); ++index) {
        const Event& event = pending[static_cast<std::size_t>(index)];
        const double rate = density_growth_rate_for_cell(
            cells_, event.slot, density_, config_, influence);
        decisions[static_cast<std::size_t>(index)] = {
            event, rate <= config_.death_growth_rate_threshold};
    }
    for (const DeathDecision& decision : decisions) {
        const Event& event = decision.event;
        if (!decision.remove) {
            const GrowthRefreshResult refresh = refresh_growth_state(
                event.slot, clock_.time_hours, cells_, density_, config_, influence);
            apply_growth_refresh(event.slot, refresh);
        } else {
            removals.push_back(event);
        }
    }

    std::vector<Vec3i> changed_sites;
    for (const Event& event : removals) {
        if (!cells_.valid(event.slot) || cells_.uid(event.slot) != event.uid) continue;
        const std::vector<Vec3i> occupied = occupied_sites_for_cell(cells_, event.slot);
        if (remove_cell(event.slot, cells_, grid_, density_)) {
            ++stats_.deaths;
            changed_sites.insert(changed_sites.end(), occupied.begin(), occupied.end());
        }
    }
    if (!changed_sites.empty()) {
        recover_neighborhood(changed_sites);
        refresh_tumor_surface(changed_sites);
        refresh_neighborhood(changed_sites);
    }
}

void Simulation3D::process_divisions(const std::vector<Event>& events) {
    struct EligibleDivision {
        Event event;
        double mother_death_before{};
        std::uint64_t priority{};
    };
    struct OrderedDivision {
        Event event;
        DivisionProposal proposal;
        double mother_death_before{};
        std::uint64_t priority{};
    };

    std::vector<EligibleDivision> eligible;
    eligible.reserve(events.size());
    const VascularInfluenceField3D* influence = config_.angiogenesis.enabled
        ? &vascular_influence_ : nullptr;
    for (const Event& event : events) {
        if (!current(event)) continue;
        const GrowthRefreshResult event_refresh = refresh_growth_state(
            event.slot, clock_.time_hours, cells_, density_, config_, influence);
        if (event_refresh.death_time_changed) {
            reschedule_event(EventKind::death, event.slot);
        }
        synchronize_migration_schedule(
            event.slot, event_refresh.migration_activation_changed);
        if (cells_.density_growth_rate(event.slot) <=
                config_.death_growth_rate_threshold ||
            cells_.division_work_remaining(event.slot) > 1.0e-5F) {
            reschedule_event(EventKind::division, event.slot);
            continue;
        }
        eligible.push_back({
            event,
            cells_.death_deadline(event.slot),
            division_conflict_priority(config_, event.time, event.uid),
        });
    }
    std::vector<OrderedDivision> ordered(eligible.size());
#ifdef _OPENMP
    const int worker_count = select_worker_count(
        cells_.alive_count(), eligible.size(), config_, available_worker_threads());
#pragma omp parallel for schedule(static) num_threads(worker_count)
#endif
    for (std::int64_t index = 0;
         index < static_cast<std::int64_t>(eligible.size()); ++index) {
        const EligibleDivision& candidate = eligible[static_cast<std::size_t>(index)];
        ordered[static_cast<std::size_t>(index)] = OrderedDivision{
            candidate.event,
            make_division_proposal(
                candidate.event.slot, cells_, grid_, config_),
            candidate.mother_death_before,
            candidate.priority,
        };
    }
    std::sort(ordered.begin(), ordered.end(), [](const OrderedDivision& lhs,
                                                  const OrderedDivision& rhs) {
        if (lhs.priority != rhs.priority) return lhs.priority > rhs.priority;
        return lhs.event.uid < rhs.event.uid;
    });

    std::unordered_set<Vec3i, Vec3iHash> reserved;
    std::unordered_set<Vec3i, Vec3iHash> locked_colocation_groups;
    std::vector<Vec3i> changed_sites;
    bool needs_opportunistic_recovery = false;

    // Every proposal and priority is frozen against one occupancy snapshot.
    // Resolve complete voxel footprints and co-location groups before any
    // biological fallback can run. A loser only retries; it cannot die as an
    // r-cell or create a K-cell co-location merely because an earlier winner
    // mutated the grid.
    for (const OrderedDivision& contender : ordered) {
        const DivisionProposal& proposal = contender.proposal;
        bool conflict = std::any_of(
            proposal.reserved_sites.begin(), proposal.reserved_sites.end(),
            [&reserved](Vec3i site) { return reserved.contains(site); });
        if (!conflict && proposal.locks_colocation_group &&
            locked_colocation_groups.contains(proposal.mother_anchor)) {
            conflict = true;
        }
        if (conflict) {
            ++stats_.conflict_rejections;
            if (cells_.valid(contender.event.slot)) {
                cells_.set_next_division_time(
                    contender.event.slot,
                    clock_.time_hours + config_.division_timing.retry_delay_hours);
                reschedule_event(EventKind::division, contender.event.slot);
            }
            continue;
        }
        reserved.insert(proposal.reserved_sites.begin(),
                        proposal.reserved_sites.end());
        if (proposal.locks_colocation_group) {
            locked_colocation_groups.insert(proposal.mother_anchor);
        }

        DivisionResult result = commit_division_proposal(
            proposal, clock_.time_hours, next_uid_, cells_, grid_, density_,
            config_, lineage_, influence);
        if (result.changed) {
            needs_opportunistic_recovery =
                needs_opportunistic_recovery || !result.stage_recovery;
            changed_sites.insert(changed_sites.end(), result.changed_sites.begin(),
                                 result.changed_sites.end());
            if (result.daughter != kEmptySlot) {
                ++stats_.divisions;
                const double rate = effective_migration_rate(
                    result.daughter, cells_, config_);
                cells_.set_next_migration_time(
                    result.daughter,
                    migration_allowed_for_cell(result.daughter, cells_, config_) &&
                            rate > 0.0
                        ? clock_.time_hours + 1.0 / rate
                        : 0.0);
                schedule_cell(result.daughter);
                if (cells_.valid(contender.event.slot)) {
                    reschedule_event(EventKind::division, contender.event.slot);
                    if (!same_time(
                            contender.mother_death_before,
                            cells_.death_deadline(contender.event.slot))) {
                        reschedule_event(EventKind::death, contender.event.slot);
                    }
                    const double mother_rate = effective_migration_rate(
                        contender.event.slot, cells_, config_);
                    cells_.set_next_migration_time(
                        contender.event.slot,
                        mother_rate > 0.0
                            ? clock_.time_hours + 1.0 / mother_rate
                            : 0.0);
                    reschedule_event(EventKind::migration, contender.event.slot);
                }
            } else if (result.mother_removed) {
                ++stats_.deaths;
            } else if (cells_.valid(contender.event.slot)) {
                cells_.set_next_division_time(
                    contender.event.slot,
                    clock_.time_hours + config_.division_timing.retry_delay_hours);
                reschedule_event(EventKind::division, contender.event.slot);
            }
        } else if (cells_.valid(contender.event.slot)) {
            if (proposal.action != DivisionAction::none) {
                ++stats_.conflict_rejections;
            }
            cells_.set_next_division_time(
                contender.event.slot,
                clock_.time_hours + config_.division_timing.retry_delay_hours);
            reschedule_event(EventKind::division, contender.event.slot);
        }
    }
    if (!changed_sites.empty()) {
        if (needs_opportunistic_recovery) recover_neighborhood(changed_sites);
        refresh_tumor_surface(changed_sites);
        refresh_neighborhood(changed_sites);
    }
}

void Simulation3D::process_migrations(const std::vector<Event>& events) {
    struct Pending { Event event; std::uint64_t sequence{}; };
    std::vector<Pending> pending;
    pending.reserve(events.size());
    for (const Event& event : events) {
        if (!current(event)) continue;
        if (refresh_migration_activation_state(
                event.slot, clock_.time_hours, cells_, density_, config_)) {
            const double rate = effective_migration_rate(
                event.slot, cells_, config_);
            cells_.set_next_migration_time(
                event.slot, rate > 0.0 ? clock_.time_hours + 1.0 / rate : 0.0);
            reschedule_event(EventKind::migration, event.slot);
            continue;
        }
        pending.push_back({event, cells_.consume_event_sequence(event.slot)});
    }
    std::vector<MoveProposal> proposals(pending.size());
    const std::uint64_t time_bucket = config_.conflict_bucket_hours > 0.0
        ? static_cast<std::uint64_t>(std::floor(clock_.time_hours / config_.conflict_bucket_hours))
        : std::bit_cast<std::uint64_t>(clock_.time_hours);
#ifdef _OPENMP
    const int worker_count = select_worker_count(
        cells_.alive_count(), pending.size(), config_, available_worker_threads());
#pragma omp parallel for schedule(static) num_threads(worker_count)
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
            if (cells_.valid(proposal.slot)) cells_.set_last_direction(proposal.slot, kStayDirection);
        } else {
            const bool conflict = std::any_of(
                proposal.reserved_sites.begin(), proposal.reserved_sites.end(),
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
            const double rate = effective_migration_rate(
                proposal.slot, cells_, config_);
            cells_.set_next_migration_time(proposal.slot,
                migration_allowed_for_cell(proposal.slot, cells_, config_) && rate > 0.0
                    ? clock_.time_hours + 1.0 / rate : 0.0);
            reschedule_event(EventKind::migration, proposal.slot);
        }
    }
    if (!changed_sites.empty()) {
        recover_neighborhood(changed_sites);
        refresh_tumor_surface(changed_sites);
        refresh_neighborhood(changed_sites);
    }
}

void Simulation3D::rebuild_tumor_surface() {
    if (!config_.angiogenesis.enabled) return;
    tumor_surface_.rebuild_from_visitor(
        [this](auto&& visitor) {
            grid_.for_each_occupied_site(
                [&](Vec3i site, Slot) { visitor(site); });
        },
        [this](Vec3i site) { return grid_.owner(site) != kEmptySlot; });
}

void Simulation3D::rebuild_lesion_index() {
    if (!config_.angiogenesis.enabled) return;
    lesion_index_.rebuild_from(
        cells_, grid_, lesion_biological_volumes(config_));
    last_lesion_refresh_time_hours_ = clock_.time_hours;
}

void Simulation3D::mark_lesion_dirty(
    std::span<const Vec3i> changed_sites) {
    if (!config_.angiogenesis.enabled || changed_sites.empty()) return;
    std::vector<Vec3i> expanded;
    expanded.reserve(changed_sites.size() * 8U);
    for (const Vec3i site : changed_sites) {
        for (int dx = 0; dx <= 1; ++dx) {
            for (int dy = 0; dy <= 1; ++dy) {
                for (int dz = 0; dz <= (config_.thin_layer ? 0 : 1); ++dz) {
                    expanded.push_back(site + Vec3i{dx, dy, dz});
                }
            }
        }
    }
    lesion_index_.mark_dirty_sites(expanded);
}

bool Simulation3D::refresh_lesion_index(bool force) {
    if (!config_.angiogenesis.enabled) return false;
    const double interval = config_.angiogenesis.lesion_refresh_interval_hours;
    const bool due = force ||
        clock_.time_hours + 1e-12 >= last_lesion_refresh_time_hours_ + interval;
    if (!due) return false;

    next_lesion_refresh_time_hours_ = 0.0;
    ++lesion_refresh_schedule_generation_;

    bool changed = false;
    if (lesion_index_.observations_dirty()) {
        changed = lesion_index_.refresh_dirty_blocks_from(
            cells_, grid_, lesion_biological_volumes(config_)) != 0U;
    }
    if (lesion_index_.topology_dirty() || lesion_index_.statistics_dirty()) {
        std::vector<LesionId> old_lesion_ids;
        old_lesion_ids.reserve(lesion_index_.lesions().size());
        for (const LesionSummary3D& lesion : lesion_index_.lesions()) {
            old_lesion_ids.push_back(lesion.id);
        }
        std::sort(old_lesion_ids.begin(), old_lesion_ids.end());

        const LesionTopologyDelta3D delta = lesion_index_.refresh_topology();

        std::unordered_set<LesionId> current_lesions;
        current_lesions.reserve(lesion_index_.lesions().size());
        for (const LesionSummary3D& lesion : lesion_index_.lesions()) {
            current_lesions.insert(lesion.id);
        }

        // Assign every former lesion to at most one successor.  A retained
        // child owns its process even when the same predecessor also appears
        // in a simultaneous merge.  Removed split parents use the smallest
        // resulting child, and remaining merge-only parents use the smallest
        // merge result.  This deterministic one-to-one transfer prevents both
        // lost state and duplicated counters.
        std::unordered_map<LesionId, LesionId> successors;
        successors.reserve(old_lesion_ids.size());
        for (const LesionId old_id : old_lesion_ids) {
            successors.emplace(
                old_id, current_lesions.contains(old_id) ? old_id : kNoLesionId);
        }
        for (const LesionSplit3D& split : delta.splits) {
            auto found = successors.find(split.predecessor);
            if (found == successors.end() || found->second != kNoLesionId) continue;
            for (const LesionId child : split.children) {
                if (current_lesions.contains(child) &&
                    (found->second == kNoLesionId || child < found->second)) {
                    found->second = child;
                }
            }
        }
        for (const LesionMerge3D& merge : delta.merges) {
            if (!current_lesions.contains(merge.result)) continue;
            for (const LesionId predecessor : merge.predecessors) {
                auto found = successors.find(predecessor);
                if (found != successors.end() &&
                    found->second == kNoLesionId) {
                    found->second = merge.result;
                }
            }
        }

        std::unordered_map<LesionId,
                           std::vector<LesionAngiogenesisState3D>> incoming;
        for (const LesionId old_id : old_lesion_ids) {
            const auto process = lesion_angiogenesis_processes_.find(old_id);
            if (process == lesion_angiogenesis_processes_.end()) continue;
            const LesionId successor = successors.at(old_id);
            if (successor == old_id) continue;
            if (successor == kNoLesionId) {
                if (process->second.state().eligible) {
                    process->second.stop(clock_.time_hours);
                }
                continue;
            }
            incoming[successor].push_back({old_id, process->second.state()});
            lesion_angiogenesis_processes_.erase(process);
        }

        std::vector<LesionId> destinations;
        destinations.reserve(incoming.size());
        for (const auto& [destination, states] : incoming) {
            (void)states;
            destinations.push_back(destination);
        }
        std::sort(destinations.begin(), destinations.end());
        for (const LesionId destination : destinations) {
            std::vector<LesionAngiogenesisState3D>& states = incoming.at(destination);
            const auto existing = lesion_angiogenesis_processes_.find(destination);
            if (existing != lesion_angiogenesis_processes_.end()) {
                states.push_back({destination, existing->second.state()});
            }
            const AngiogenesisProcessState3D merged =
                merge_angiogenesis_process_states(states, clock_.time_hours);
            auto [found, inserted] =
                lesion_angiogenesis_processes_.try_emplace(
                    destination, config_.seed, destination);
            (void)inserted;
            found->second.restore(merged);
            if (merged.eligible && merged.next_seed_time_hours > 0.0) {
                schedule_seed_event(destination);
            }
        }

        update_lesion_source_ownership(successors, current_lesions);
        changed = true;
    }
    last_lesion_refresh_time_hours_ = clock_.time_hours;
    return changed || force;
}

void Simulation3D::schedule_lesion_refresh_event() {
    if (!config_.angiogenesis.enabled ||
        !lesion_index_.refresh_needed() ||
        next_lesion_refresh_time_hours_ > 0.0) {
        return;
    }
    const double due = last_lesion_refresh_time_hours_ +
                       config_.angiogenesis.lesion_refresh_interval_hours;
    next_lesion_refresh_time_hours_ = std::max(due, clock_.time_hours);
    ++lesion_refresh_schedule_generation_;
    schedule(EventKind::lesion_refresh, kEmptySlot, kNoLesionId,
             next_lesion_refresh_time_hours_,
             lesion_refresh_schedule_generation_);
}

void Simulation3D::refresh_tumor_surface(std::span<const Vec3i> changed_sites) {
    if (!config_.angiogenesis.enabled || changed_sites.empty()) return;
    std::vector<Vec3i> expanded;
    expanded.reserve(changed_sites.size() * 8U);
    for (const Vec3i anchor : changed_sites) {
        for (int dx = 0; dx <= 1; ++dx) {
            for (int dy = 0; dy <= 1; ++dy) {
                for (int dz = 0; dz <= (config_.thin_layer ? 0 : 1); ++dz) {
                    expanded.push_back(anchor + Vec3i{dx, dy, dz});
                }
            }
        }
    }
    std::sort(expanded.begin(), expanded.end());
    expanded.erase(std::unique(expanded.begin(), expanded.end()), expanded.end());
    tumor_surface_.refresh(expanded,
        [this](Vec3i site) { return grid_.owner(site) != kEmptySlot; });
    lesion_index_.mark_dirty_sites(expanded);
}

double Simulation3D::biological_tumor_volume() const noexcept {
    return static_cast<double>(cells_.stage_count(CellStage::large)) *
               config_.angiogenesis.stage0_biological_volume_voxels3 +
           static_cast<double>(cells_.stage_count(CellStage::small)) *
               config_.angiogenesis.stage1_biological_volume_voxels3 +
           static_cast<double>(cells_.stage_count(CellStage::ultrasmall)) *
               config_.angiogenesis.stage2_biological_volume_voxels3;
}

void Simulation3D::recompute_aggregate_angiogenesis_state() {
    std::vector<LesionAngiogenesisState3D> processes;
    processes.reserve(lesion_angiogenesis_processes_.size());
    for (const auto& [lesion_id, process] : lesion_angiogenesis_processes_) {
        processes.push_back({lesion_id, process.state()});
    }
    aggregate_angiogenesis_state_ = aggregate_angiogenesis_process_states(
        processes, clock_.time_hours);
}

LesionId Simulation3D::current_lesion_for_source(
    LesionId source_lesion_id) const noexcept {
    const auto alias = lesion_source_ownership_.find(source_lesion_id);
    if (alias != lesion_source_ownership_.end()) return alias->second;
    return lesion_index_.find_lesion(source_lesion_id) != nullptr
        ? source_lesion_id : kNoLesionId;
}

void Simulation3D::update_lesion_source_ownership(
    const std::unordered_map<LesionId, LesionId>& successors,
    const std::unordered_set<LesionId>& current_lesions) {
    for (auto& [source, owner] : lesion_source_ownership_) {
        (void)source;
        if (owner == kNoLesionId) continue;
        const auto successor = successors.find(owner);
        if (successor != successors.end()) {
            owner = successor->second;
        } else if (!current_lesions.contains(owner)) {
            owner = kNoLesionId;
        }
    }
    for (const auto& [source, owner] : successors) {
        if (owner == source && current_lesions.contains(source)) {
            lesion_source_ownership_.erase(source);
        } else {
            lesion_source_ownership_[source] = owner;
        }
    }
    // Current IDs always own themselves.  All other targets must be direct
    // current owners (or the explicit retired sentinel), never alias chains.
    for (const LesionId current : current_lesions) {
        lesion_source_ownership_.erase(current);
    }
    for (auto& [source, owner] : lesion_source_ownership_) {
        (void)source;
        if (owner != kNoLesionId && !current_lesions.contains(owner)) {
            owner = kNoLesionId;
        }
    }
}

void Simulation3D::sync_angiogenesis_eligibility(bool force_refresh) {
    if (!initialized_ || !config_.angiogenesis.enabled) return;
    const bool refreshed = refresh_lesion_index(force_refresh);
    if (!refreshed && !force_refresh) {
        schedule_lesion_refresh_event();
        return;
    }

    std::unordered_set<LesionId> current_lesions;
    current_lesions.reserve(lesion_index_.lesions().size());
    for (const LesionSummary3D& lesion : lesion_index_.lesions()) {
        current_lesions.insert(lesion.id);
    }
    for (auto& [lesion_id, process] : lesion_angiogenesis_processes_) {
        if (!current_lesions.contains(lesion_id) && process.state().eligible) {
            process.stop(clock_.time_hours);
        }
    }

    recompute_aggregate_angiogenesis_state();
    const bool global_limit_reached =
        aggregate_angiogenesis_state_.committed_roots >=
        config_.angiogenesis.max_total_roots;
    for (const LesionSummary3D& lesion : lesion_index_.lesions()) {
        auto [found, inserted] = lesion_angiogenesis_processes_.try_emplace(
            lesion.id, config_.seed, lesion.id);
        (void)inserted;
        AngiogenesisProcess3D& process = found->second;
        const bool geometry_eligible =
            lesion.core_blocks.size() >=
            config_.angiogenesis.trigger_minimum_core_blocks;
        if (!geometry_eligible) {
            if (process.state().eligible) process.stop(clock_.time_hours);
            continue;
        }
        const bool changed = process.update_volume(
            clock_.time_hours, lesion.biological_volume,
            config_.angiogenesis.trigger_activation_volume_voxels3,
            config_.angiogenesis.trigger_deactivation_volume_voxels3,
            config_.angiogenesis.trigger_delay_hours,
            config_.angiogenesis.seed_rate_sites_per_30_days);
        if ((global_limit_reached ||
             process.state().committed_roots >=
                 config_.angiogenesis.max_roots_per_lesion) &&
            process.state().eligible) {
            process.stop(clock_.time_hours);
        } else if (changed && process.state().eligible) {
            schedule_seed_event(lesion.id);
        }
    }
    if (global_limit_reached) {
        for (auto& [lesion_id, process] : lesion_angiogenesis_processes_) {
            (void)lesion_id;
            if (process.state().eligible) process.stop(clock_.time_hours);
        }
    }
    recompute_aggregate_angiogenesis_state();
}

void Simulation3D::schedule_seed_event(LesionId lesion_id) {
    const auto found = lesion_angiogenesis_processes_.find(lesion_id);
    if (found == lesion_angiogenesis_processes_.end()) return;
    const auto& state = found->second.state();
    schedule(EventKind::angiogenesis_seed, kEmptySlot, lesion_id,
             state.next_seed_time_hours, state.schedule_generation);
}

bool Simulation3D::process_seed_event(const Event& event) {
    if (!current(event)) return false;
    sync_angiogenesis_eligibility(true);
    if (!current(event)) return false;
    auto process_found = lesion_angiogenesis_processes_.find(event.uid);
    const LesionSummary3D* lesion = lesion_index_.find_lesion(event.uid);
    if (process_found == lesion_angiogenesis_processes_.end() || lesion == nullptr) {
        return false;
    }
    AngiogenesisProcess3D& process = process_found->second;
    ++stats_.angiogenesis_seed_attempts;
    bool committed = false;
    if (active_vessel_tip_count() + 2U <= config_.angiogenesis.max_active_tips &&
        active_vessel_tip_count(event.uid) + 2U <=
            config_.angiogenesis.max_active_tips_per_lesion &&
        aggregate_angiogenesis_state_.committed_roots <
            config_.angiogenesis.max_total_roots &&
        process.state().committed_roots <
            config_.angiogenesis.max_roots_per_lesion &&
        !tumor_surface_.empty()) {
        const auto candidates =
            tumor_surface_.sample_external_subset_without_replacement(
            std::min<std::size_t>(config_.angiogenesis.surface_max_sampling_attempts,
                                  tumor_surface_.size()),
            0.0, splitmix64(config_.seed ^ event.uid),
            process.state().attempted_events,
            [this, lesion_id = event.uid](const ExposedFace3D& face) {
                const auto owner = lesion_index_.lesion_for_face(
                    face, cells_, grid_);
                return owner.has_value() && *owner == lesion_id;
            });
        const std::int64_t minimum_squared =
            static_cast<std::int64_t>(config_.angiogenesis.surface_min_separation_voxels) *
            config_.angiogenesis.surface_min_separation_voxels;
        std::vector<Vec3i> existing_roots;
        for (const VesselNodeSlot slot : vessel_nodes_.alive_slots()) {
            if (vessel_nodes_.role(slot) == VesselBranchRole::root) {
                existing_roots.push_back(vessel_nodes_.position(slot));
            }
        }
        const Vec3i inward_target = rounded_lesion_centroid(*lesion);
        for (const ExposedFace3D& face : candidates) {
            const Vec3i root = face.inside;
            if (grid_.owner(face.inside) == kEmptySlot ||
                grid_.owner(face.outside()) != kEmptySlot) continue;
            const auto face_lesion = lesion_index_.lesion_for_face(
                face, cells_, grid_);
            if (!face_lesion || *face_lesion != event.uid) continue;
            if (!vessel_grid_.in_domain(root) || centerline_nodes_.contains(root)) continue;
            bool separated = true;
            for (const Vec3i existing : existing_roots) {
                if (squared_distance(root, existing) < minimum_squared) {
                    separated = false;
                    break;
                }
            }
            if (separated && create_vessel_root(
                    event.uid, face, inward_target,
                    process.state().attempted_events)) {
                committed = true;
                break;
            }
        }
    }
    process.consume_event(
        clock_.time_hours, committed, config_.angiogenesis.seed_rate_sites_per_30_days);
    if (committed) {
        ++stats_.angiogenesis_roots;
    } else {
        ++stats_.angiogenesis_seed_rejections;
    }
    recompute_aggregate_angiogenesis_state();
    if (aggregate_angiogenesis_state_.committed_roots >=
            config_.angiogenesis.max_total_roots ||
        process.state().committed_roots >=
            config_.angiogenesis.max_roots_per_lesion) {
        process.stop(clock_.time_hours);
    } else {
        schedule_seed_event(event.uid);
    }
    recompute_aggregate_angiogenesis_state();
    return committed;
}

bool Simulation3D::root_has_local_support(
    LesionId source_lesion_id,
    std::span<const Vec3i> root_capsule) const {
    std::unordered_set<Slot> displaced;
    bool intersects_source = false;
    for (const Vec3i site : root_capsule) {
        for (const Slot slot : grid_.occupants(site)) {
            if (!cells_.valid(slot)) continue;
            displaced.insert(slot);
            const auto owner = lesion_index_.lesion_for_anchor(cells_.anchor(slot));
            intersects_source = intersects_source ||
                (owner.has_value() && *owner == source_lesion_id);
        }
    }
    if (!intersects_source) return false;

    std::unordered_set<Slot> support;
    for (const Vec3i site : root_capsule) {
        for (const Vec3i normal : kAxisFaceNormals3D) {
            for (const Slot slot : grid_.occupants(site + normal)) {
                if (!cells_.valid(slot) || displaced.contains(slot)) continue;
                const auto owner =
                    lesion_index_.lesion_for_anchor(cells_.anchor(slot));
                if (owner.has_value() && *owner == source_lesion_id) {
                    support.insert(slot);
                }
            }
        }
    }
    return support.size() >= config_.angiogenesis.surface_min_local_cells;
}

bool Simulation3D::create_vessel_root(LesionId source_lesion_id,
                                      const ExposedFace3D& face,
                                      Vec3i inward_target,
                                      std::uint64_t seed_event_sequence) {
    (void)seed_event_sequence;
    const auto lesion = lesion_index_.lesion_for_face(face, cells_, grid_);
    if (!lesion || *lesion != source_lesion_id ||
        grid_.owner(face.inside) == kEmptySlot ||
        grid_.owner(face.outside()) != kEmptySlot) return false;
    const Vec3i root_position = face.inside;
    const float diameter = static_cast<float>(config_.angiogenesis.diameter_voxels);
    const std::vector<Vec3i> root_capsule =
        rasterize_capsule(root_position, root_position, diameter);
    if (!vessel_grid_.all_in_domain(root_capsule)) return false;
    // A seed event is allowed to displace biological cells, but it must not
    // silently lay a second root over an existing vascular capsule.  Perform
    // this preflight before consuming an id or mutating either occupancy
    // layer so rejected Poisson events leave the model state unchanged.
    if (vessel_grid_.any_occupied(root_capsule)) return false;
    if (!root_has_local_support(source_lesion_id, root_capsule)) return false;

    const VesselId vessel_id = next_vessel_id_;
    const bool immediate = config_.angiogenesis.influence_activation == "immediate";
    std::vector<Vec3i> displaced = displace_cells_for_vessel(root_capsule);
    if (!vessel_grid_.add_sites(
            root_capsule, VesselBranchRole::root, immediate, vessel_id).placed) {
        throw std::logic_error("preflighted vessel root placement failed");
    }
    ++next_vessel_id_;

    VesselNodeInit3D root;
    root.position = root_position;
    root.uid = next_vessel_node_uid_++;
    root.vessel_id = vessel_id;
    root.source_lesion_id = source_lesion_id;
    root.role = VesselBranchRole::root;
    root.perfused = immediate;
    root.diameter_voxels = diameter;
    root.created_time_hours = clock_.time_hours;
    const VesselNodeSlot root_slot = vessel_nodes_.create(root);
    centerline_nodes_.emplace(root_position, root_slot);
    if (immediate) {
        perfused_vessels_.insert(vessel_id);
        vascular_influence_.add_sources(root_capsule);
    }

    Vec3i inward_axis = inward_target - root_position;
    if (squared_length(inward_axis) == 0) {
        inward_axis = {-face.outward_normal.x, -face.outward_normal.y,
                       -face.outward_normal.z};
    }
    VesselTipInit3D inward;
    inward.position = root_position;
    inward.bias_axis = inward_axis;
    inward.target = inward_target;
    inward.uid = next_vessel_tip_uid_++;
    inward.vessel_id = vessel_id;
    inward.source_lesion_id = source_lesion_id;
    inward.current_node_uid = root.uid;
    inward.current_node_slot = root_slot;
    inward.role = VesselBranchRole::inward;
    inward.status = VesselTipStatus::active;
    inward.perfused = immediate;
    inward.diameter_voxels = diameter;
    inward.speed_voxels_per_hour =
        static_cast<float>(config_.angiogenesis.inward_speed_voxels_per_hour);
    inward.max_length_voxels =
        static_cast<float>(config_.angiogenesis.inward_max_length_voxels);
    const VesselTipSlot inward_slot = vessel_tips_.create(inward);

    VesselTipInit3D outward = inward;
    outward.bias_axis = face.outward_normal;
    outward.target = root_position + Vec3i{
        face.outward_normal.x * config_.angiogenesis.outward_max_length_voxels,
        face.outward_normal.y * config_.angiogenesis.outward_max_length_voxels,
        face.outward_normal.z * config_.angiogenesis.outward_max_length_voxels};
    outward.uid = next_vessel_tip_uid_++;
    outward.role = VesselBranchRole::outward;
    outward.speed_voxels_per_hour =
        static_cast<float>(config_.angiogenesis.outward_speed_voxels_per_hour);
    outward.max_length_voxels =
        static_cast<float>(config_.angiogenesis.outward_max_length_voxels);
    const VesselTipSlot outward_slot = vessel_tips_.create(outward);

    schedule_vessel_tip(inward_slot);
    schedule_vessel_tip(outward_slot);
    if (!displaced.empty()) {
        recover_neighborhood(displaced);
        refresh_tumor_surface(displaced);
        refresh_neighborhood(displaced);
    }
    if (immediate) refresh_growth_near_vessel(root_capsule);
    return true;
}

std::vector<DirectionId> Simulation3D::feasible_vessel_directions(
    VesselTipSlot slot) const {
    std::vector<DirectionId> feasible;
    if (!vessel_tips_.valid(slot) ||
        vessel_tips_.status(slot) != VesselTipStatus::active) return feasible;
    const Vec3i from = vessel_tips_.position(slot);
    const float diameter = vessel_tips_.diameter_voxels(slot);
    const float grown = vessel_tips_.grown_length_voxels(slot);
    const float maximum = vessel_tips_.max_length_voxels(slot);
    for (DirectionId direction = 1; direction <= 26; ++direction) {
        const Vec3i delta = direction_vector(direction);
        if (config_.thin_layer && delta.z != 0) continue;
        const double length = std::sqrt(static_cast<double>(squared_length(delta)));
        if (static_cast<double>(grown) + length > static_cast<double>(maximum) + 1e-6) continue;
        const Vec3i to = from + delta;
        const std::vector<Vec3i> capsule = rasterize_capsule(from, to, diameter);
        if (!vessel_grid_.all_in_domain(capsule)) continue;
        if (vessel_tips_.role(slot) == VesselBranchRole::outward &&
            std::any_of(capsule.begin(), capsule.end(),
                        [this](Vec3i site) { return grid_.owner(site) != kEmptySlot; })) {
            continue;
        }
        feasible.push_back(direction);
    }
    return feasible;
}

void Simulation3D::schedule_vessel_tip(VesselTipSlot slot) {
    if (!vessel_tips_.valid(slot) ||
        vessel_tips_.status(slot) != VesselTipStatus::active) return;
    const Vec3i position = vessel_tips_.position(slot);
    if (vessel_tips_.role(slot) == VesselBranchRole::inward) {
        const Vec3i bias = vessel_tips_.target(slot) - position;
        vessel_tips_.set_bias_axis(slot, bias);
        if (segment_length(position, vessel_tips_.target(slot)) <=
            config_.angiogenesis.inward_target_tolerance_voxels) {
            vessel_tips_.set_status(slot, VesselTipStatus::reached_target);
            vessel_tips_.set_pending_direction(slot, kStayDirection);
            vessel_tips_.set_next_growth_time(slot, 0.0);
            vessel_tips_.bump_schedule_generation(slot);
            return;
        }
    }
    if (vessel_tips_.grown_length_voxels(slot) >=
        vessel_tips_.max_length_voxels(slot) - 1e-6F) {
        vessel_tips_.set_status(slot, VesselTipStatus::max_length);
        vessel_tips_.set_pending_direction(slot, kStayDirection);
        vessel_tips_.set_next_growth_time(slot, 0.0);
        vessel_tips_.bump_schedule_generation(slot);
        return;
    }

    const std::vector<DirectionId> feasible = feasible_vessel_directions(slot);
    VesselDirectionParameters3D parameters;
    parameters.forward_half_angle_degrees =
        config_.angiogenesis.direction_half_angle_degrees;
    parameters.turn_half_angle_degrees =
        config_.angiogenesis.direction_turn_half_angle_degrees;
    parameters.persistence_probability =
        config_.angiogenesis.direction_persistence_probability;
    parameters.forward_bias = config_.angiogenesis.direction_forward_bias;
    parameters.distance_weight_exponent =
        config_.angiogenesis.direction_distance_weight_exponent;
    const std::uint64_t sequence = vessel_tips_.consume_event_sequence(slot);
    const DirectionId selected = select_vessel_growth_direction(
        feasible, vessel_tips_.bias_axis(slot), vessel_tips_.last_direction(slot),
        parameters, config_.seed, vessel_tips_.uid(slot), sequence);
    if (selected == kStayDirection) {
        vessel_tips_.set_status(slot, VesselTipStatus::blocked);
        vessel_tips_.set_pending_direction(slot, kStayDirection);
        vessel_tips_.set_next_growth_time(slot, 0.0);
        vessel_tips_.bump_schedule_generation(slot);
        return;
    }
    const double duration = std::sqrt(
        static_cast<double>(squared_length(direction_vector(selected)))) /
        vessel_tips_.speed_voxels_per_hour(slot);
    vessel_tips_.set_pending_direction(slot, selected);
    vessel_tips_.set_next_growth_time(slot, clock_.time_hours + duration);
    const std::uint32_t generation = vessel_tips_.bump_schedule_generation(slot);
    schedule(EventKind::vessel_growth, slot, vessel_tips_.uid(slot),
             vessel_tips_.next_growth_time(slot), generation);
}

void Simulation3D::restore_vessel_tip_event(VesselTipSlot slot) {
    if (!vessel_tips_.valid(slot) ||
        vessel_tips_.status(slot) != VesselTipStatus::active) return;
    if (vessel_tips_.pending_direction(slot) == kStayDirection ||
        !(vessel_tips_.next_growth_time(slot) > clock_.time_hours)) {
        throw std::runtime_error(
            "active restored vessel tip has no valid pending event");
    }
    schedule(EventKind::vessel_growth, slot, vessel_tips_.uid(slot),
             vessel_tips_.next_growth_time(slot),
             vessel_tips_.schedule_generation(slot));
}

Simulation3D::VesselGrowthProposal Simulation3D::make_vessel_growth_proposal(
    const Event& event) const {
    VesselGrowthProposal proposal;
    proposal.event = event;
    if (!current(event)) return proposal;
    const VesselTipSlot slot = event.slot;
    proposal.direction = vessel_tips_.pending_direction(slot);
    proposal.from = vessel_tips_.position(slot);
    proposal.to = proposal.from + direction_vector(proposal.direction);
    proposal.capsule = rasterize_capsule(
        proposal.from, proposal.to, vessel_tips_.diameter_voxels(slot));
    if (!vessel_grid_.all_in_domain(proposal.capsule)) return proposal;
    if (vessel_tips_.role(slot) == VesselBranchRole::outward &&
        std::any_of(proposal.capsule.begin(), proposal.capsule.end(),
                    [this](Vec3i site) { return grid_.owner(site) != kEmptySlot; })) {
        return proposal;
    }

    // Consecutive rasterized capsules necessarily overlap at their shared
    // node (and overlap more for thick vessels).  That overlap is part of the
    // same continuous branch and is not an anastomosis.  Any occupied voxel
    // elsewhere in the new capsule is a true side/crossing collision, even
    // when the proposed endpoint is not an existing centreline node.
    const VesselNodeSlot current_node = vessel_tips_.current_node_slot(slot);
    const VesselNodeSlot previous_node = vessel_nodes_.parent_node_slot(current_node);
    const std::vector<Vec3i> allowed_previous_capsule =
        previous_node == kEmptyVesselNodeSlot
            ? rasterize_capsule(proposal.from, proposal.from,
                                vessel_tips_.diameter_voxels(slot))
            : rasterize_capsule(vessel_nodes_.position(previous_node), proposal.from,
                                vessel_tips_.diameter_voxels(slot));
    proposal.anastomosis = std::any_of(
        proposal.capsule.begin(), proposal.capsule.end(),
        [this, &allowed_previous_capsule](Vec3i site) {
            return vessel_grid_.occupied(site) &&
                   std::find(allowed_previous_capsule.begin(),
                             allowed_previous_capsule.end(), site) ==
                       allowed_previous_capsule.end();
        });
    const std::uint64_t time_key = config_.conflict_bucket_hours > 0.0
        ? static_cast<std::uint64_t>(
              std::floor(event.time / config_.conflict_bucket_hours))
        : double_bits(event.time);
    proposal.priority = rng_word(
        config_.seed, event.uid, kVesselConflictEvent, time_key, 0);
    proposal.valid = true;
    return proposal;
}

void Simulation3D::process_vessel_growth(const std::vector<Event>& events) {
    std::vector<VesselGrowthProposal> proposals(events.size());
#ifdef _OPENMP
    const int worker_count = select_worker_count(
        cells_.alive_count(), events.size(), config_, available_worker_threads());
#pragma omp parallel for schedule(static) num_threads(worker_count)
#endif
    for (std::int64_t index = 0;
         index < static_cast<std::int64_t>(events.size()); ++index) {
        proposals[static_cast<std::size_t>(index)] =
            make_vessel_growth_proposal(events[static_cast<std::size_t>(index)]);
    }
    std::sort(proposals.begin(), proposals.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.priority != rhs.priority) return lhs.priority > rhs.priority;
        return lhs.event.uid < rhs.event.uid;
    });

    std::unordered_set<Vec3i, Vec3iHash> reserved;
    for (const VesselGrowthProposal& proposal : proposals) {
        ++stats_.vessel_growth_attempts;
        bool conflict = !proposal.valid;
        if (!conflict) {
            for (const Vec3i site : proposal.capsule) {
                // `reserved` contains voxels that were empty at the start of
                // this same-time batch and claimed by an earlier winner.  A
                // winner has already committed by the time later proposals
                // are inspected, so consulting the now-mutated vessel grid
                // here would incorrectly hide the reservation.
                if (reserved.contains(site)) {
                    conflict = true;
                    break;
                }
            }
        }
        if (conflict) {
            ++stats_.conflict_rejections;
            if (vessel_tips_.valid(proposal.event.slot) &&
                vessel_tips_.status(proposal.event.slot) == VesselTipStatus::active) {
                vessel_tips_.set_pending_direction(proposal.event.slot, kStayDirection);
                vessel_tips_.set_next_growth_time(proposal.event.slot, 0.0);
                schedule_vessel_tip(proposal.event.slot);
            }
            continue;
        }
        for (const Vec3i site : proposal.capsule) {
            if (!vessel_grid_.occupied(site)) reserved.insert(site);
        }
        if (commit_vessel_growth(proposal)) ++stats_.vessel_growth_commits;
    }
}

bool Simulation3D::commit_vessel_growth(const VesselGrowthProposal& proposal) {
    if (!proposal.valid || !current(proposal.event)) return false;
    const VesselTipSlot tip_slot = proposal.event.slot;
    const VesselBranchRole role = vessel_tips_.role(tip_slot);
    if (role == VesselBranchRole::outward &&
        std::any_of(proposal.capsule.begin(), proposal.capsule.end(),
                    [this](Vec3i site) { return grid_.owner(site) != kEmptySlot; })) {
        return false;
    }

    const auto collision = centerline_nodes_.find(proposal.to);
    const bool anastomosis = proposal.anastomosis ||
        (collision != centerline_nodes_.end() &&
         collision->second != vessel_tips_.current_node_slot(tip_slot));
    std::vector<Vec3i> displaced;
    if (role == VesselBranchRole::inward) {
        displaced = displace_cells_for_vessel(proposal.capsule);
    }
    const VesselId vessel_id = vessel_tips_.vessel_id(tip_slot);
    const bool perfused_before = vessel_is_perfused(vessel_id);
    if (!vessel_grid_.add_sites(
            proposal.capsule, role, perfused_before, vessel_id).placed) {
        throw std::logic_error("preflighted vessel segment placement failed");
    }

    const VesselNodeSlot parent_slot = vessel_tips_.current_node_slot(tip_slot);
    VesselNodeInit3D node;
    node.position = proposal.to;
    node.uid = next_vessel_node_uid_++;
    node.parent_uid = vessel_nodes_.uid(parent_slot);
    node.parent_node_slot = parent_slot;
    node.vessel_id = vessel_id;
    node.source_lesion_id = vessel_tips_.source_lesion_id(tip_slot);
    node.role = role;
    node.perfused = perfused_before;
    node.diameter_voxels = vessel_tips_.diameter_voxels(tip_slot);
    node.created_time_hours = clock_.time_hours;
    const VesselNodeSlot node_slot = vessel_nodes_.create(node);
    if (!anastomosis) centerline_nodes_[proposal.to] = node_slot;

    const float length = static_cast<float>(segment_length(proposal.from, proposal.to));
    vessel_tips_.set_position(tip_slot, proposal.to);
    vessel_tips_.set_current_node_uid(tip_slot, node.uid);
    vessel_tips_.set_current_node_slot(tip_slot, node_slot);
    vessel_tips_.set_last_direction(tip_slot, proposal.direction);
    vessel_tips_.set_pending_direction(tip_slot, kStayDirection);
    vessel_tips_.set_next_growth_time(tip_slot, 0.0);
    vessel_tips_.set_grown_length_voxels(
        tip_slot, vessel_tips_.grown_length_voxels(tip_slot) + length);

    if (perfused_before) {
        vascular_influence_.add_sources(proposal.capsule);
        refresh_growth_near_vessel(proposal.capsule);
    }
    if (role == VesselBranchRole::outward && !perfused_before &&
        vessel_tips_.grown_length_voxels(tip_slot) >=
            config_.angiogenesis.outward_external_connection_distance_voxels) {
        activate_vessel_perfusion(vessel_id);
    }

    if (anastomosis) {
        vessel_tips_.set_status(tip_slot, VesselTipStatus::merged);
        ++stats_.vessel_anastomoses;
    } else if (role == VesselBranchRole::inward &&
               segment_length(proposal.to, vessel_tips_.target(tip_slot)) <=
                   config_.angiogenesis.inward_target_tolerance_voxels) {
        vessel_tips_.set_status(tip_slot, VesselTipStatus::reached_target);
    } else if (vessel_tips_.grown_length_voxels(tip_slot) >=
               vessel_tips_.max_length_voxels(tip_slot) - 1e-6F) {
        vessel_tips_.set_status(tip_slot, VesselTipStatus::max_length);
    }

    if (!displaced.empty()) {
        recover_neighborhood(displaced);
        refresh_tumor_surface(displaced);
        refresh_neighborhood(displaced);
    }
    if (vessel_tips_.status(tip_slot) == VesselTipStatus::active) {
        schedule_vessel_tip(tip_slot);
    } else {
        vessel_tips_.bump_schedule_generation(tip_slot);
    }
    return true;
}

bool Simulation3D::vessel_is_perfused(VesselId vessel_id) const {
    return perfused_vessels_.contains(vessel_id);
}

void Simulation3D::activate_vessel_perfusion(VesselId vessel_id) {
    if (!perfused_vessels_.insert(vessel_id).second) return;
    std::vector<Vec3i> affected;
    for (const VesselNodeSlot slot : vessel_nodes_.alive_slots()) {
        if (vessel_nodes_.vessel_id(slot) != vessel_id) continue;
        vessel_nodes_.set_perfused(slot, true);
        std::vector<Vec3i> capsule;
        if (vessel_nodes_.parent_node_slot(slot) == kEmptyVesselNodeSlot) {
            capsule = rasterize_capsule(vessel_nodes_.position(slot),
                                        vessel_nodes_.position(slot),
                                        vessel_nodes_.diameter_voxels(slot));
        } else {
            capsule = rasterize_capsule(
                vessel_nodes_.position(vessel_nodes_.parent_node_slot(slot)),
                vessel_nodes_.position(slot), vessel_nodes_.diameter_voxels(slot));
        }
        vessel_grid_.mark_perfused(capsule);
        vascular_influence_.add_sources(capsule);
        affected.insert(affected.end(), capsule.begin(), capsule.end());
    }
    for (const VesselTipSlot slot : vessel_tips_.alive_slots()) {
        if (vessel_tips_.vessel_id(slot) == vessel_id) vessel_tips_.set_perfused(slot, true);
    }
    std::sort(affected.begin(), affected.end());
    affected.erase(std::unique(affected.begin(), affected.end()), affected.end());
    refresh_growth_near_vessel(affected);
}

std::vector<Vec3i> Simulation3D::displace_cells_for_vessel(
    std::span<const Vec3i> capsule) {
    std::unordered_set<Slot> unique;
    for (const Vec3i site : capsule) {
        for (const Slot slot : grid_.occupants(site)) {
            if (cells_.valid(slot)) unique.insert(slot);
        }
    }
    std::vector<Slot> slots(unique.begin(), unique.end());
    std::sort(slots.begin(), slots.end(), [this](Slot lhs, Slot rhs) {
        return cells_.uid(lhs) < cells_.uid(rhs);
    });
    std::vector<Vec3i> changed;
    for (const Slot slot : slots) {
        if (!cells_.valid(slot)) continue;
        const std::vector<Vec3i> occupied = occupied_sites_for_cell(cells_, slot);
        changed.insert(changed.end(), occupied.begin(), occupied.end());
        if (remove_cell(slot, cells_, grid_, density_)) ++stats_.vascular_displacements;
    }
    std::sort(changed.begin(), changed.end());
    changed.erase(std::unique(changed.begin(), changed.end()), changed.end());
    return changed;
}

void Simulation3D::refresh_growth_near_vessel(std::span<const Vec3i> vessel_sites) {
    if (vessel_sites.empty()) return;

    std::vector<Vec3i> sources(vessel_sites.begin(), vessel_sites.end());
    std::sort(sources.begin(), sources.end());
    sources.erase(std::unique(sources.begin(), sources.end()), sources.end());

    const double cutoff = vascular_influence_.cutoff_radius_voxels();
    const double cutoff_squared = cutoff * cutoff;
    const std::int64_t radius = static_cast<std::int64_t>(std::ceil(cutoff));
    const std::int64_t block_edge = density_.block_edge();
    const std::int64_t minimum_coordinate =
        std::numeric_limits<std::int32_t>::min();
    const std::int64_t maximum_coordinate =
        std::numeric_limits<std::int32_t>::max();

    // Associate every density block with only the vessel sources whose cutoff
    // boxes overlap it. A block is visited once even when many adjacent
    // capsule voxels cover it; the source list then provides the exact
    // Euclidean membership test for each anchor in that block.
    std::unordered_map<Vec3i, std::vector<Vec3i>, Vec3iHash> block_sources;
    for (const Vec3i source : sources) {
        const auto clamp_coordinate = [&](std::int32_t value,
                                          std::int64_t offset) noexcept {
            return std::clamp(static_cast<std::int64_t>(value) + offset,
                              minimum_coordinate, maximum_coordinate);
        };
        const std::int64_t first_x = floor_div_positive(
            clamp_coordinate(source.x, -radius), block_edge);
        const std::int64_t last_x = floor_div_positive(
            clamp_coordinate(source.x, radius), block_edge);
        const std::int64_t first_y = floor_div_positive(
            clamp_coordinate(source.y, -radius), block_edge);
        const std::int64_t last_y = floor_div_positive(
            clamp_coordinate(source.y, radius), block_edge);
        const std::int64_t z_radius = config_.thin_layer ? 0 : radius;
        const std::int64_t first_z = floor_div_positive(
            clamp_coordinate(source.z, -z_radius), block_edge);
        const std::int64_t last_z = floor_div_positive(
            clamp_coordinate(source.z, z_radius), block_edge);
        for (std::int64_t x = first_x; x <= last_x; ++x) {
            for (std::int64_t y = first_y; y <= last_y; ++y) {
                for (std::int64_t z = first_z; z <= last_z; ++z) {
                    block_sources[{static_cast<std::int32_t>(x),
                                   static_cast<std::int32_t>(y),
                                   static_cast<std::int32_t>(z)}]
                        .push_back(source);
                }
            }
        }
    }

    vascular_refresh_diagnostics_.queried_density_blocks += block_sources.size();
    std::unordered_set<Slot> unique;
    for (const auto& [block, relevant_sources] : block_sources) {
        density_.for_each_slot_in_block(block, [&](Slot slot) {
            ++vascular_refresh_diagnostics_.visited_density_slots;
            if (!cells_.valid(slot)) return;
            const Vec3i anchor = cells_.anchor(slot);
            if (!influenced_by_any_source(
                    anchor, relevant_sources, cutoff_squared)) return;
            // A zero maximum relief or an already-truncated influence profile
            // must not cause a biological refresh merely because geometry was
            // inside a coarse query block.
            if (!(vascular_influence_.relief(anchor) > 0.0F)) return;
            unique.insert(slot);
        });
    }
    std::vector<Slot> slots(unique.begin(), unique.end());
    std::sort(slots.begin(), slots.end(), [this](Slot lhs, Slot rhs) {
        return cells_.uid(lhs) < cells_.uid(rhs);
    });
    for (const Slot slot : slots) {
        ++vascular_refresh_diagnostics_.refreshed_cell_slots;
        const GrowthRefreshResult refresh = refresh_growth_state(
            slot, clock_.time_hours, cells_, density_, config_, &vascular_influence_);
        apply_growth_refresh(slot, refresh);
    }
}

std::vector<Slot> Simulation3D::nearby_slots(const std::vector<Vec3i>& sites,
                                              int radius) const {
    std::unordered_set<Slot> unique;
    for (const Vec3i center : sites) {
        const Vec3i delta{radius, radius, config_.thin_layer ? 0 : radius};
        density_.for_each_slot_in_box(center - delta, center + delta, [&](Slot slot) {
            if (cells_.valid(slot)) unique.insert(slot);
        });
    }
    return {unique.begin(), unique.end()};
}

void Simulation3D::recover_neighborhood(std::vector<Vec3i>& changed_sites) {
    struct RecoveryCandidate {
        StageRecoveryProposal proposal;
        CellUid uid{};
        std::uint64_t priority{};
    };
    std::vector<RecoveryCandidate> candidates;
    for (const Slot slot : nearby_slots(changed_sites, 2)) {
        if (!cells_.valid(slot) || cells_.stage(slot) == CellStage::large) continue;
        const CellUid uid = cells_.uid(slot);
        const std::uint64_t sequence = cells_.event_sequence(slot);
        StageRecoveryProposal proposal = make_stage_recovery_proposal(
            slot, cells_, grid_, config_, sequence);
        if (proposal.action == StageRecoveryAction::none) continue;
        candidates.push_back({
            std::move(proposal),
            uid,
            stage_recovery_conflict_priority(config_, clock_.time_hours, uid),
        });
    }
    std::sort(candidates.begin(), candidates.end(), [](const RecoveryCandidate& lhs,
                                                        const RecoveryCandidate& rhs) {
        if (lhs.priority != rhs.priority) return lhs.priority > rhs.priority;
        return lhs.uid < rhs.uid;
    });
    std::unordered_set<Vec3i, Vec3iHash> reserved;
    std::unordered_set<Vec3i, Vec3iHash> locked_groups;
    for (const RecoveryCandidate& candidate : candidates) {
        const StageRecoveryProposal& proposal = candidate.proposal;
        bool conflict = std::any_of(
            proposal.reserved_sites.begin(), proposal.reserved_sites.end(),
            [&reserved](Vec3i site) { return reserved.contains(site); });
        if (!conflict && proposal.locks_colocation_group &&
            locked_groups.contains(proposal.from)) {
            conflict = true;
        }
        if (conflict) {
            ++stats_.conflict_rejections;
            continue;
        }
        reserved.insert(proposal.reserved_sites.begin(),
                        proposal.reserved_sites.end());
        if (proposal.locks_colocation_group) locked_groups.insert(proposal.from);
        if (!commit_stage_recovery_proposal(proposal, cells_, grid_)) {
            ++stats_.conflict_rejections;
            continue;
        }
        cells_.set_event_sequence(
            proposal.slot, proposal.event_sequence + 1);
        density_.move(proposal.from, proposal.target,
                      cells_.type(proposal.slot), proposal.slot);
        changed_sites.push_back(proposal.from);
        changed_sites.push_back(proposal.target);
    }
    std::sort(changed_sites.begin(), changed_sites.end());
    changed_sites.erase(
        std::unique(changed_sites.begin(), changed_sites.end()),
        changed_sites.end());
}

void Simulation3D::refresh_neighborhood(const std::vector<Vec3i>& changed_sites) {
    std::vector<Slot> slots = nearby_slots(
        changed_sites, config_.growth_density_window_edge / 2 + 1);
    std::sort(slots.begin(), slots.end(), [this](Slot lhs, Slot rhs) {
        return cells_.uid(lhs) < cells_.uid(rhs);
    });
    const VascularInfluenceField3D* influence = config_.angiogenesis.enabled
        ? &vascular_influence_ : nullptr;
    for (const Slot slot : slots) {
        const GrowthRefreshResult refresh = refresh_growth_state(
            slot, clock_.time_hours, cells_, density_, config_, influence);
        apply_growth_refresh(slot, refresh);
    }
    refresh_migration_activation_near(changed_sites);
}

Vec3i Simulation3D::migration_activation_query_block(Vec3i anchor) const noexcept {
    const std::int64_t edge = config_.migration_activation_block_edge;
    const auto floor_div = [edge](std::int32_t value) {
        const std::int64_t wide = value;
        std::int64_t quotient = wide / edge;
        if (wide % edge < 0) --quotient;
        return static_cast<std::int32_t>(quotient);
    };
    return {floor_div(anchor.x), floor_div(anchor.y), floor_div(anchor.z)};
}

std::uint8_t Simulation3D::migration_activation_class(Vec3i query_block) const {
    const std::int64_t edge = config_.migration_activation_block_edge;
    const auto representative = [edge](std::int32_t coordinate) {
        return static_cast<std::int32_t>(std::clamp<std::int64_t>(
            static_cast<std::int64_t>(coordinate) * edge,
            std::numeric_limits<std::int32_t>::min(),
            std::numeric_limits<std::int32_t>::max()));
    };
    const Vec3i anchor{representative(query_block.x),
                       representative(query_block.y),
                       representative(query_block.z)};
    const DensityCounts3D counts = density_.estimate_quantized_box(
        anchor, config_.migration_activation_window_edge,
        config_.migration_activation_block_edge, config_.thin_layer);
    const double window_edge = config_.migration_activation_window_edge;
    const double small_capacity = config_.thin_layer
        ? window_edge * window_edge
        : window_edge * window_edge * window_edge;
    const double large_capacity = small_capacity /
        (config_.thin_layer ? 4.0 : 8.0);
    std::uint8_t classes = 0;
    if (static_cast<double>(counts.total()) / small_capacity >=
        config_.migration_activation_threshold) {
        classes |= kSmallMigrationActivationClass;
    }
    if (static_cast<double>(counts.total()) / large_capacity >=
        config_.migration_activation_threshold) {
        classes |= kLargeMigrationActivationClass;
    }
    return classes;
}

void Simulation3D::rebuild_migration_activation_class_cache() {
    migration_activation_class_cache_.clear();
    migration_activation_bulk_slot_visits_ = 0;
    migration_activation_direct_slot_visits_ = 0;
    migration_activation_class_recomputes_ = 0;
    if (!config_.migration_activation_enabled) return;

    std::unordered_set<Vec3i, Vec3iHash> unique;
    for (const Slot slot : cells_.alive_slots()) {
        unique.insert(migration_activation_query_block(cells_.anchor(slot)));
    }
    std::vector<Vec3i> blocks(unique.begin(), unique.end());
    std::sort(blocks.begin(), blocks.end());
    migration_activation_class_cache_.reserve(blocks.size());
    for (const Vec3i block : blocks) {
        migration_activation_class_cache_.emplace(
            block, migration_activation_class(block));
    }
}

void Simulation3D::apply_migration_activation_class(
    Slot slot, std::uint8_t classes) {
    if (!cells_.valid(slot)) return;
    const std::uint8_t stage_class = cells_.stage(slot) == CellStage::large
        ? kLargeMigrationActivationClass
        : kSmallMigrationActivationClass;
    const bool high_density = (classes & stage_class) != 0;
    const bool already_active =
        (cells_.flags(slot) & static_cast<std::uint8_t>(kMigrationActive)) != 0;
    // Density is a trigger, not a continuously clamped gate. Falling density
    // does not terminate an already active interval; only its scheduled end
    // event does. After expiry, a later high-density refresh can activate it
    // again.
    if (!high_density || already_active) return;
    synchronize_migration_schedule(
        slot, refresh_migration_activation_state(
                  slot, clock_.time_hours, cells_, density_, config_));
}

void Simulation3D::refresh_migration_activation_near(
    const std::vector<Vec3i>& changed_sites) {
    if (!config_.migration_activation_enabled || changed_sites.empty()) return;
    const std::int64_t query_edge = config_.migration_activation_block_edge;
    const std::int64_t lower =
        (config_.migration_activation_window_edge - 1) / 2;
    const std::int64_t upper =
        config_.migration_activation_window_edge - lower - 1;
    const auto floor_div = [](std::int64_t value, std::int64_t divisor) {
        std::int64_t quotient = value / divisor;
        if (value % divisor < 0) --quotient;
        return quotient;
    };
    const auto ceil_div = [&](std::int64_t value, std::int64_t divisor) {
        return -floor_div(-value, divisor);
    };
    const auto affected_range = [&](std::int32_t value) {
        const std::int64_t shifted =
            static_cast<std::int64_t>(value) - query_edge / 2;
        return std::pair<std::int64_t, std::int64_t>{
            ceil_div(shifted - upper, query_edge),
            floor_div(shifted + lower, query_edge)};
    };

    // Direct anchors are always refreshed. This is necessary when a cell moves
    // between two blocks whose aggregate classes both remain unchanged, and
    // when a newly born cell is the first occupant of a previously unseen
    // query block.
    std::unordered_set<Slot> direct_slots;
    std::unordered_set<Vec3i, Vec3iHash> direct_blocks;
    for (const Vec3i site : changed_sites) {
        density_.for_each_slot_in_box(site, site, [&](Slot slot) {
            if (!cells_.valid(slot)) return;
            direct_slots.insert(slot);
            direct_blocks.insert(
                migration_activation_query_block(cells_.anchor(slot)));
        });
    }

    // A changed anchor affects only a bounded number of quantized query blocks
    // (typically 2-3 per axis for 70/32). Recompute classes only for blocks that
    // already have occupants/cache state, plus new blocks containing a direct
    // moved/born anchor. Empty surrounding blocks are not materialized.
    std::unordered_set<Vec3i, Vec3iHash> affected_blocks;
    for (const Vec3i site : changed_sites) {
        const auto [first_x, last_x] = affected_range(site.x);
        const auto [first_y, last_y] = affected_range(site.y);
        const auto [first_z, last_z] = affected_range(site.z);
        for (std::int64_t x = first_x; x <= last_x; ++x) {
            for (std::int64_t y = first_y; y <= last_y; ++y) {
                for (std::int64_t z = first_z; z <= last_z; ++z) {
                    if (x < std::numeric_limits<std::int32_t>::min() ||
                        x > std::numeric_limits<std::int32_t>::max() ||
                        y < std::numeric_limits<std::int32_t>::min() ||
                        y > std::numeric_limits<std::int32_t>::max() ||
                        z < std::numeric_limits<std::int32_t>::min() ||
                        z > std::numeric_limits<std::int32_t>::max()) continue;
                    const Vec3i block{static_cast<std::int32_t>(x),
                                      static_cast<std::int32_t>(y),
                                      static_cast<std::int32_t>(z)};
                    if (migration_activation_class_cache_.contains(block) ||
                        direct_blocks.contains(block)) {
                        affected_blocks.insert(block);
                    }
                }
            }
        }
    }
    affected_blocks.insert(direct_blocks.begin(), direct_blocks.end());
    std::vector<Vec3i> blocks(affected_blocks.begin(), affected_blocks.end());
    std::sort(blocks.begin(), blocks.end());

    std::unordered_set<Slot> bulk_slots;
    for (const Vec3i block : blocks) {
        const std::uint8_t new_classes = migration_activation_class(block);
        ++migration_activation_class_recomputes_;
        const auto existing = migration_activation_class_cache_.find(block);
        if (existing == migration_activation_class_cache_.end()) {
            migration_activation_class_cache_.emplace(block, new_classes);
            continue;
        }
        const std::uint8_t changed_classes = existing->second ^ new_classes;
        existing->second = new_classes;
        if (changed_classes == 0) continue;

        const auto block_component = [query_edge](std::int32_t coordinate,
                                                   bool maximum) {
            const std::int64_t value = maximum
                ? (static_cast<std::int64_t>(coordinate) + 1) * query_edge - 1
                : static_cast<std::int64_t>(coordinate) * query_edge;
            return static_cast<std::int32_t>(std::clamp<std::int64_t>(
                value, std::numeric_limits<std::int32_t>::min(),
                std::numeric_limits<std::int32_t>::max()));
        };
        const Vec3i minimum{block_component(block.x, false),
                            block_component(block.y, false),
                            block_component(block.z, false)};
        const Vec3i maximum{block_component(block.x, true),
                            block_component(block.y, true),
                            block_component(block.z, true)};
        density_.for_each_slot_in_box(minimum, maximum, [&](Slot slot) {
            if (!cells_.valid(slot)) return;
            const std::uint8_t stage_class =
                cells_.stage(slot) == CellStage::large
                    ? kLargeMigrationActivationClass
                    : kSmallMigrationActivationClass;
            if ((changed_classes & stage_class) != 0) bulk_slots.insert(slot);
        });
    }

    for (const Slot slot : direct_slots) bulk_slots.erase(slot);
    std::vector<Slot> direct(direct_slots.begin(), direct_slots.end());
    std::vector<Slot> bulk(bulk_slots.begin(), bulk_slots.end());
    const auto by_uid = [this](Slot lhs, Slot rhs) {
        return cells_.uid(lhs) < cells_.uid(rhs);
    };
    std::sort(direct.begin(), direct.end(), by_uid);
    std::sort(bulk.begin(), bulk.end(), by_uid);
    for (const Slot slot : direct) {
        if (!cells_.valid(slot)) continue;
        const auto found = migration_activation_class_cache_.find(
            migration_activation_query_block(cells_.anchor(slot)));
        if (found == migration_activation_class_cache_.end()) {
            throw std::logic_error(
                "migration activation cache is missing a direct cell block");
        }
        ++migration_activation_direct_slot_visits_;
        apply_migration_activation_class(slot, found->second);
    }
    for (const Slot slot : bulk) {
        if (!cells_.valid(slot)) continue;
        const auto found = migration_activation_class_cache_.find(
            migration_activation_query_block(cells_.anchor(slot)));
        if (found == migration_activation_class_cache_.end()) {
            throw std::logic_error(
                "migration activation cache is missing a bulk cell block");
        }
        ++migration_activation_bulk_slot_visits_;
        apply_migration_activation_class(slot, found->second);
    }
}

std::size_t Simulation3D::active_vessel_tip_count() const {
    std::size_t count = 0;
    for (const VesselTipSlot slot : vessel_tips_.alive_slots()) {
        if (vessel_tips_.status(slot) == VesselTipStatus::active) ++count;
    }
    return count;
}

AngiogenesisProcessState3D Simulation3D::angiogenesis_state() const {
    std::vector<LesionAngiogenesisState3D> processes;
    processes.reserve(lesion_angiogenesis_processes_.size());
    for (const auto& [lesion_id, process] : lesion_angiogenesis_processes_) {
        processes.push_back({lesion_id, process.state()});
    }
    return aggregate_angiogenesis_process_states(processes, clock_.time_hours);
}

std::size_t Simulation3D::active_vessel_tip_count(
    LesionId source_lesion_id) const {
    std::size_t count = 0;
    for (const VesselTipSlot slot : vessel_tips_.alive_slots()) {
        if (vessel_tips_.status(slot) == VesselTipStatus::active &&
            current_lesion_for_source(
                vessel_tips_.source_lesion_id(slot)) == source_lesion_id) {
            ++count;
        }
    }
    return count;
}

VasculatureState3D Simulation3D::snapshot_vasculature() const {
    VasculatureState3D state;
    state.lesions.next_lesion_id = lesion_index_.next_lesion_id();
    state.lesions.last_refresh_time_hours = last_lesion_refresh_time_hours_;
    state.lesions.next_refresh_time_hours = next_lesion_refresh_time_hours_;
    state.lesions.refresh_schedule_generation =
        lesion_refresh_schedule_generation_;
    state.lesions.core_identity = lesion_index_.snapshot_core_identity();
    state.lesions.dirty_blocks = lesion_index_.snapshot_dirty_block_state();
    state.lesions.processes.reserve(lesion_angiogenesis_processes_.size());
    for (const auto& [lesion_id, process] : lesion_angiogenesis_processes_) {
        state.lesions.processes.push_back({lesion_id, process.state()});
    }
    std::sort(state.lesions.processes.begin(), state.lesions.processes.end(),
              [](const auto& lhs, const auto& rhs) {
                  return lhs.lesion_id < rhs.lesion_id;
              });
    state.process = aggregate_angiogenesis_process_states(
        state.lesions.processes, clock_.time_hours);
    state.lesions.source_ownership.reserve(
        lesion_source_ownership_.size());
    for (const auto& [source, owner] : lesion_source_ownership_) {
        state.lesions.source_ownership.push_back({source, owner});
    }
    std::sort(state.lesions.source_ownership.begin(),
              state.lesions.source_ownership.end(),
              [](const auto& lhs, const auto& rhs) {
                  return lhs.source_lesion_id < rhs.source_lesion_id;
              });
    for (const VesselNodeSlot slot : vessel_nodes_.alive_slots()) {
        state.nodes.push_back(vessel_nodes_.snapshot(slot));
    }
    for (const VesselTipSlot slot : vessel_tips_.alive_slots()) {
        state.tips.push_back(vessel_tips_.snapshot(slot));
    }
    state.perfused_vessels.assign(perfused_vessels_.begin(), perfused_vessels_.end());
    std::sort(state.perfused_vessels.begin(), state.perfused_vessels.end());
    state.next_vessel_id = next_vessel_id_;
    state.next_node_uid = next_vessel_node_uid_;
    state.next_tip_uid = next_vessel_tip_uid_;
    return state;
}

std::uint64_t Simulation3D::state_checksum() const {
    const std::vector<Slot> slots = snapshot_cell_slots();
    std::uint64_t checksum = splitmix64(slots.size());
    checksum = hash_combine(checksum, cells_.slot_count());
    checksum = hash_combine(checksum, cells_.free_slots().size());
    for (const Slot slot : cells_.free_slots()) {
        checksum = hash_combine(checksum, slot);
    }
    checksum = hash_combine(checksum, next_uid_);
    checksum = hash_combine(checksum, clock_.completed_events);
    checksum = hash_combine(checksum, double_bits(clock_.time_hours));
    checksum = hash_combine(checksum, stats_.migration_attempts);
    checksum = hash_combine(checksum, stats_.migration_commits);
    checksum = hash_combine(checksum, stats_.divisions);
    checksum = hash_combine(checksum, stats_.deaths);
    checksum = hash_combine(checksum, stats_.conflict_rejections);
    checksum = hash_combine(checksum, stats_.angiogenesis_seed_attempts);
    checksum = hash_combine(checksum, stats_.angiogenesis_roots);
    checksum = hash_combine(checksum, stats_.angiogenesis_seed_rejections);
    checksum = hash_combine(checksum, stats_.vessel_growth_attempts);
    checksum = hash_combine(checksum, stats_.vessel_growth_commits);
    checksum = hash_combine(checksum, stats_.vessel_anastomoses);
    checksum = hash_combine(checksum, stats_.vascular_displacements);
    for (const Slot slot : slots) {
        const CellInit cell = cells_.snapshot(slot);
        checksum = hash_combine(checksum, slot);
        checksum = hash_combine(checksum, cell.uid);
        checksum = hash_combine(checksum, cell.parent_uid);
        checksum = hash_combine(checksum, cell.clone_id);
        checksum = hash_combine(checksum, static_cast<std::uint32_t>(cell.anchor.x));
        checksum = hash_combine(checksum, static_cast<std::uint32_t>(cell.anchor.y));
        checksum = hash_combine(checksum, static_cast<std::uint32_t>(cell.anchor.z));
        checksum = hash_combine(checksum, static_cast<std::uint8_t>(cell.type));
        checksum = hash_combine(checksum, static_cast<std::uint8_t>(cell.stage));
        checksum = hash_combine(checksum, cell.viability);
        checksum = hash_combine(checksum, cell.flags);
        checksum = hash_combine(checksum, cell.last_direction);
        checksum = hash_combine(checksum, float_bits(cell.inherent_growth_rate));
        checksum = hash_combine(checksum, float_bits(cell.density_growth_rate));
        checksum = hash_combine(checksum, float_bits(cell.migration_rate));
        checksum = hash_combine(checksum, float_bits(cell.normal_migration_rate));
        checksum = hash_combine(
            checksum, double_bits(cell.migration_activation_end_time));
        checksum = hash_combine(checksum, float_bits(cell.division_work_remaining));
        checksum = hash_combine(checksum, double_bits(cell.next_migration_time));
        checksum = hash_combine(checksum, double_bits(cell.next_division_time));
        checksum = hash_combine(checksum, double_bits(cell.death_deadline));
        checksum = hash_combine(checksum, double_bits(cell.last_update_time));
        checksum = hash_combine(checksum, cell.event_sequence);
        checksum = hash_combine(checksum, cell.migration_schedule_generation);
        checksum = hash_combine(checksum, cell.division_schedule_generation);
        checksum = hash_combine(checksum, cell.death_schedule_generation);
    }
    for (const LineageEdge& edge : lineage_) {
        checksum = hash_combine(checksum, double_bits(edge.birth_time));
        checksum = hash_combine(checksum, edge.child_uid);
        checksum = hash_combine(checksum, edge.parent_uid);
        checksum = hash_combine(checksum, edge.clone_id);
        checksum = hash_combine(checksum, static_cast<std::uint8_t>(edge.type));
    }

    const VasculatureState3D vascular = snapshot_vasculature();
    checksum = hash_combine(checksum, vascular.next_vessel_id);
    checksum = hash_combine(checksum, vascular.next_node_uid);
    checksum = hash_combine(checksum, vascular.next_tip_uid);
    checksum = hash_combine(checksum, vascular.process.eligible);
    checksum = hash_combine(checksum, double_bits(vascular.process.next_seed_time_hours));
    checksum = hash_combine(checksum, double_bits(vascular.process.eligibility_started_hours));
    checksum = hash_combine(checksum, double_bits(vascular.process.accumulated_eligible_hours));
    checksum = hash_combine(checksum, vascular.process.event_sequence);
    checksum = hash_combine(checksum, vascular.process.schedule_generation);
    checksum = hash_combine(checksum, vascular.process.attempted_events);
    checksum = hash_combine(checksum, vascular.process.committed_roots);
    checksum = hash_combine(checksum, vascular.process.rejected_events);
    checksum = hash_combine(checksum, vascular.lesions.next_lesion_id);
    checksum = hash_combine(
        checksum, double_bits(vascular.lesions.last_refresh_time_hours));
    checksum = hash_combine(
        checksum, double_bits(vascular.lesions.next_refresh_time_hours));
    checksum = hash_combine(
        checksum, vascular.lesions.refresh_schedule_generation);
    for (const LesionCoreIdentity3D& core : vascular.lesions.core_identity) {
        checksum = hash_combine(checksum,
                                static_cast<std::uint32_t>(core.block.x));
        checksum = hash_combine(checksum,
                                static_cast<std::uint32_t>(core.block.y));
        checksum = hash_combine(checksum,
                                static_cast<std::uint32_t>(core.block.z));
        checksum = hash_combine(checksum, core.lesion_id);
    }
    for (const LesionDirtyBlockState3D& block :
         vascular.lesions.dirty_blocks) {
        checksum = hash_combine(checksum,
                                static_cast<std::uint32_t>(block.block.x));
        checksum = hash_combine(checksum,
                                static_cast<std::uint32_t>(block.block.y));
        checksum = hash_combine(checksum,
                                static_cast<std::uint32_t>(block.block.z));
        checksum = hash_combine(checksum, block.exists);
        checksum = hash_combine(checksum, block.cell_count);
        checksum = hash_combine(checksum, block.occupied_voxel_count);
        checksum = hash_combine(checksum,
                                double_bits(block.biological_volume));
        checksum = hash_combine(
            checksum, static_cast<std::uint64_t>(block.cell_coordinate_sum_x));
        checksum = hash_combine(
            checksum, static_cast<std::uint64_t>(block.cell_coordinate_sum_y));
        checksum = hash_combine(
            checksum, static_cast<std::uint64_t>(block.cell_coordinate_sum_z));
        checksum = hash_combine(
            checksum,
            static_cast<std::uint64_t>(block.occupied_coordinate_sum_x));
        checksum = hash_combine(
            checksum,
            static_cast<std::uint64_t>(block.occupied_coordinate_sum_y));
        checksum = hash_combine(
            checksum,
            static_cast<std::uint64_t>(block.occupied_coordinate_sum_z));
    }
    for (const LesionAngiogenesisState3D& lesion :
         vascular.lesions.processes) {
        const AngiogenesisProcessState3D& process = lesion.process;
        checksum = hash_combine(checksum, lesion.lesion_id);
        checksum = hash_combine(checksum, process.eligible);
        checksum = hash_combine(
            checksum, double_bits(process.next_seed_time_hours));
        checksum = hash_combine(
            checksum, double_bits(process.eligibility_started_hours));
        checksum = hash_combine(
            checksum, double_bits(process.accumulated_eligible_hours));
        checksum = hash_combine(checksum, process.event_sequence);
        checksum = hash_combine(checksum, process.schedule_generation);
        checksum = hash_combine(checksum, process.attempted_events);
        checksum = hash_combine(checksum, process.committed_roots);
        checksum = hash_combine(checksum, process.rejected_events);
    }
    for (const LesionSourceOwnership3D& ownership :
         vascular.lesions.source_ownership) {
        checksum = hash_combine(checksum, ownership.source_lesion_id);
        checksum = hash_combine(checksum, ownership.current_lesion_id);
    }
    for (const VesselNodeInit3D& node : vascular.nodes) {
        checksum = hash_combine(checksum, node.uid);
        checksum = hash_combine(checksum, node.parent_uid);
        checksum = hash_combine(checksum, node.parent_node_slot);
        checksum = hash_combine(checksum, node.vessel_id);
        checksum = hash_combine(checksum, node.source_lesion_id);
        checksum = hash_combine(checksum, static_cast<std::uint32_t>(node.position.x));
        checksum = hash_combine(checksum, static_cast<std::uint32_t>(node.position.y));
        checksum = hash_combine(checksum, static_cast<std::uint32_t>(node.position.z));
        checksum = hash_combine(checksum, static_cast<std::uint8_t>(node.role));
        checksum = hash_combine(checksum, node.perfused);
        checksum = hash_combine(checksum, float_bits(node.diameter_voxels));
        checksum = hash_combine(checksum, double_bits(node.created_time_hours));
    }
    for (const VesselTipInit3D& tip : vascular.tips) {
        checksum = hash_combine(checksum, tip.uid);
        checksum = hash_combine(checksum, tip.vessel_id);
        checksum = hash_combine(checksum, tip.source_lesion_id);
        checksum = hash_combine(checksum, tip.current_node_uid);
        checksum = hash_combine(checksum, tip.current_node_slot);
        checksum = hash_combine(checksum, static_cast<std::uint32_t>(tip.position.x));
        checksum = hash_combine(checksum, static_cast<std::uint32_t>(tip.position.y));
        checksum = hash_combine(checksum, static_cast<std::uint32_t>(tip.position.z));
        checksum = hash_combine(checksum, static_cast<std::uint32_t>(tip.bias_axis.x));
        checksum = hash_combine(checksum, static_cast<std::uint32_t>(tip.bias_axis.y));
        checksum = hash_combine(checksum, static_cast<std::uint32_t>(tip.bias_axis.z));
        checksum = hash_combine(checksum, static_cast<std::uint32_t>(tip.target.x));
        checksum = hash_combine(checksum, static_cast<std::uint32_t>(tip.target.y));
        checksum = hash_combine(checksum, static_cast<std::uint32_t>(tip.target.z));
        checksum = hash_combine(checksum, static_cast<std::uint8_t>(tip.role));
        checksum = hash_combine(checksum, static_cast<std::uint8_t>(tip.status));
        checksum = hash_combine(checksum, tip.perfused);
        checksum = hash_combine(checksum, tip.last_direction);
        checksum = hash_combine(checksum, tip.pending_direction);
        checksum = hash_combine(checksum, float_bits(tip.diameter_voxels));
        checksum = hash_combine(checksum, float_bits(tip.speed_voxels_per_hour));
        checksum = hash_combine(checksum, float_bits(tip.max_length_voxels));
        checksum = hash_combine(checksum, float_bits(tip.grown_length_voxels));
        checksum = hash_combine(checksum, double_bits(tip.next_growth_time));
        checksum = hash_combine(checksum, tip.event_sequence);
        checksum = hash_combine(checksum, tip.schedule_generation);
    }
    for (const VesselId vessel_id : vascular.perfused_vessels) {
        checksum = hash_combine(checksum, vessel_id);
    }
    return checksum;
}

std::vector<CellInit> Simulation3D::snapshot_cells() const {
    std::vector<CellInit> result;
    result.reserve(cells_.alive_count());
    for (const Slot slot : snapshot_cell_slots()) result.push_back(cells_.snapshot(slot));
    return result;
}

std::vector<Slot> Simulation3D::snapshot_cell_slots() const {
    std::vector<Slot> slots = cells_.alive_slots();
    std::sort(slots.begin(), slots.end(), [this](Slot lhs, Slot rhs) {
        return cells_.uid(lhs) < cells_.uid(rhs);
    });
    return slots;
}

}  // namespace atcg3d
