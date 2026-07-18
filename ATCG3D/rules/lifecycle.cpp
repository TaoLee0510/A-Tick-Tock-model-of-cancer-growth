#include "rules/lifecycle.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_set>

#include "core/stateless_rng.hpp"
#include "geometry/directions.hpp"
#include "geometry/footprint.hpp"
#include "rules/density.hpp"
#include "rules/initial_rates.hpp"
#include "rules/migration.hpp"

namespace atcg3d {
namespace {

constexpr float kDivisionWorkTolerance = 1.0e-5F;
constexpr std::uint64_t kNormalMigrationRateDomain = 0x4e4f524d4d494752ULL;
constexpr std::uint64_t kActivationDurationDomain = 0x4143544455524154ULL;
constexpr std::uint64_t kDivisionMigrationRateDomain = 0x4449564d49475241ULL;

bool same_schedule_time(double lhs, double rhs) noexcept {
    return std::abs(lhs - rhs) <=
           1.0e-10 * std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

double float_time_tolerance(double lhs, double rhs) noexcept {
    return 2.0 * static_cast<double>(std::numeric_limits<float>::epsilon()) *
           std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

double representable_activation_end(double now, double sampled_end,
                                    double latest_end) {
    if (!std::isfinite(now) || !std::isfinite(sampled_end) ||
        !std::isfinite(latest_end) || latest_end <= now) {
        return 0.0;
    }

    // Find the first stored float strictly after now and the last stored
    // float not beyond the unfinished division cycle. If the interval is
    // shorter than one float ULP, activation is safely skipped.
    float earliest = static_cast<float>(now);
    if (!std::isfinite(earliest)) return 0.0;
    if (static_cast<double>(earliest) <= now) {
        earliest = std::nextafter(
            earliest, std::numeric_limits<float>::infinity());
    }
    float latest = static_cast<float>(latest_end);
    if (!std::isfinite(latest)) return 0.0;
    if (static_cast<double>(latest) > latest_end) {
        latest = std::nextafter(
            latest, -std::numeric_limits<float>::infinity());
    }
    if (!std::isfinite(earliest) || earliest > latest) return 0.0;

    float stored = static_cast<float>(sampled_end);
    if (!std::isfinite(stored)) return 0.0;
    if (static_cast<double>(stored) < sampled_end) {
        stored = std::nextafter(
            stored, std::numeric_limits<float>::infinity());
    }
    stored = std::clamp(stored, earliest, latest);
    return static_cast<double>(stored);
}

CellInit daughter_init(const CellStore3D& cells,
                       Slot mother,
                       CellUid uid,
                       Vec3i anchor,
                       CellStage stage,
                       double now,
                       float inherited_growth_rate) {
    CellInit daughter;
    daughter.anchor = anchor;
    daughter.uid = uid;
    daughter.parent_uid = cells.uid(mother);
    daughter.clone_id = cells.clone_id(mother);
    daughter.type = cells.type(mother);
    daughter.stage = stage;
    daughter.viability = 1;
    daughter.flags = kDirtyDensity;
    daughter.last_direction = kStayDirection;
    daughter.inherent_growth_rate = inherited_growth_rate;
    daughter.density_growth_rate = inherited_growth_rate;
    daughter.migration_rate = cells.migration_rate(mother);
    daughter.normal_migration_rate = cells.normal_migration_rate(mother);
    daughter.last_update_time = now;
    daughter.event_sequence = 0;
    return daughter;
}

float sample_cycle_inherent_migration_rate(CellType type,
                                           const Model3DConfig& config,
                                           CellUid uid,
                                           std::uint64_t event_sequence) {
    if (type == CellType::K && config.initial_K_migration_rate_model == "fixed") {
        return static_cast<float>(config.initial_K_migration_rate);
    }
    const BetaRateConfig& parameters = type == CellType::r
        ? config.activated_r_migration_beta : config.initial_K_migration_beta;
    double rate = parameters.scale * sample_stateless_beta_fraction(
        parameters.alpha, parameters.beta, config.seed, uid,
        kDivisionMigrationRateDomain + (type == CellType::r ? 0U : 16U),
        event_sequence);
    if (parameters.lower_clamp_enabled && rate <= parameters.lower_clamp_threshold) {
        rate = parameters.lower_clamp_value;
    }
    return static_cast<float>(rate);
}

void reset_migration_for_new_cycle(Slot slot,
                                   CellStore3D& cells,
                                   const Model3DConfig& config,
                                   std::uint64_t event_sequence) {
    const float inherent = sample_cycle_inherent_migration_rate(
        cells.type(slot), config, cells.uid(slot), event_sequence);
    cells.set_migration_rate(slot, inherent);
    cells.set_normal_migration_rate(
        slot, sample_normal_migration_rate(
                  cells.type(slot), inherent, config, cells.uid(slot),
                  event_sequence));
    cells.set_migration_activation_end_time(slot, 0.0);
    cells.set_flags(slot, cells.flags(slot) &
        static_cast<std::uint8_t>(~kMigrationActive));
    cells.set_last_direction(slot, kStayDirection);
}

void update_colocation_stages(Vec3i site,
                              CellStore3D& cells,
                              const SparseChunkGrid3D& grid) {
    const std::vector<Slot> occupants = grid.occupants(site);
    for (const Slot slot : occupants) {
        if (cells.valid(slot)) {
            cells.set_stage(slot, occupants.size() <= 1 ? CellStage::small : CellStage::ultrasmall);
        }
    }
}

std::vector<DirectionId> free_small_directions(Vec3i anchor,
                                                const SparseChunkGrid3D& grid,
                                                bool thin_layer) {
    std::vector<DirectionId> result;
    for (DirectionId direction = 1; direction <= 26; ++direction) {
        const Vec3i delta = direction_vector(direction);
        if ((!thin_layer || delta.z == 0) && grid.available(anchor + delta)) {
            result.push_back(direction);
        }
    }
    return result;
}

}  // namespace

namespace {

std::uint64_t conflict_time_key(const Model3DConfig& config,
                                double event_time) noexcept {
    return config.conflict_bucket_hours > 0.0
        ? static_cast<std::uint64_t>(
              std::floor(event_time / config.conflict_bucket_hours))
        : std::bit_cast<std::uint64_t>(event_time);
}

std::uint64_t cell_conflict_priority(const Model3DConfig& config,
                                     double event_time,
                                     CellUid uid,
                                     RngEventKind event_kind) noexcept {
    return rng_word(config.seed, uid, static_cast<std::uint64_t>(event_kind),
                    conflict_time_key(config, event_time), 0);
}

}  // namespace

std::uint64_t division_conflict_priority(const Model3DConfig& config,
                                         double event_time,
                                         CellUid uid) {
    return cell_conflict_priority(
        config, event_time, uid, RngEventKind::division_conflict_priority);
}

std::uint64_t stage_recovery_conflict_priority(
    const Model3DConfig& config,
    double event_time,
    CellUid uid) {
    return cell_conflict_priority(
        config, event_time, uid, RngEventKind::stage_recovery_conflict_priority);
}

double sample_division_delay(double density_growth_rate,
                             std::uint64_t seed,
                             CellUid uid,
                             std::uint64_t event_sequence) {
    return sample_division_delay(density_growth_rate, DivisionTimingConfig{},
                                 seed, uid, event_sequence);
}

double sample_division_delay(double density_growth_rate,
                             const DivisionTimingConfig& timing,
                             std::uint64_t seed,
                             CellUid uid,
                             std::uint64_t event_sequence) {
    if (!(density_growth_rate > 0.0)) {
        return 0.0;
    }
    const double expected = timing.base_cycle_hours / density_growth_rate;
    const double minimum = timing.minimum_fraction * expected;
    const double range = timing.stochastic_tail_fraction * expected;
    const double probability = std::clamp(
        timing.stochastic_time_quantum_hours /
            std::max(range, timing.stochastic_time_quantum_hours),
        1e-12, 1.0);
    if (probability >= 1.0) return minimum + timing.stochastic_time_quantum_hours;
    const double uniform = std::clamp(
        rng_unit(seed, uid, static_cast<std::uint64_t>(RngEventKind::division_timing), event_sequence),
        1e-15, 1.0 - 1e-15);
    const double geometric =
        1.0 + std::floor(std::log1p(-uniform) / std::log1p(-probability));
    return minimum + geometric * timing.stochastic_time_quantum_hours;
}

double sample_death_delay(double mean_hours,
                          std::uint64_t seed,
                          CellUid uid,
                          std::uint64_t event_sequence) {
    if (!(mean_hours > 0.0) || !std::isfinite(mean_hours)) return 0.0;
    const double probability = std::clamp(1.0 / mean_hours, 1e-12, 1.0);
    if (probability >= 1.0) return 1.0;
    const double uniform = std::clamp(
        rng_unit(seed, uid, static_cast<std::uint64_t>(RngEventKind::death_timing),
                 event_sequence),
        1e-15, 1.0 - 1e-15);
    return 1.0 + std::floor(
        std::log1p(-uniform) / std::log1p(-probability));
}

float sample_normal_migration_rate(CellType type,
                                   float inherent_rate,
                                   const Model3DConfig& config,
                                   CellUid uid,
                                   std::uint64_t event_sequence) {
    if (type == CellType::K) return inherent_rate;
    return static_cast<float>(config.normal_r_migration_beta.scale *
        sample_stateless_beta_fraction(
            config.normal_r_migration_beta.alpha,
            config.normal_r_migration_beta.beta,
            config.seed, uid, kNormalMigrationRateDomain, event_sequence));
}

double effective_migration_rate(Slot slot,
                                const CellStore3D& cells,
                                const Model3DConfig&) {
    if (!cells.valid(slot)) return 0.0;
    const bool active =
        (cells.flags(slot) & static_cast<std::uint8_t>(kMigrationActive)) != 0;
    if (active) return cells.migration_rate(slot);
    // Zero is retained as a compatibility sentinel for manually constructed
    // test/benchmark cells created before normal_migration_rate was introduced.
    return cells.normal_migration_rate(slot) > 0.0F
        ? cells.normal_migration_rate(slot) : cells.migration_rate(slot);
}

bool refresh_migration_activation_state(Slot slot,
                                        double now,
                                        CellStore3D& cells,
                                        const BlockDensityIndex3D& density,
                                        const Model3DConfig& config) {
    if (!cells.valid(slot) || !config.migration_activation_enabled ||
        (cells.flags(slot) & static_cast<std::uint8_t>(kMigrationActive)) != 0) {
        return false;
    }
    const double migration_density = migration_activation_density(
        density, cells.anchor(slot), cells.stage(slot),
        config.migration_activation_window_edge,
        config.migration_activation_block_edge, config.thin_layer);
    if (migration_density < config.migration_activation_threshold) return false;
    const double density_rate = cells.density_growth_rate(slot);
    if (!(density_rate > config.death_growth_rate_threshold)) return false;
    const float quantized_now = static_cast<float>(now);
    if (!std::isfinite(quantized_now)) return false;
    const double elapsed_since_growth_refresh = std::max(
        0.0, static_cast<double>(quantized_now) -
                 cells.last_update_time(slot));
    const double remaining_work = std::max(
        0.0, static_cast<double>(cells.division_work_remaining(slot)) -
                 density_rate * elapsed_since_growth_refresh);
    const double remaining_cycle_hours = remaining_work / density_rate;
    if (!(remaining_cycle_hours > 0.0) || !std::isfinite(remaining_cycle_hours)) {
        return false;
    }
    const double latest_end = now + remaining_cycle_hours;
    // Check that a strictly-future stored time exists before consuming an RNG
    // sequence. Very short intervals can be smaller than one float ULP.
    if (representable_activation_end(now, now, latest_end) == 0.0) {
        return false;
    }
    const std::uint64_t sequence = cells.consume_event_sequence(slot);
    const double duration_fraction = sample_stateless_beta_fraction(
        config.migration_activation_duration_alpha,
        config.migration_activation_duration_beta,
        config.seed, cells.uid(slot), kActivationDurationDomain, sequence);
    const double sampled_end = now + duration_fraction * remaining_cycle_hours;
    const double end_time = representable_activation_end(
        now, sampled_end, latest_end);
    if (!(end_time > now) || end_time > latest_end) {
        throw std::runtime_error("sampled migration activation duration is invalid");
    }
    cells.set_flags(slot, cells.flags(slot) |
        static_cast<std::uint8_t>(kMigrationActive));
    const double stored_end =
        cells.set_migration_activation_end_time(slot, end_time);
    if (!(stored_end > now) || stored_end > latest_end) {
        throw std::logic_error(
            "migration activation end changed outside its representable interval");
    }
    cells.set_last_direction(slot, kStayDirection);
    return true;
}

bool expire_migration_activation_state(Slot slot,
                                       double now,
                                       CellStore3D& cells,
                                       const Model3DConfig& config) {
    if (!cells.valid(slot) ||
        (cells.flags(slot) & static_cast<std::uint8_t>(kMigrationActive)) == 0) {
        return false;
    }
    const double end_time = cells.migration_activation_end_time(slot);
    if (!(end_time > 0.0) || now + 1.0e-10 < end_time) return false;
    const std::uint64_t sequence = cells.consume_event_sequence(slot);
    cells.set_normal_migration_rate(
        slot, sample_normal_migration_rate(
                  cells.type(slot), cells.migration_rate(slot), config,
                  cells.uid(slot), sequence));
    cells.set_flags(slot, cells.flags(slot) &
        static_cast<std::uint8_t>(~kMigrationActive));
    cells.set_migration_activation_end_time(slot, 0.0);
    cells.set_last_direction(slot, kStayDirection);
    return true;
}

bool migration_allowed_for_cell(Slot slot,
                                const CellStore3D& cells,
                                const Model3DConfig& config) {
    (void)config;
    return cells.valid(slot);
}

double r_to_K_division_density(const BlockDensityIndex3D& density,
                               Vec3i anchor,
                               CellStage stage,
                               const RToKConversionConfig& config,
                               bool thin_layer) {
    return migration_activation_density(
        density, anchor, stage, config.density_window_edge,
        config.query_block_edge, thin_layer);
}

bool should_convert_r_daughter(CellType mother_type,
                               double local_density,
                               const RToKConversionConfig& config,
                               std::uint64_t seed,
                               CellUid mother_uid,
                               std::uint64_t division_event_sequence) {
    if (!config.enabled || mother_type != CellType::r ||
        !std::isfinite(local_density) ||
        local_density < config.density_threshold ||
        config.probability_per_division <= 0.0) {
        return false;
    }
    if (config.probability_per_division >= 1.0) {
        return true;
    }
    return rng_unit(
               seed, mother_uid,
               static_cast<std::uint64_t>(RngEventKind::division_type_conversion),
               division_event_sequence) < config.probability_per_division;
}

GrowthRefreshResult refresh_growth_state(
    Slot slot,
    double now,
    CellStore3D& cells,
    const BlockDensityIndex3D& density,
    const Model3DConfig& config,
    const VascularInfluenceField3D* vascular_influence) {
    GrowthRefreshResult result;
    if (!cells.valid(slot)) {
        return result;
    }
    const double old_update_time = cells.last_update_time(slot);
    if (!std::isfinite(now) || now < 0.0 ||
        now + float_time_tolerance(now, old_update_time) < old_update_time) {
        throw std::invalid_argument("growth refresh time is invalid or moved backwards");
    }
    const double old_division_time = cells.next_division_time(slot);
    const double old_death_time = cells.death_deadline(slot);
    const float old_density_rate = cells.density_growth_rate(slot);
    const float old_work = cells.division_work_remaining(slot);
    if (!std::isfinite(old_work) || old_work < 0.0F) {
        throw std::logic_error("cell division work is invalid");
    }
    // Work against quantized observation times, so intermediate refresh errors
    // telescope instead of accumulating one rounding bias per refresh.
    const double stored_now = cells.set_last_update_time(slot, now);
    const double elapsed = std::max(0.0, stored_now - old_update_time);
    const double completed_work =
        std::max(0.0, static_cast<double>(old_density_rate)) * elapsed;
    const float remaining = static_cast<float>(std::max(
        0.0, static_cast<double>(old_work) - completed_work));
    cells.set_division_work_remaining(
        slot, remaining <= kDivisionWorkTolerance ? 0.0F : remaining);

    const double rate = density_growth_rate_for_cell(
        cells, slot, density, config, vascular_influence);
    if (!std::isfinite(rate)) {
        throw std::runtime_error("density growth rate is not finite");
    }
    const float quantized_rate = static_cast<float>(rate);
    cells.set_density_growth_rate(slot, quantized_rate);
    result.migration_activation_changed =
        refresh_migration_activation_state(slot, now, cells, density, config);

    const double effective_rate = static_cast<double>(quantized_rate);
    if (effective_rate > config.death_growth_rate_threshold) {
        cells.set_death_deadline(slot, 0.0);
        const float work = cells.division_work_remaining(slot);
        if (work > kDivisionWorkTolerance) {
            if (!(old_density_rate == quantized_rate &&
                  old_division_time > now &&
                  static_cast<double>(old_density_rate) >
                      config.death_growth_rate_threshold)) {
                cells.set_next_division_time(
                    slot, now + static_cast<double>(work) / effective_rate);
            }
        } else if (old_work > kDivisionWorkTolerance || old_division_time <= now) {
            // Work became complete during this lazy update. Equal-time events
            // are legal and are handled on the next scheduler turn if this
            // refresh was triggered by another event in the current batch.
            cells.set_next_division_time(slot, now);
        }
    } else {
        cells.set_next_division_time(slot, 0.0);
        if (cells.death_deadline(slot) <= now) {
            const double delay = cells.type(slot) == CellType::r
                ? config.r_death_delay_hours : config.K_death_delay_hours;
            const std::uint64_t sequence = cells.consume_event_sequence(slot);
            cells.set_death_deadline(
                slot, now + sample_death_delay(delay, config.seed,
                                               cells.uid(slot), sequence));
        }
    }
    result.division_time_changed =
        !same_schedule_time(old_division_time, cells.next_division_time(slot));
    result.death_time_changed =
        !same_schedule_time(old_death_time, cells.death_deadline(slot));
    return result;
}

void initialize_division_cycle(Slot slot,
                               double now,
                               CellStore3D& cells,
                               const Model3DConfig& config) {
    if (!cells.valid(slot)) return;
    if (!std::isfinite(now) || now < 0.0) {
        throw std::invalid_argument("division cycle start time must be finite and nonnegative");
    }
    const double reference_rate =
        std::max(1.0e-12, static_cast<double>(cells.inherent_growth_rate(slot)));
    const std::uint64_t sequence = cells.consume_event_sequence(slot);
    const double delay = sample_division_delay(
        reference_rate, config.division_timing, config.seed, cells.uid(slot), sequence);
    const double work = reference_rate * delay;
    if (!(work > 0.0) || !std::isfinite(work) ||
        work > std::numeric_limits<float>::max()) {
        throw std::runtime_error("sampled division work is invalid");
    }
    const float stored_work = static_cast<float>(work);
    cells.set_division_work_remaining(slot, stored_work);
    cells.set_last_update_time(slot, now);
    const double density_rate = cells.density_growth_rate(slot);
    cells.set_next_division_time(
        slot, density_rate > config.death_growth_rate_threshold
                  ? now + static_cast<double>(stored_work) / density_rate
                  : 0.0);
}

bool remove_cell(Slot slot,
                 CellStore3D& cells,
                 SparseChunkGrid3D& grid,
                 BlockDensityIndex3D& density) {
    if (!cells.valid(slot)) {
        return false;
    }
    const Vec3i anchor = cells.anchor(slot);
    const CellStage stage = cells.stage(slot);
    const CellType type = cells.type(slot);
    if (stage == CellStage::large) {
        grid.remove_large(anchor, slot);
    } else {
        grid.remove(anchor, slot);
    }
    density.remove(anchor, type, slot);
    cells.erase(slot);
    if (stage != CellStage::large) {
        update_colocation_stages(anchor, cells, grid);
    }
    return true;
}

StageRecoveryProposal make_stage_recovery_proposal(
    Slot slot,
    const CellStore3D& cells,
    const SparseChunkGrid3D& grid,
    const Model3DConfig& config,
    std::uint64_t event_sequence) {
    StageRecoveryProposal proposal;
    if (!cells.valid(slot) || cells.stage(slot) == CellStage::large) {
        return proposal;
    }
    const Vec3i site = cells.anchor(slot);
    proposal.slot = slot;
    proposal.uid = cells.uid(slot);
    proposal.event_sequence = event_sequence;
    proposal.from = site;
    const std::vector<Slot> group = grid.occupants(site);
    if (group.size() > 1) {
        const auto free_directions = free_small_directions(site, grid, config.thin_layer);
        if (free_directions.empty()) {
            return proposal;
        }
        const auto index = rng_bounded(config.seed, proposal.uid,
                                       static_cast<std::uint64_t>(RngEventKind::division_location),
                                       event_sequence, free_directions.size());
        proposal.target = site + direction_vector(
            free_directions[static_cast<std::size_t>(index)]);
        proposal.action = StageRecoveryAction::separate_colocated;
        proposal.reserved_sites = {proposal.target};
        proposal.locks_colocation_group = true;
        return proposal;
    }

    std::vector<Vec3i> candidates;
    for (const Vec3i anchor : stage_recovery_anchors(site)) {
        bool available = true;
        for (const Vec3i voxel : large_footprint(anchor)) {
            const Slot owner = grid.owner(voxel);
            if (owner != slot && !grid.available(voxel)) {
                available = false;
                break;
            }
        }
        if (available) {
            candidates.push_back(anchor);
        }
    }
    if (candidates.empty()) {
        return proposal;
    }
    const auto index = rng_bounded(config.seed, proposal.uid,
                                   static_cast<std::uint64_t>(RngEventKind::division_location),
                                   event_sequence, candidates.size());
    proposal.target = candidates[static_cast<std::size_t>(index)];
    proposal.action = StageRecoveryAction::restore_large;
    for (const Vec3i voxel : large_footprint(proposal.target)) {
        if (voxel != site) proposal.reserved_sites.push_back(voxel);
    }
    return proposal;
}

bool commit_stage_recovery_proposal(const StageRecoveryProposal& proposal,
                                    CellStore3D& cells,
                                    SparseChunkGrid3D& grid) {
    if (proposal.action == StageRecoveryAction::none ||
        !cells.valid(proposal.slot) || cells.uid(proposal.slot) != proposal.uid ||
        cells.anchor(proposal.slot) != proposal.from ||
        cells.stage(proposal.slot) == CellStage::large) {
        return false;
    }
    if (proposal.action == StageRecoveryAction::separate_colocated) {
        if (!grid.available(proposal.target) ||
            !grid.remove(proposal.from, proposal.slot)) {
            return false;
        }
        if (!grid.place_single(proposal.target, proposal.slot)) {
            grid.add_colocated(proposal.from, proposal.slot);
            return false;
        }
        cells.set_anchor(proposal.slot, proposal.target);
        cells.set_stage(proposal.slot, CellStage::small);
        update_colocation_stages(proposal.from, cells, grid);
        return true;
    }
    for (const Vec3i voxel : large_footprint(proposal.target)) {
        const Slot owner = grid.owner(voxel);
        if (owner != proposal.slot && !grid.available(voxel)) return false;
    }
    if (!grid.remove(proposal.from, proposal.slot)) return false;
    if (!grid.place_large(proposal.target, proposal.slot)) {
        grid.place_single(proposal.from, proposal.slot);
        return false;
    }
    cells.set_anchor(proposal.slot, proposal.target);
    cells.set_stage(proposal.slot, CellStage::large);
    return true;
}

bool try_stage_recovery(Slot slot,
                        CellStore3D& cells,
                        SparseChunkGrid3D& grid,
                        const Model3DConfig& config,
                        std::uint64_t event_sequence) {
    return commit_stage_recovery_proposal(
        make_stage_recovery_proposal(
            slot, cells, grid, config, event_sequence),
        cells, grid);
}

DivisionProposal make_division_proposal(Slot mother,
                                        const CellStore3D& cells,
                                        const SparseChunkGrid3D& grid,
                                        const Model3DConfig& config) {
    DivisionProposal proposal;
    if (!cells.valid(mother)) return proposal;
    proposal.mother = mother;
    proposal.mother_uid = cells.uid(mother);
    proposal.event_sequence = cells.event_sequence(mother);
    proposal.mother_anchor = cells.anchor(mother);
    proposal.mother_stage = cells.stage(mother);

    if (proposal.mother_stage == CellStage::large) {
        std::vector<Vec3i> candidates;
        for (const Vec3i anchor : chebyshev_shell(
                 proposal.mother_anchor, config.division_shell_radius)) {
            if (grid.can_place_large(anchor)) candidates.push_back(anchor);
        }
        if (!candidates.empty()) {
            const auto index = rng_bounded(
                config.seed, proposal.mother_uid,
                static_cast<std::uint64_t>(RngEventKind::division_location),
                proposal.event_sequence, candidates.size());
            proposal.action = DivisionAction::large_daughter;
            proposal.daughter_target = candidates[static_cast<std::size_t>(index)];
            const auto footprint = large_footprint(proposal.daughter_target);
            proposal.reserved_sites.assign(footprint.begin(), footprint.end());
            return proposal;
        }
        if (!config.allow_shape_reduction) return proposal;
        std::vector<Vec3i> sites;
        for (const Vec3i site : shape_reduction_sites(proposal.mother_anchor)) {
            if (config.thin_layer && site.z != proposal.mother_anchor.z) continue;
            const Slot owner = grid.owner(site);
            if (!grid.blocked_by_vessel(site) &&
                (owner == mother || grid.available(site))) {
                sites.push_back(site);
            }
        }
        if (sites.size() < 2) return proposal;
        deterministic_shuffle(
            sites.begin(), sites.end(), config.seed, proposal.mother_uid,
            static_cast<std::uint64_t>(RngEventKind::division_location),
            proposal.event_sequence);
        proposal.action = DivisionAction::shape_reduction;
        proposal.mother_target = sites[0];
        proposal.daughter_target = sites[1];
        proposal.reserved_sites = {proposal.mother_target,
                                   proposal.daughter_target};
        return proposal;
    }

    if (proposal.mother_stage == CellStage::ultrasmall) {
        proposal.stage_recovery = make_stage_recovery_proposal(
            mother, cells, grid, config, proposal.event_sequence);
        if (proposal.stage_recovery.action != StageRecoveryAction::none) {
            proposal.action = DivisionAction::stage_recovery;
            proposal.reserved_sites = proposal.stage_recovery.reserved_sites;
            proposal.locks_colocation_group =
                proposal.stage_recovery.locks_colocation_group;
        }
        return proposal;
    }

    const std::vector<DirectionId> directions = free_small_directions(
        proposal.mother_anchor, grid, config.thin_layer);
    if (!directions.empty()) {
        const auto index = rng_bounded(
            config.seed, proposal.mother_uid,
            static_cast<std::uint64_t>(RngEventKind::division_location),
            proposal.event_sequence, directions.size());
        proposal.action = DivisionAction::small_daughter;
        proposal.daughter_target = proposal.mother_anchor + direction_vector(
            directions[static_cast<std::size_t>(index)]);
        proposal.reserved_sites = {proposal.daughter_target};
        return proposal;
    }
    if (cells.type(mother) == CellType::r) {
        proposal.action = DivisionAction::remove_r_mother;
    } else if (config.ultrasmall_enabled) {
        proposal.action = DivisionAction::colocated_daughter;
        proposal.daughter_target = proposal.mother_anchor;
        proposal.locks_colocation_group = true;
    }
    return proposal;
}

DivisionResult commit_division_proposal(
    const DivisionProposal& proposal,
    double now,
    CellUid& next_uid,
    CellStore3D& cells,
    SparseChunkGrid3D& grid,
    BlockDensityIndex3D& density,
    const Model3DConfig& config,
    std::vector<LineageEdge>& lineage,
    const VascularInfluenceField3D* vascular_influence) {
    DivisionResult result;
    const Slot mother = proposal.mother;
    if (proposal.action == DivisionAction::none || !cells.valid(mother) ||
        cells.uid(mother) != proposal.mother_uid ||
        cells.anchor(mother) != proposal.mother_anchor ||
        cells.stage(mother) != proposal.mother_stage ||
        cells.event_sequence(mother) != proposal.event_sequence) {
        return result;
    }
    const Vec3i mother_anchor = cells.anchor(mother);
    const CellUid mother_uid = cells.uid(mother);
    const CellType mother_type = cells.type(mother);
    const CellStage mother_stage = cells.stage(mother);
    // Division attempts are transactional. Read the sequence for stateless
    // proposal RNG, but advance it only when an actual biological/geometric
    // result commits. A failed retry must not change future evolution.
    const std::uint64_t sequence = cells.event_sequence(mother);
    const bool conversion_eligible = config.r_to_K_conversion.enabled &&
        mother_type == CellType::r;
    const double conversion_density = conversion_eligible
        ? r_to_K_division_density(
              density, mother_anchor, mother_stage, config.r_to_K_conversion,
              config.thin_layer)
        : 0.0;
    const double variation = config.division_timing.inherited_growth_multiplier_min +
        (config.division_timing.inherited_growth_multiplier_max -
         config.division_timing.inherited_growth_multiplier_min) *
        rng_unit(config.seed, mother_uid,
                 static_cast<std::uint64_t>(RngEventKind::division_timing), sequence, 1);
    const double mother_growth_cap = mother_type == CellType::r
        ? config.division_timing.r_max_inherent_growth_rate
        : config.division_timing.K_max_inherent_growth_rate;
    const float inherited_growth = static_cast<float>(std::min(
        static_cast<double>(cells.inherent_growth_rate(mother)) * variation,
        mother_growth_cap));

    auto finish_daughter = [&](Slot daughter) {
        result.changed = true;
        result.daughter = daughter;
        cells.set_event_sequence(mother, sequence + 1);
        cells.set_inherent_growth_rate(mother, inherited_growth);
        const bool convert = conversion_eligible &&
            should_convert_r_daughter(
                mother_type, conversion_density,
                config.r_to_K_conversion, config.seed, mother_uid, sequence);
        if (convert) {
            const InitialCellRates3D rates = sample_initial_cell_rates(
                config, CellType::K, config.seed, cells.uid(daughter));
            const float growth = static_cast<float>(std::min(
                rates.inherent_growth_rate,
                config.division_timing.K_max_inherent_growth_rate));
            cells.set_type(daughter, CellType::K);
            cells.set_inherent_growth_rate(daughter, growth);
            cells.set_density_growth_rate(daughter, growth);
            (void)rates.migration_rate;
        }
        // A successful biological division starts fresh migration state for
        // both cells. Rates are keyed by immutable UID and the committed
        // division sequence; an active mother is never blindly copied to its
        // daughter, and an r->K daughter receives K rates.
        reset_migration_for_new_cycle(mother, cells, config, sequence);
        reset_migration_for_new_cycle(daughter, cells, config, sequence);
        lineage.push_back({now, cells.uid(daughter), mother_uid, cells.clone_id(daughter), cells.type(daughter)});
        density.add(cells.anchor(daughter), cells.type(daughter), daughter);
        refresh_growth_state(mother, now, cells, density, config, vascular_influence);
        refresh_growth_state(daughter, now, cells, density, config, vascular_influence);
        initialize_division_cycle(mother, now, cells, config);
        initialize_division_cycle(daughter, now, cells, config);
        (void)refresh_migration_activation_state(
            mother, now, cells, density, config);
        (void)refresh_migration_activation_state(
            daughter, now, cells, density, config);
    };

    if (mother_stage == CellStage::large) {
        if (proposal.action == DivisionAction::large_daughter) {
            if (!grid.can_place_large(proposal.daughter_target)) return result;
            CellInit daughter = daughter_init(
                cells, mother, next_uid, proposal.daughter_target,
                CellStage::large, now, inherited_growth);
            const Slot daughter_slot = cells.create(daughter);
            if (!grid.place_large(daughter.anchor, daughter_slot)) {
                cells.erase(daughter_slot);
                return result;
            }
            ++next_uid;
            result.changed_sites = {mother_anchor, daughter.anchor};
            finish_daughter(daughter_slot);
            return result;
        }
        if (proposal.action != DivisionAction::shape_reduction ||
            proposal.mother_target == proposal.daughter_target) {
            return result;
        }
        for (const Vec3i site : {proposal.mother_target,
                                 proposal.daughter_target}) {
            const Slot owner = grid.owner(site);
            if (grid.blocked_by_vessel(site) ||
                (owner != mother && !grid.available(site))) {
                return result;
            }
        }
        CellInit daughter = daughter_init(cells, mother, next_uid,
                                          proposal.daughter_target,
                                          CellStage::small, now, inherited_growth);
        // Allocate the daughter slot before changing occupancy so allocation
        // failure cannot leave the mother half-converted.
        const Slot daughter_slot = cells.create(daughter);
        grid.remove_large(mother_anchor, mother);
        if (!grid.place_single(proposal.mother_target, mother)) {
            cells.erase(daughter_slot);
            if (!grid.place_large(mother_anchor, mother)) {
                throw std::logic_error("shape-reduction rollback could not restore mother footprint");
            }
            return result;
        }
        if (!grid.place_single(proposal.daughter_target, daughter_slot)) {
            cells.erase(daughter_slot);
            if (!grid.remove(proposal.mother_target, mother) ||
                !grid.place_large(mother_anchor, mother)) {
                throw std::logic_error("shape-reduction rollback could not restore mother footprint");
            }
            return result;
        }
        ++next_uid;
        cells.set_anchor(mother, proposal.mother_target);
        cells.set_stage(mother, CellStage::small);
        density.move(mother_anchor, proposal.mother_target,
                     cells.type(mother), mother);
        result.changed_sites = {mother_anchor, proposal.mother_target,
                                proposal.daughter_target};
        finish_daughter(daughter_slot);
        return result;
    }

    if (mother_stage == CellStage::ultrasmall) {
        if (proposal.action == DivisionAction::stage_recovery &&
            commit_stage_recovery_proposal(
                proposal.stage_recovery, cells, grid)) {
            cells.set_event_sequence(mother, sequence + 1);
            density.move(mother_anchor, cells.anchor(mother), cells.type(mother), mother);
            result.changed = true;
            result.stage_recovery = true;
            result.changed_sites = {mother_anchor, cells.anchor(mother)};
        }
        // A co-located cell that cannot separate remains in the explicit
        // stage-2 group. Lack of space is a retryable crowding condition, not
        // the r-cell failed-division viability rule used below for stage 1.
        return result;
    }

    if (proposal.action == DivisionAction::small_daughter) {
        if (!grid.available(proposal.daughter_target)) return result;
        CellInit daughter = daughter_init(cells, mother, next_uid,
                                          proposal.daughter_target,
                                          CellStage::small, now, inherited_growth);
        const Slot daughter_slot = cells.create(daughter);
        if (!grid.place_single(proposal.daughter_target, daughter_slot)) {
            cells.erase(daughter_slot);
            return result;
        }
        ++next_uid;
        result.changed_sites = {mother_anchor, proposal.daughter_target};
        finish_daughter(daughter_slot);
        return result;
    }

    if (proposal.action == DivisionAction::remove_r_mother &&
        cells.type(mother) == CellType::r) {
        result.mother_removed = remove_cell(mother, cells, grid, density);
        result.changed = result.mother_removed;
        result.changed_sites = {mother_anchor};
        return result;
    }
    if (proposal.action != DivisionAction::colocated_daughter ||
        cells.type(mother) != CellType::K || !config.ultrasmall_enabled) {
        return result;
    }
    CellInit daughter = daughter_init(cells, mother, next_uid, mother_anchor,
                                      CellStage::ultrasmall, now, inherited_growth);
    const Slot daughter_slot = cells.create(daughter);
    if (!grid.add_colocated(mother_anchor, daughter_slot)) {
        cells.erase(daughter_slot);
        return result;
    }
    ++next_uid;
    cells.set_stage(mother, CellStage::ultrasmall);
    result.changed_sites = {mother_anchor};
    finish_daughter(daughter_slot);
    return result;
}

DivisionResult divide_cell(Slot mother,
                           double now,
                           CellUid& next_uid,
                           CellStore3D& cells,
                           SparseChunkGrid3D& grid,
                           BlockDensityIndex3D& density,
                           const Model3DConfig& config,
                           std::vector<LineageEdge>& lineage,
                           const VascularInfluenceField3D* vascular_influence) {
    return commit_division_proposal(
        make_division_proposal(mother, cells, grid, config), now, next_uid,
        cells, grid, density, config, lineage, vascular_influence);
}

}  // namespace atcg3d
