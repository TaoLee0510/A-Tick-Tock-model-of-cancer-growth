#include "rules/lifecycle.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_set>

#include "core/stateless_rng.hpp"
#include "geometry/directions.hpp"
#include "geometry/footprint.hpp"
#include "rules/density.hpp"
#include "rules/migration.hpp"

namespace atcg3d {
namespace {

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
    daughter.last_update_time = now;
    daughter.event_sequence = 0;
    return daughter;
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

double sample_division_delay(double density_growth_rate,
                             std::uint64_t seed,
                             CellUid uid,
                             std::uint64_t event_sequence) {
    if (!(density_growth_rate > 0.0)) {
        return 0.0;
    }
    const double expected = 24.0 / density_growth_rate;
    const double minimum = 0.9 * expected;
    const double range = 0.1 * expected;
    const double probability = std::clamp(1.0 / std::max(range, 1.0), 1e-12, 1.0);
    if (probability >= 1.0) {
        return minimum;
    }
    const double uniform = std::clamp(
        rng_unit(seed, uid, static_cast<std::uint64_t>(RngEventKind::division_timing), event_sequence),
        1e-15, 1.0 - 1e-15);
    const double geometric = std::floor(std::log1p(-uniform) / std::log1p(-probability));
    return minimum + geometric;
}

void refresh_growth_state(Slot slot,
                          double now,
                          CellStore3D& cells,
                          const BlockDensityIndex3D& density,
                          const Model3DConfig& config) {
    if (!cells.valid(slot)) {
        return;
    }
    const double rate = density_growth_rate_for_cell(cells, slot, density, config);
    cells.set_density_growth_rate(slot, static_cast<float>(rate));
    cells.set_last_update_time(slot, now);
    if (config.migration_activation_enabled) {
        const double migration_density = migration_activation_density(
            density, cells.anchor(slot), cells.stage(slot),
            config.migration_activation_window_edge,
            config.migration_activation_block_edge);
        std::uint8_t flags = cells.flags(slot) & static_cast<std::uint8_t>(~kMigrationActive);
        if (migration_density >= config.migration_activation_threshold) {
            flags |= kMigrationActive;
        }
        cells.set_flags(slot, flags);
    }
    const std::uint64_t sequence = cells.consume_event_sequence(slot);
    if (rate > 0.0) {
        cells.set_death_deadline(slot, 0.0);
        cells.set_next_division_time(
            slot, now + sample_division_delay(rate, config.seed, cells.uid(slot), sequence));
    } else {
        cells.set_next_division_time(slot, 0.0);
        if (cells.death_deadline(slot) <= now) {
            const double delay = cells.type(slot) == CellType::r
                ? config.r_death_delay_hours : config.K_death_delay_hours;
            cells.set_death_deadline(slot, now + delay);
        }
    }
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
    density.remove(anchor, type);
    cells.erase(slot);
    if (stage != CellStage::large) {
        update_colocation_stages(anchor, cells, grid);
    }
    return true;
}

bool try_stage_recovery(Slot slot,
                        CellStore3D& cells,
                        SparseChunkGrid3D& grid,
                        const Model3DConfig& config,
                        std::uint64_t event_sequence) {
    if (!cells.valid(slot) || cells.stage(slot) == CellStage::large) {
        return false;
    }
    const Vec3i site = cells.anchor(slot);
    const std::vector<Slot> group = grid.occupants(site);
    if (group.size() > 1) {
        const auto free_directions = free_small_directions(site, grid, config.thin_layer);
        if (free_directions.empty()) {
            return false;
        }
        const auto index = rng_bounded(config.seed, cells.uid(slot),
                                       static_cast<std::uint64_t>(RngEventKind::division_location),
                                       event_sequence, free_directions.size());
        const Vec3i target = site + direction_vector(free_directions[static_cast<std::size_t>(index)]);
        if (!grid.remove(site, slot) || !grid.place_single(target, slot)) {
            grid.add_colocated(site, slot);
            return false;
        }
        cells.set_anchor(slot, target);
        cells.set_stage(slot, CellStage::small);
        update_colocation_stages(site, cells, grid);
        return true;
    }

    std::vector<Vec3i> candidates;
    for (const Vec3i anchor : stage_recovery_anchors(site)) {
        bool available = true;
        for (const Vec3i voxel : large_footprint(anchor)) {
            const Slot owner = grid.owner(voxel);
            if (owner != kEmptySlot && owner != slot) {
                available = false;
                break;
            }
        }
        if (available) {
            candidates.push_back(anchor);
        }
    }
    if (candidates.empty()) {
        return false;
    }
    const auto index = rng_bounded(config.seed, cells.uid(slot),
                                   static_cast<std::uint64_t>(RngEventKind::division_location),
                                   event_sequence, candidates.size());
    const Vec3i anchor = candidates[static_cast<std::size_t>(index)];
    grid.remove(site, slot);
    if (!grid.place_large(anchor, slot)) {
        grid.place_single(site, slot);
        return false;
    }
    cells.set_anchor(slot, anchor);
    cells.set_stage(slot, CellStage::large);
    return true;
}

DivisionResult divide_cell(Slot mother,
                           double now,
                           CellUid& next_uid,
                           CellStore3D& cells,
                           SparseChunkGrid3D& grid,
                           BlockDensityIndex3D& density,
                           const Model3DConfig& config,
                           std::vector<LineageEdge>& lineage) {
    DivisionResult result;
    if (!cells.valid(mother)) {
        return result;
    }
    const Vec3i mother_anchor = cells.anchor(mother);
    const CellUid mother_uid = cells.uid(mother);
    const std::uint64_t sequence = cells.consume_event_sequence(mother);
    const double variation = 0.95 + 0.10 * rng_unit(
        config.seed, mother_uid, static_cast<std::uint64_t>(RngEventKind::division_timing), sequence, 1);
    const float inherited_growth = static_cast<float>(cells.inherent_growth_rate(mother) * variation);
    cells.set_inherent_growth_rate(mother, inherited_growth);

    auto finish_daughter = [&](Slot daughter) {
        result.changed = true;
        result.daughter = daughter;
        lineage.push_back({now, cells.uid(daughter), mother_uid, cells.clone_id(daughter), cells.type(daughter)});
        density.add(cells.anchor(daughter), cells.type(daughter));
        refresh_growth_state(mother, now, cells, density, config);
        refresh_growth_state(daughter, now, cells, density, config);
    };

    if (cells.stage(mother) == CellStage::large) {
        std::vector<Vec3i> candidates;
        for (const Vec3i anchor : chebyshev_shell(mother_anchor, config.division_shell_radius)) {
            if (grid.can_place_large(anchor)) {
                candidates.push_back(anchor);
            }
        }
        if (!candidates.empty()) {
            const auto index = rng_bounded(config.seed, mother_uid,
                                           static_cast<std::uint64_t>(RngEventKind::division_location),
                                           sequence, candidates.size());
            CellInit daughter = daughter_init(cells, mother, next_uid++,
                                              candidates[static_cast<std::size_t>(index)],
                                              CellStage::large, now, inherited_growth);
            const Slot daughter_slot = cells.create(daughter);
            if (!grid.place_large(daughter.anchor, daughter_slot)) {
                cells.erase(daughter_slot);
                return result;
            }
            result.changed_sites.push_back(mother_anchor);
            result.changed_sites.push_back(daughter.anchor);
            finish_daughter(daughter_slot);
            return result;
        }

        if (!config.allow_shape_reduction) {
            return result;
        }
        std::vector<Vec3i> sites;
        for (const Vec3i site : shape_reduction_sites(mother_anchor)) {
            const Slot owner = grid.owner(site);
            if (owner == kEmptySlot || owner == mother) {
                sites.push_back(site);
            }
        }
        if (sites.size() < 2) {
            return result;
        }
        deterministic_shuffle(sites.begin(), sites.end(), config.seed, mother_uid,
                              static_cast<std::uint64_t>(RngEventKind::division_location), sequence);
        const Vec3i mother_site = sites[0];
        const Vec3i daughter_site = sites[1];
        grid.remove_large(mother_anchor, mother);
        cells.set_anchor(mother, mother_site);
        cells.set_stage(mother, CellStage::small);
        grid.place_single(mother_site, mother);
        density.move(mother_anchor, mother_site, cells.type(mother));
        CellInit daughter = daughter_init(cells, mother, next_uid++, daughter_site,
                                          CellStage::small, now, inherited_growth);
        const Slot daughter_slot = cells.create(daughter);
        if (!grid.place_single(daughter_site, daughter_slot)) {
            cells.erase(daughter_slot);
            return result;
        }
        result.changed_sites = {mother_anchor, mother_site, daughter_site};
        finish_daughter(daughter_slot);
        return result;
    }

    if (cells.stage(mother) == CellStage::ultrasmall) {
        if (try_stage_recovery(mother, cells, grid, config, sequence)) {
            density.move(mother_anchor, cells.anchor(mother), cells.type(mother));
            result.changed = true;
            result.changed_sites = {mother_anchor, cells.anchor(mother)};
            refresh_growth_state(mother, now, cells, density, config);
        } else {
            result.mother_removed = remove_cell(mother, cells, grid, density);
            result.changed = result.mother_removed;
            result.changed_sites = {mother_anchor};
        }
        return result;
    }

    const std::vector<DirectionId> directions = free_small_directions(mother_anchor, grid, config.thin_layer);
    if (!directions.empty()) {
        const auto index = rng_bounded(config.seed, mother_uid,
                                       static_cast<std::uint64_t>(RngEventKind::division_location),
                                       sequence, directions.size());
        const Vec3i daughter_site = mother_anchor + direction_vector(directions[static_cast<std::size_t>(index)]);
        CellInit daughter = daughter_init(cells, mother, next_uid++, daughter_site,
                                          CellStage::small, now, inherited_growth);
        const Slot daughter_slot = cells.create(daughter);
        if (!grid.place_single(daughter_site, daughter_slot)) {
            cells.erase(daughter_slot);
            return result;
        }
        result.changed_sites = {mother_anchor, daughter_site};
        finish_daughter(daughter_slot);
        return result;
    }

    if (cells.type(mother) == CellType::r) {
        result.mother_removed = remove_cell(mother, cells, grid, density);
        result.changed = result.mother_removed;
        result.changed_sites = {mother_anchor};
        return result;
    }
    if (!config.ultrasmall_enabled) {
        return result;
    }
    CellInit daughter = daughter_init(cells, mother, next_uid++, mother_anchor,
                                      CellStage::ultrasmall, now, inherited_growth);
    const Slot daughter_slot = cells.create(daughter);
    if (!grid.add_colocated(mother_anchor, daughter_slot)) {
        cells.erase(daughter_slot);
        return result;
    }
    cells.set_stage(mother, CellStage::ultrasmall);
    result.changed_sites = {mother_anchor};
    finish_daughter(daughter_slot);
    return result;
}

}  // namespace atcg3d
