#include "rules/migration.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_set>

#include "core/stateless_rng.hpp"
#include "geometry/directions.hpp"
#include "geometry/footprint.hpp"

namespace atcg3d {
namespace {

DirectionId weighted_choice(const std::vector<DirectionId>& candidates,
                            double distance_weight_exponent,
                            std::uint64_t seed,
                            CellUid uid,
                            std::uint64_t event_sequence,
                            std::uint64_t draw) {
    if (candidates.empty()) {
        return kStayDirection;
    }
    if (distance_weight_exponent == 0.0) {
        const auto index = rng_bounded(seed, uid,
                                       static_cast<std::uint64_t>(RngEventKind::migration_direction),
                                       event_sequence, candidates.size(), draw);
        return candidates[static_cast<std::size_t>(index)];
    }
    double total = 0.0;
    std::vector<double> cumulative;
    cumulative.reserve(candidates.size());
    for (const DirectionId direction : candidates) {
        const double length = std::sqrt(static_cast<double>(squared_length(direction_vector(direction))));
        total += std::pow(length, -distance_weight_exponent);
        cumulative.push_back(total);
    }
    const double target = rng_unit(seed, uid,
                                   static_cast<std::uint64_t>(RngEventKind::migration_direction),
                                   event_sequence, draw) * total;
    const auto iterator = std::lower_bound(cumulative.begin(), cumulative.end(), target);
    const auto index = static_cast<std::size_t>(std::distance(cumulative.begin(), iterator));
    return candidates[std::min(index, candidates.size() - 1)];
}

std::vector<DirectionId> density_filtered(const std::vector<DirectionId>& directions,
                                          Vec3i anchor,
                                          const BlockDensityIndex3D& density,
                                          const Model3DConfig& config) {
    std::vector<DirectionId> result;
    result.reserve(directions.size());
    for (const DirectionId direction : directions) {
        if (density.estimate_directional_density(anchor, direction,
                                                 config.direction_density_radius,
                                                 config.direction_density_half_angle_degrees) <=
            config.direction_density_threshold) {
            result.push_back(direction);
        }
    }
    return result;
}

}  // namespace

std::vector<DirectionId> feasible_directions(Slot slot,
                                             const CellStore3D& cells,
                                             const SparseChunkGrid3D& grid,
                                             bool thin_layer) {
    std::vector<DirectionId> result;
    if (!cells.valid(slot)) {
        return result;
    }
    result.reserve(thin_layer ? 8 : 26);
    const Vec3i anchor = cells.anchor(slot);
    for (DirectionId direction = 1; direction <= 26; ++direction) {
        const Vec3i displacement = direction_vector(direction);
        if (thin_layer && displacement.z != 0) {
            continue;
        }
        const bool allowed = cells.stage(slot) == CellStage::large
            ? grid.can_move_large(anchor, displacement)
            : grid.available(anchor + displacement);
        if (allowed) {
            result.push_back(direction);
        }
    }
    return result;
}

DirectionId select_migration_direction(Slot slot,
                                       const CellStore3D& cells,
                                       const SparseChunkGrid3D& grid,
                                       const BlockDensityIndex3D& density,
                                       const Model3DConfig& config,
                                       std::uint64_t event_sequence) {
    std::vector<DirectionId> feasible = feasible_directions(slot, cells, grid, config.thin_layer);
    if (feasible.empty()) {
        return kStayDirection;
    }
    const CellUid uid = cells.uid(slot);
    if (cells.type(slot) == CellType::K) {
        return weighted_choice(feasible, config.distance_weight_exponent, config.seed, uid, event_sequence, 0);
    }

    const DirectionId previous = cells.last_direction(slot);
    if (previous == kStayDirection) {
        feasible = density_filtered(feasible, cells.anchor(slot), density, config);
        return weighted_choice(feasible, config.distance_weight_exponent, config.seed, uid, event_sequence, 0);
    }

    if (config.persistence_uses_density) {
        feasible = density_filtered(feasible, cells.anchor(slot), density, config);
    }
    const bool forward_available = std::find(feasible.begin(), feasible.end(), previous) != feasible.end();
    std::vector<DirectionId> turns;
    for (const DirectionId direction : feasible) {
        if (direction != previous &&
            direction_angle_degrees(previous, direction) <= config.turn_half_angle_degrees + 1e-10) {
            turns.push_back(direction);
        }
    }
    if (forward_available && turns.empty()) {
        return previous;
    }
    if (forward_available) {
        if (rng_unit(config.seed, uid,
                     static_cast<std::uint64_t>(RngEventKind::migration_direction),
                     event_sequence, 0) < config.continue_probability) {
            return previous;
        }
        return weighted_choice(turns, config.distance_weight_exponent, config.seed, uid, event_sequence, 1);
    }
    return weighted_choice(turns, config.distance_weight_exponent, config.seed, uid, event_sequence, 0);
}

MoveProposal make_move_proposal(Slot slot,
                                const CellStore3D& cells,
                                const SparseChunkGrid3D& grid,
                                const BlockDensityIndex3D& density,
                                const Model3DConfig& config,
                                std::uint64_t event_sequence,
                                std::uint64_t time_bucket) {
    MoveProposal proposal;
    if (!cells.valid(slot)) {
        return proposal;
    }
    proposal.slot = slot;
    proposal.uid = cells.uid(slot);
    proposal.direction = select_migration_direction(slot, cells, grid, density, config, event_sequence);
    proposal.from = cells.anchor(slot);
    proposal.to = proposal.from + direction_vector(proposal.direction);
    proposal.priority = rng_word(config.seed, proposal.uid,
                                 static_cast<std::uint64_t>(RngEventKind::conflict_priority),
                                 event_sequence, time_bucket);
    if (proposal.direction == kStayDirection) {
        proposal.to = proposal.from;
        return proposal;
    }
    if (cells.stage(slot) == CellStage::large) {
        proposal.reserved_sites = entering_voxels(proposal.from, direction_vector(proposal.direction));
    } else {
        proposal.reserved_sites.push_back(proposal.to);
    }
    return proposal;
}

bool commit_move(const MoveProposal& proposal,
                 CellStore3D& cells,
                 SparseChunkGrid3D& grid,
                 BlockDensityIndex3D& density) {
    if (proposal.direction == kStayDirection || !cells.valid(proposal.slot) ||
        cells.uid(proposal.slot) != proposal.uid || cells.anchor(proposal.slot) != proposal.from) {
        if (cells.valid(proposal.slot)) {
            cells.set_last_direction(proposal.slot, kStayDirection);
        }
        return false;
    }
    if (!std::all_of(proposal.reserved_sites.begin(), proposal.reserved_sites.end(),
                     [&grid](Vec3i site) { return grid.available(site); })) {
        return false;
    }
    const CellType type = cells.type(proposal.slot);
    const CellStage stage = cells.stage(proposal.slot);
    if (stage == CellStage::large) {
        grid.remove_large(proposal.from, proposal.slot);
        if (!grid.place_large(proposal.to, proposal.slot)) {
            grid.place_large(proposal.from, proposal.slot);
            return false;
        }
    } else {
        const std::vector<Slot> former_group = grid.occupants(proposal.from);
        if (!grid.remove(proposal.from, proposal.slot) || !grid.place_single(proposal.to, proposal.slot)) {
            grid.add_colocated(proposal.from, proposal.slot);
            return false;
        }
        cells.set_stage(proposal.slot, CellStage::small);
        const std::vector<Slot> remaining = grid.occupants(proposal.from);
        for (const Slot other : remaining) {
            if (cells.valid(other)) {
                cells.set_stage(other, remaining.size() == 1 ? CellStage::small : CellStage::ultrasmall);
            }
        }
        (void)former_group;
    }
    density.move(proposal.from, proposal.to, type);
    cells.set_anchor(proposal.slot, proposal.to);
    cells.set_last_direction(proposal.slot, proposal.direction);
    return true;
}

}  // namespace atcg3d
