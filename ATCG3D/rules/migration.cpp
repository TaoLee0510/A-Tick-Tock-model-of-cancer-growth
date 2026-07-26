#include "rules/migration.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <vector>

#include "core/stateless_rng.hpp"
#include "geometry/directions.hpp"
#include "geometry/footprint.hpp"

namespace atcg3d {
namespace {

DirectionId weighted_choice(const DirectionCandidates3D& candidates,
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
    for (const DirectionId direction : candidates) {
        const double length = std::sqrt(static_cast<double>(squared_length(direction_vector(direction))));
        total += std::pow(length, -distance_weight_exponent);
    }
    const double target = rng_unit(seed, uid,
                                   static_cast<std::uint64_t>(RngEventKind::migration_direction),
                                   event_sequence, draw) * total;
    double cumulative = 0.0;
    for (const DirectionId direction : candidates) {
        const double length = std::sqrt(
            static_cast<double>(squared_length(direction_vector(direction))));
        cumulative += std::pow(length, -distance_weight_exponent);
        if (target <= cumulative) return direction;
    }
    return candidates[candidates.size() - 1];
}

DirectionCandidates3D density_filtered(
    const DirectionCandidates3D& directions,
    Vec3i anchor,
    const BlockDensityIndex3D& density,
    const Model3DConfig& config) {
    DirectionCandidates3D result;
    const auto directional_densities =
        density.estimate_all_directional_densities(
            anchor, config.direction_density_radius,
            config.direction_density_half_angle_degrees, config.thin_layer);
    for (const DirectionId direction : directions) {
        if (directional_densities[direction] <=
            config.direction_density_threshold) {
            result.push_back(direction);
        }
    }
    return result;
}

DirectionId select_from_candidates(
    Slot slot,
    DirectionCandidates3D candidates,
    const CellStore3D& cells,
    const Model3DConfig& config,
    std::uint64_t event_sequence) {
    if (candidates.empty()) return kStayDirection;
    const CellUid uid = cells.uid(slot);
    const bool activated = config.migration_activation_enabled &&
        (cells.flags(slot) & static_cast<std::uint8_t>(kMigrationActive)) != 0;
    const DirectionId previous = cells.last_direction(slot);
    if (cells.type(slot) == CellType::K || !activated ||
        previous == kStayDirection) {
        return weighted_choice(candidates, config.distance_weight_exponent,
                               config.seed, uid, event_sequence, 0);
    }

    const bool forward_available =
        std::find(candidates.begin(), candidates.end(), previous) !=
        candidates.end();
    DirectionCandidates3D turns;
    for (const DirectionId direction : candidates) {
        if (direction != previous &&
            direction_angle_degrees(previous, direction) <=
                config.turn_half_angle_degrees + 1e-10) {
            turns.push_back(direction);
        }
    }
    if (forward_available && turns.empty()) return previous;
    if (forward_available &&
        rng_unit(config.seed, uid,
                 static_cast<std::uint64_t>(
                     RngEventKind::migration_direction),
                 event_sequence, 0) < config.continue_probability) {
        return previous;
    }
    return weighted_choice(turns, config.distance_weight_exponent,
                           config.seed, uid, event_sequence,
                           forward_available ? 1 : 0);
}

}  // namespace

DirectionCandidates3D feasible_directions(
    Slot slot,
    const CellStore3D& cells,
    const SparseChunkGrid3D& grid,
    bool thin_layer) {
    DirectionCandidates3D result;
    if (!cells.valid(slot)) {
        return result;
    }
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
    DirectionCandidates3D feasible =
        feasible_directions(slot, cells, grid, config.thin_layer);
    if (feasible.empty()) {
        return kStayDirection;
    }
    const CellUid uid = cells.uid(slot);
    const bool activated = config.migration_activation_enabled &&
        (cells.flags(slot) & static_cast<std::uint8_t>(kMigrationActive)) != 0;
    if (cells.type(slot) == CellType::K || !activated) {
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
    DirectionCandidates3D turns;
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

DirectionCandidates3D feasible_crowding_swap_directions(
    Slot slot,
    const CellStore3D& cells,
    const SparseChunkGrid3D& grid,
    bool thin_layer) {
    DirectionCandidates3D result;
    if (!cells.valid(slot) || cells.stage(slot) != CellStage::small) {
        return result;
    }
    const Vec3i anchor = cells.anchor(slot);
    for (DirectionId direction = 1; direction <= 26; ++direction) {
        const Vec3i displacement = direction_vector(direction);
        if (thin_layer && displacement.z != 0) continue;
        const Vec3i target = anchor + displacement;
        if (grid.blocked_by_vessel(target)) continue;
        const Slot partner = grid.owner(target);
        if (partner == kEmptySlot || partner == slot ||
            !cells.valid(partner) ||
            cells.stage(partner) != CellStage::small ||
            cells.anchor(partner) != target) {
            continue;
        }
        const std::vector<Slot> occupants = grid.occupants(target);
        if (occupants.size() == 1 && occupants.front() == partner) {
            result.push_back(direction);
        }
    }
    return result;
}

DirectionId select_crowding_swap_direction(
    Slot slot,
    const CellStore3D& cells,
    const SparseChunkGrid3D& grid,
    const Model3DConfig& config,
    std::uint64_t event_sequence) {
    return select_from_candidates(
        slot,
        feasible_crowding_swap_directions(
            slot, cells, grid, config.thin_layer),
        cells, config, event_sequence);
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
                                 time_bucket, 0);
    if (proposal.direction == kStayDirection) {
        proposal.to = proposal.from;
        return proposal;
    }
    if (cells.stage(slot) == CellStage::large) {
        const std::vector<Vec3i> entering =
            entering_voxels(
                proposal.from, direction_vector(proposal.direction));
        proposal.reserved_sites.assign(entering.begin(), entering.end());
    } else {
        proposal.reserved_sites.push_back(proposal.to);
    }
    return proposal;
}

MoveProposal make_crowding_swap_proposal(
    Slot slot,
    DirectionId direction,
    const CellStore3D& cells,
    const SparseChunkGrid3D& grid,
    const Model3DConfig& config,
    std::uint64_t time_bucket) {
    MoveProposal proposal;
    if (!cells.valid(slot) || cells.stage(slot) != CellStage::small ||
        direction == kStayDirection || direction > 26) {
        return proposal;
    }
    proposal.slot = slot;
    proposal.uid = cells.uid(slot);
    proposal.direction = direction;
    proposal.from = cells.anchor(slot);
    proposal.to = proposal.from + direction_vector(direction);
    proposal.priority = rng_word(
        config.seed, proposal.uid,
        static_cast<std::uint64_t>(RngEventKind::conflict_priority),
        time_bucket, 0);
    const std::vector<Slot> occupants = grid.occupants(proposal.to);
    if (occupants.size() != 1 || occupants.front() == slot ||
        !cells.valid(occupants.front()) ||
        cells.stage(occupants.front()) != CellStage::small ||
        cells.anchor(occupants.front()) != proposal.to ||
        grid.blocked_by_vessel(proposal.to)) {
        proposal.direction = kStayDirection;
        proposal.to = proposal.from;
        return proposal;
    }
    proposal.swaps_anchors = true;
    proposal.swap_partner = occupants.front();
    proposal.swap_partner_uid = cells.uid(proposal.swap_partner);
    proposal.reserved_sites = {proposal.from, proposal.to};
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
    density.move(proposal.from, proposal.to, type, proposal.slot);
    cells.set_anchor(proposal.slot, proposal.to);
    cells.set_last_direction(proposal.slot, proposal.direction);
    return true;
}

bool commit_crowding_swap(const MoveProposal& proposal,
                          CellStore3D& cells,
                          SparseChunkGrid3D& grid,
                          BlockDensityIndex3D& density) {
    if (!proposal.swaps_anchors ||
        proposal.direction == kStayDirection ||
        !cells.valid(proposal.slot) ||
        !cells.valid(proposal.swap_partner) ||
        cells.uid(proposal.slot) != proposal.uid ||
        cells.uid(proposal.swap_partner) != proposal.swap_partner_uid ||
        cells.stage(proposal.slot) != CellStage::small ||
        cells.stage(proposal.swap_partner) != CellStage::small ||
        cells.anchor(proposal.slot) != proposal.from ||
        cells.anchor(proposal.swap_partner) != proposal.to ||
        grid.blocked_by_vessel(proposal.from) ||
        grid.blocked_by_vessel(proposal.to)) {
        return false;
    }
    const std::vector<Slot> from_occupants = grid.occupants(proposal.from);
    const std::vector<Slot> to_occupants = grid.occupants(proposal.to);
    if (from_occupants.size() != 1 ||
        from_occupants.front() != proposal.slot ||
        to_occupants.size() != 1 ||
        to_occupants.front() != proposal.swap_partner) {
        return false;
    }

    if (!grid.remove(proposal.from, proposal.slot) ||
        !grid.remove(proposal.to, proposal.swap_partner)) {
        if (grid.owner(proposal.from) == kEmptySlot) {
            (void)grid.place_single(proposal.from, proposal.slot);
        }
        return false;
    }
    const bool actor_placed =
        grid.place_single(proposal.to, proposal.slot);
    const bool partner_placed =
        actor_placed &&
        grid.place_single(proposal.from, proposal.swap_partner);
    if (!partner_placed) {
        if (actor_placed) {
            (void)grid.remove(proposal.to, proposal.slot);
        }
        (void)grid.place_single(proposal.from, proposal.slot);
        (void)grid.place_single(proposal.to, proposal.swap_partner);
        return false;
    }

    const CellType actor_type = cells.type(proposal.slot);
    const CellType partner_type = cells.type(proposal.swap_partner);
    density.move(proposal.from, proposal.to, actor_type, proposal.slot);
    density.move(proposal.to, proposal.from, partner_type,
                 proposal.swap_partner);
    cells.set_anchor(proposal.slot, proposal.to);
    cells.set_anchor(proposal.swap_partner, proposal.from);
    cells.set_last_direction(proposal.slot, proposal.direction);
    cells.set_last_direction(
        proposal.swap_partner, opposite_direction(proposal.direction));
    cells.clear_swap_wait(proposal.slot);
    cells.clear_swap_wait(proposal.swap_partner);
    return true;
}

}  // namespace atcg3d
