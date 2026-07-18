#include "vasculature/lesion_index.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <utility>

#include "core/cell_store.hpp"
#include "space/chunk_grid.hpp"
#include "vasculature/surface_index.hpp"

namespace atcg3d {
namespace {

bool coordinate_less(Vec3i lhs, Vec3i rhs) noexcept {
    if (lhs.x != rhs.x) return lhs.x < rhs.x;
    if (lhs.y != rhs.y) return lhs.y < rhs.y;
    return lhs.z < rhs.z;
}

Vec3i checked_offset(Vec3i value, int dx, int dy, int dz, bool& valid) noexcept {
    const std::int64_t x = static_cast<std::int64_t>(value.x) + dx;
    const std::int64_t y = static_cast<std::int64_t>(value.y) + dy;
    const std::int64_t z = static_cast<std::int64_t>(value.z) + dz;
    valid = x >= std::numeric_limits<std::int32_t>::min() &&
            x <= std::numeric_limits<std::int32_t>::max() &&
            y >= std::numeric_limits<std::int32_t>::min() &&
            y <= std::numeric_limits<std::int32_t>::max() &&
            z >= std::numeric_limits<std::int32_t>::min() &&
            z <= std::numeric_limits<std::int32_t>::max();
    return {static_cast<std::int32_t>(x), static_cast<std::int32_t>(y),
            static_cast<std::int32_t>(z)};
}

std::vector<Vec3i> connectivity_offsets(LesionConnectivity3D connectivity) {
    std::vector<Vec3i> offsets;
    offsets.reserve(connectivity == LesionConnectivity3D::face_6 ? 6U : 26U);
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                if (dx == 0 && dy == 0 && dz == 0) continue;
                if (connectivity == LesionConnectivity3D::face_6 &&
                    std::abs(dx) + std::abs(dy) + std::abs(dz) != 1) {
                    continue;
                }
                offsets.push_back({dx, dy, dz});
            }
        }
    }
    std::sort(offsets.begin(), offsets.end(), coordinate_less);
    return offsets;
}

std::int32_t clamp_lattice_coordinate(std::int64_t value) noexcept {
    return static_cast<std::int32_t>(std::clamp(
        value,
        static_cast<std::int64_t>(std::numeric_limits<std::int32_t>::min()),
        static_cast<std::int64_t>(std::numeric_limits<std::int32_t>::max())));
}

Vec3i block_minimum_site(Vec3i coordinate, int edge) noexcept {
    return {clamp_lattice_coordinate(static_cast<std::int64_t>(coordinate.x) * edge),
            clamp_lattice_coordinate(static_cast<std::int64_t>(coordinate.y) * edge),
            clamp_lattice_coordinate(static_cast<std::int64_t>(coordinate.z) * edge)};
}

Vec3i block_maximum_site(Vec3i coordinate, int edge) noexcept {
    return {clamp_lattice_coordinate(static_cast<std::int64_t>(coordinate.x) * edge +
                                     edge - 1),
            clamp_lattice_coordinate(static_cast<std::int64_t>(coordinate.y) * edge +
                                     edge - 1),
            clamp_lattice_coordinate(static_cast<std::int64_t>(coordinate.z) * edge +
                                     edge - 1)};
}

std::int64_t checked_coordinate_sum(std::int64_t current,
                                    std::int64_t delta) {
    if ((delta > 0 && current > std::numeric_limits<std::int64_t>::max() - delta) ||
        (delta < 0 && current < std::numeric_limits<std::int64_t>::min() - delta)) {
        throw std::overflow_error("lesion coordinate sum overflow");
    }
    return current + delta;
}

void validate_sum_within_block(std::int64_t sum,
                               std::uint64_t count,
                               std::int32_t minimum,
                               std::int32_t maximum,
                               const char* field) {
    if (count == 0) {
        if (sum != 0) {
            throw std::invalid_argument(std::string(field) +
                                        " must be zero when count is zero");
        }
        return;
    }
    const long double lower = static_cast<long double>(count) * minimum;
    const long double upper = static_cast<long double>(count) * maximum;
    const long double value = sum;
    if (value < lower || value > upper) {
        throw std::invalid_argument(std::string(field) +
                                    " lies outside its coarse block");
    }
}

}  // namespace

LesionIndex3D::LesionIndex3D(LesionIndexConfig3D config) : config_(config) {
    if (config_.block_edge <= 0) {
        throw std::invalid_argument("lesion block edge must be positive");
    }
    if (config_.connectivity != LesionConnectivity3D::face_6 &&
        config_.connectivity != LesionConnectivity3D::full_26) {
        throw std::invalid_argument("lesion connectivity must be 6 or 26");
    }
    if (!std::isfinite(config_.core_activation_occupied_fraction) ||
        !std::isfinite(config_.core_deactivation_occupied_fraction) ||
        config_.core_activation_occupied_fraction < 0.0 ||
        config_.core_activation_occupied_fraction > 1.0 ||
        config_.core_deactivation_occupied_fraction < 0.0 ||
        config_.core_deactivation_occupied_fraction >
            config_.core_activation_occupied_fraction) {
        throw std::invalid_argument("invalid lesion core occupied-fraction hysteresis");
    }
    if (config_.minimum_cells_per_core_block == 0) {
        throw std::invalid_argument("minimum lesion core cell count must be positive");
    }
    if (!std::isfinite(config_.minimum_biological_volume_per_core_block) ||
        config_.minimum_biological_volume_per_core_block < 0.0) {
        throw std::invalid_argument("minimum lesion biological volume must be finite and nonnegative");
    }
    if (config_.halo_blocks < 0 || config_.halo_blocks > 1024) {
        throw std::invalid_argument("lesion halo block radius must be in [0, 1024]");
    }

    const std::uint64_t edge = static_cast<std::uint64_t>(config_.block_edge);
    if (edge > std::numeric_limits<std::uint64_t>::max() / edge ||
        edge * edge > std::numeric_limits<std::uint64_t>::max() / edge) {
        throw std::invalid_argument("lesion block volume overflows uint64");
    }
    block_voxel_capacity_ = edge * edge * edge;
}

void LesionIndex3D::reset() noexcept {
    blocks_.clear();
    previous_core_identity_.clear();
    core_lesion_by_block_.clear();
    assigned_lesion_by_block_.clear();
    lesions_.clear();
    dirty_blocks_.clear();
    next_lesion_id_ = 1;
    topology_dirty_ = false;
    statistics_dirty_ = false;
}

void LesionIndex3D::begin_rebuild() {
    blocks_.clear();
    core_lesion_by_block_.clear();
    assigned_lesion_by_block_.clear();
    lesions_.clear();
    dirty_blocks_.clear();
    topology_dirty_ = true;
    statistics_dirty_ = true;
}

int LesionIndex3D::floor_div(int value, int divisor) noexcept {
    int quotient = value / divisor;
    const int remainder = value % divisor;
    if (remainder != 0 && ((remainder < 0) != (divisor < 0))) {
        --quotient;
    }
    return quotient;
}

Vec3i LesionIndex3D::block_coordinate(Vec3i site) const noexcept {
    return {floor_div(site.x, config_.block_edge),
            floor_div(site.y, config_.block_edge),
            floor_div(site.z, config_.block_edge)};
}

void LesionIndex3D::add_cell_anchor(Vec3i anchor, double biological_volume) {
    if (!std::isfinite(biological_volume) || biological_volume <= 0.0) {
        throw std::invalid_argument("cell biological volume must be finite and positive");
    }
    const Vec3i coordinate = block_coordinate(anchor);
    const auto existing = blocks_.find(coordinate);
    const bool qualified_before = existing != blocks_.end() &&
                                  qualifies_as_core(coordinate, existing->second);
    BlockAggregate& block = blocks_[coordinate];
    if (block.cell_count == std::numeric_limits<std::uint64_t>::max()) {
        throw std::overflow_error("lesion block cell count overflow");
    }
    const double next_biological_volume = block.biological_volume + biological_volume;
    if (!std::isfinite(next_biological_volume)) {
        throw std::overflow_error("lesion block biological volume overflow");
    }
    const std::int64_t next_x = checked_coordinate_sum(block.cell_sum_x, anchor.x);
    const std::int64_t next_y = checked_coordinate_sum(block.cell_sum_y, anchor.y);
    const std::int64_t next_z = checked_coordinate_sum(block.cell_sum_z, anchor.z);
    ++block.cell_count;
    block.biological_volume = next_biological_volume;
    block.cell_sum_x = next_x;
    block.cell_sum_y = next_y;
    block.cell_sum_z = next_z;
    topology_dirty_ = topology_dirty_ ||
                      qualified_before != qualifies_as_core(coordinate, block);
    statistics_dirty_ = true;
}

void LesionIndex3D::remove_cell_anchor(Vec3i anchor, double biological_volume) {
    if (!std::isfinite(biological_volume) || biological_volume <= 0.0) {
        throw std::invalid_argument("cell biological volume must be finite and positive");
    }
    const Vec3i coordinate = block_coordinate(anchor);
    auto found = blocks_.find(coordinate);
    if (found == blocks_.end() || found->second.cell_count == 0 ||
        found->second.biological_volume + 1e-12 < biological_volume) {
        throw std::logic_error("lesion cell removal underflow");
    }
    const bool qualified_before = qualifies_as_core(coordinate, found->second);
    BlockAggregate& block = found->second;
    const std::int64_t next_x = checked_coordinate_sum(
        block.cell_sum_x, -static_cast<std::int64_t>(anchor.x));
    const std::int64_t next_y = checked_coordinate_sum(
        block.cell_sum_y, -static_cast<std::int64_t>(anchor.y));
    const std::int64_t next_z = checked_coordinate_sum(
        block.cell_sum_z, -static_cast<std::int64_t>(anchor.z));
    --block.cell_count;
    block.biological_volume -= biological_volume;
    if (block.cell_count == 0 || std::abs(block.biological_volume) < 1e-12) {
        block.biological_volume = 0.0;
    }
    block.cell_sum_x = next_x;
    block.cell_sum_y = next_y;
    block.cell_sum_z = next_z;
    erase_block_if_empty(coordinate);
    const auto after = blocks_.find(coordinate);
    const bool qualified_after = after != blocks_.end() &&
                                 qualifies_as_core(coordinate, after->second);
    topology_dirty_ = topology_dirty_ || qualified_before != qualified_after;
    statistics_dirty_ = true;
}

void LesionIndex3D::add_occupied_site(Vec3i site) {
    const Vec3i coordinate = block_coordinate(site);
    const auto existing = blocks_.find(coordinate);
    const bool qualified_before = existing != blocks_.end() &&
                                  qualifies_as_core(coordinate, existing->second);
    BlockAggregate& block = blocks_[coordinate];
    if (block.occupied_voxel_count == std::numeric_limits<std::uint64_t>::max()) {
        throw std::overflow_error("lesion block occupied-voxel count overflow");
    }
    if (block.occupied_voxel_count >= block_voxel_capacity_) {
        throw std::logic_error("occupied site was added more than once to lesion block");
    }
    const std::int64_t next_x = checked_coordinate_sum(block.occupied_sum_x, site.x);
    const std::int64_t next_y = checked_coordinate_sum(block.occupied_sum_y, site.y);
    const std::int64_t next_z = checked_coordinate_sum(block.occupied_sum_z, site.z);
    ++block.occupied_voxel_count;
    block.occupied_sum_x = next_x;
    block.occupied_sum_y = next_y;
    block.occupied_sum_z = next_z;
    topology_dirty_ = topology_dirty_ ||
                      qualified_before != qualifies_as_core(coordinate, block);
    statistics_dirty_ = true;
}

void LesionIndex3D::remove_occupied_site(Vec3i site) {
    const Vec3i coordinate = block_coordinate(site);
    auto found = blocks_.find(coordinate);
    if (found == blocks_.end() || found->second.occupied_voxel_count == 0) {
        throw std::logic_error("lesion occupied-site removal underflow");
    }
    const bool qualified_before = qualifies_as_core(coordinate, found->second);
    BlockAggregate& block = found->second;
    const std::int64_t next_x = checked_coordinate_sum(
        block.occupied_sum_x, -static_cast<std::int64_t>(site.x));
    const std::int64_t next_y = checked_coordinate_sum(
        block.occupied_sum_y, -static_cast<std::int64_t>(site.y));
    const std::int64_t next_z = checked_coordinate_sum(
        block.occupied_sum_z, -static_cast<std::int64_t>(site.z));
    --block.occupied_voxel_count;
    block.occupied_sum_x = next_x;
    block.occupied_sum_y = next_y;
    block.occupied_sum_z = next_z;
    erase_block_if_empty(coordinate);
    const auto after = blocks_.find(coordinate);
    const bool qualified_after = after != blocks_.end() &&
                                 qualifies_as_core(coordinate, after->second);
    topology_dirty_ = topology_dirty_ || qualified_before != qualified_after;
    statistics_dirty_ = true;
}

void LesionIndex3D::rebuild_from(const CellStore3D& cells,
                                 const SparseChunkGrid3D& grid,
                                 LesionBiologicalVolumes3D biological_volumes) {
    const std::array<double, 3> stage_volumes{{
        biological_volumes.stage0_large,
        biological_volumes.stage1_small,
        biological_volumes.stage2_ultrasmall,
    }};
    if (std::any_of(stage_volumes.begin(), stage_volumes.end(), [](double value) {
            return !std::isfinite(value) || value <= 0.0;
        })) {
        throw std::invalid_argument("lesion stage biological volumes must be finite and positive");
    }
    begin_rebuild();
    for (std::size_t index = 0; index < cells.slot_count(); ++index) {
        const Slot slot = static_cast<Slot>(index);
        if (!cells.valid(slot)) continue;
        const std::size_t stage = static_cast<std::size_t>(cells.stage(slot));
        if (stage >= stage_volumes.size()) {
            throw std::logic_error("cell has invalid stage during lesion rebuild");
        }
        add_cell_anchor(cells.anchor(slot), stage_volumes[stage]);
    }
    grid.for_each_occupied_site(
        [this](Vec3i site, Slot) { add_occupied_site(site); });
    refresh_topology();
}

std::size_t LesionIndex3D::mark_dirty_sites(std::span<const Vec3i> changed_sites) {
    const std::size_t before = dirty_blocks_.size();
    dirty_blocks_.reserve(dirty_blocks_.size() + changed_sites.size());
    for (const Vec3i site : changed_sites) {
        dirty_blocks_.insert(block_coordinate(site));
    }
    return dirty_blocks_.size() - before;
}

std::size_t LesionIndex3D::refresh_dirty_blocks_from(
    const CellStore3D& cells,
    const SparseChunkGrid3D& grid,
    LesionBiologicalVolumes3D biological_volumes) {
    const std::array<double, 3> stage_volumes{{
        biological_volumes.stage0_large,
        biological_volumes.stage1_small,
        biological_volumes.stage2_ultrasmall,
    }};
    if (std::any_of(stage_volumes.begin(), stage_volumes.end(), [](double value) {
            return !std::isfinite(value) || value <= 0.0;
        })) {
        throw std::invalid_argument("lesion stage biological volumes must be finite and positive");
    }

    std::vector<Vec3i> coordinates(dirty_blocks_.begin(), dirty_blocks_.end());
    std::sort(coordinates.begin(), coordinates.end(), coordinate_less);
    for (const Vec3i coordinate : coordinates) {
        const auto old = blocks_.find(coordinate);
        const bool qualified_before = old != blocks_.end() &&
                                      qualifies_as_core(coordinate, old->second);

        BlockAggregate replacement;
        std::unordered_set<Slot> unique_slots;
        const Vec3i minimum = block_minimum_site(coordinate, config_.block_edge);
        const Vec3i maximum = block_maximum_site(coordinate, config_.block_edge);
        for (std::int64_t z = minimum.z; z <= maximum.z; ++z) {
            for (std::int64_t y = minimum.y; y <= maximum.y; ++y) {
                for (std::int64_t x = minimum.x; x <= maximum.x; ++x) {
                    const Vec3i site{static_cast<std::int32_t>(x),
                                     static_cast<std::int32_t>(y),
                                     static_cast<std::int32_t>(z)};
                    if (grid.owner(site) == kEmptySlot) continue;
                    ++replacement.occupied_voxel_count;
                    replacement.occupied_sum_x = checked_coordinate_sum(
                        replacement.occupied_sum_x, site.x);
                    replacement.occupied_sum_y = checked_coordinate_sum(
                        replacement.occupied_sum_y, site.y);
                    replacement.occupied_sum_z = checked_coordinate_sum(
                        replacement.occupied_sum_z, site.z);
                    for (const Slot slot : grid.occupants(site)) {
                        if (!cells.valid(slot) ||
                            block_coordinate(cells.anchor(slot)) != coordinate) {
                            continue;
                        }
                        unique_slots.insert(slot);
                    }
                }
            }
        }

        std::vector<Slot> ordered_slots(unique_slots.begin(), unique_slots.end());
        std::sort(ordered_slots.begin(), ordered_slots.end());
        for (const Slot slot : ordered_slots) {
            const Vec3i anchor = cells.anchor(slot);
            const std::size_t stage = static_cast<std::size_t>(cells.stage(slot));
            if (stage >= stage_volumes.size()) {
                throw std::logic_error("cell has invalid stage during lesion block refresh");
            }
            ++replacement.cell_count;
            replacement.biological_volume += stage_volumes[stage];
            replacement.cell_sum_x = checked_coordinate_sum(
                replacement.cell_sum_x, anchor.x);
            replacement.cell_sum_y = checked_coordinate_sum(
                replacement.cell_sum_y, anchor.y);
            replacement.cell_sum_z = checked_coordinate_sum(
                replacement.cell_sum_z, anchor.z);
        }

        if (replacement.cell_count == 0 && replacement.occupied_voxel_count == 0) {
            blocks_.erase(coordinate);
        } else {
            blocks_.insert_or_assign(coordinate, replacement);
        }
        const auto after = blocks_.find(coordinate);
        const bool qualified_after = after != blocks_.end() &&
                                     qualifies_as_core(coordinate, after->second);
        topology_dirty_ = topology_dirty_ || qualified_before != qualified_after;
    }
    dirty_blocks_.clear();
    if (!coordinates.empty()) statistics_dirty_ = true;
    return coordinates.size();
}

void LesionIndex3D::erase_block_if_empty(Vec3i coordinate) {
    const auto found = blocks_.find(coordinate);
    if (found != blocks_.end() && found->second.cell_count == 0 &&
        found->second.occupied_voxel_count == 0) {
        blocks_.erase(found);
    }
}

bool LesionIndex3D::qualifies_as_core(Vec3i coordinate,
                                      const BlockAggregate& block) const noexcept {
    if (block.cell_count < config_.minimum_cells_per_core_block ||
        block.biological_volume < config_.minimum_biological_volume_per_core_block) {
        return false;
    }
    const bool was_core = previous_core_identity_.contains(coordinate);
    const double threshold = was_core
        ? config_.core_deactivation_occupied_fraction
        : config_.core_activation_occupied_fraction;
    const double occupied_fraction = static_cast<double>(block.occupied_voxel_count) /
                                     static_cast<double>(block_voxel_capacity_);
    return occupied_fraction + 1e-15 >= threshold;
}

LesionId LesionIndex3D::allocate_lesion_id() {
    if (next_lesion_id_ == kNoLesionId) {
        throw std::overflow_error("lesion ID space exhausted");
    }
    return next_lesion_id_++;
}

LesionTopologyDelta3D LesionIndex3D::refresh_topology() {
    if (observations_dirty()) {
        throw std::logic_error(
            "lesion topology refresh requires dirty observations to be scanned first");
    }
    std::vector<Vec3i> core_coordinates;
    core_coordinates.reserve(blocks_.size());
    for (const auto& [coordinate, block] : blocks_) {
        if (qualifies_as_core(coordinate, block)) {
            core_coordinates.push_back(coordinate);
        }
    }
    std::sort(core_coordinates.begin(), core_coordinates.end(), coordinate_less);

    const std::unordered_set<Vec3i, Vec3iHash> core_set(core_coordinates.begin(),
                                                       core_coordinates.end());
    std::unordered_set<Vec3i, Vec3iHash> visited;
    visited.reserve(core_coordinates.size());
    const std::vector<Vec3i> offsets = connectivity_offsets(config_.connectivity);
    std::vector<std::vector<Vec3i>> components;

    for (const Vec3i start : core_coordinates) {
        if (!visited.insert(start).second) continue;
        std::vector<Vec3i> component;
        std::queue<Vec3i> pending;
        pending.push(start);
        while (!pending.empty()) {
            const Vec3i current = pending.front();
            pending.pop();
            component.push_back(current);
            for (const Vec3i offset : offsets) {
                bool valid = false;
                const Vec3i neighbor = checked_offset(
                    current, offset.x, offset.y, offset.z, valid);
                if (!valid || !core_set.contains(neighbor) ||
                    !visited.insert(neighbor).second) {
                    continue;
                }
                pending.push(neighbor);
            }
        }
        std::sort(component.begin(), component.end(), coordinate_less);
        components.push_back(std::move(component));
    }

    // Build every possible previous-ID overlap edge, then greedily match the
    // globally greatest overlaps. Component size and lexicographic minimum
    // coordinate provide deterministic split tie-breaking; the smaller old ID
    // wins a deterministic merge tie.
    struct MatchEdge {
        std::size_t overlap{};
        std::size_t component_index{};
        std::size_t component_size{};
        Vec3i component_minimum{};
        LesionId old_id{};
    };
    std::vector<MatchEdge> match_edges;
    std::vector<std::vector<LesionId>> component_predecessors(components.size());
    for (std::size_t component_index = 0; component_index < components.size();
         ++component_index) {
        std::unordered_map<LesionId, std::size_t> overlaps;
        for (const Vec3i coordinate : components[component_index]) {
            const auto old = previous_core_identity_.find(coordinate);
            if (old != previous_core_identity_.end()) {
                ++overlaps[old->second];
            }
        }
        for (const auto& [old_id, overlap] : overlaps) {
            match_edges.push_back({overlap, component_index,
                                   components[component_index].size(),
                                   components[component_index].front(), old_id});
            component_predecessors[component_index].push_back(old_id);
        }
        std::sort(component_predecessors[component_index].begin(),
                  component_predecessors[component_index].end());
    }
    std::sort(match_edges.begin(), match_edges.end(), [](const MatchEdge& lhs,
                                                         const MatchEdge& rhs) {
        if (lhs.overlap != rhs.overlap) return lhs.overlap > rhs.overlap;
        if (lhs.component_size != rhs.component_size) {
            return lhs.component_size > rhs.component_size;
        }
        if (lhs.component_minimum != rhs.component_minimum) {
            return coordinate_less(lhs.component_minimum, rhs.component_minimum);
        }
        return lhs.old_id < rhs.old_id;
    });

    std::vector<LesionId> component_ids(components.size(), kNoLesionId);
    std::unordered_set<LesionId> claimed_old_ids;
    for (const MatchEdge& edge : match_edges) {
        if (component_ids[edge.component_index] != kNoLesionId ||
            claimed_old_ids.contains(edge.old_id)) {
            continue;
        }
        component_ids[edge.component_index] = edge.old_id;
        claimed_old_ids.insert(edge.old_id);
    }
    for (LesionId& id : component_ids) {
        if (id == kNoLesionId) id = allocate_lesion_id();
    }

    LesionTopologyDelta3D delta;
    std::set<LesionId> old_ids;
    for (const auto& [coordinate, id] : previous_core_identity_) {
        (void)coordinate;
        old_ids.insert(id);
    }
    const std::set<LesionId> new_ids(component_ids.begin(), component_ids.end());
    std::set_difference(new_ids.begin(), new_ids.end(), old_ids.begin(), old_ids.end(),
                        std::back_inserter(delta.created));
    std::set_difference(old_ids.begin(), old_ids.end(), new_ids.begin(), new_ids.end(),
                        std::back_inserter(delta.removed));

    std::map<LesionId, std::vector<LesionId>> predecessor_children;
    for (std::size_t component_index = 0; component_index < components.size();
         ++component_index) {
        const std::vector<LesionId>& predecessors =
            component_predecessors[component_index];
        if (predecessors.size() > 1) {
            delta.merges.push_back({component_ids[component_index], predecessors});
        }
        for (const LesionId predecessor : predecessors) {
            predecessor_children[predecessor].push_back(component_ids[component_index]);
        }
    }
    for (auto& [predecessor, children] : predecessor_children) {
        std::sort(children.begin(), children.end());
        children.erase(std::unique(children.begin(), children.end()), children.end());
        if (children.size() > 1) {
            const LesionId retained =
                std::binary_search(children.begin(), children.end(), predecessor)
                    ? predecessor
                    : kNoLesionId;
            delta.splits.push_back({predecessor, retained, children});
        }
    }
    std::sort(delta.merges.begin(), delta.merges.end(), [](const LesionMerge3D& lhs,
                                                           const LesionMerge3D& rhs) {
        return lhs.result < rhs.result;
    });

    core_lesion_by_block_.clear();
    core_lesion_by_block_.reserve(core_coordinates.size());
    for (std::size_t component_index = 0; component_index < components.size();
         ++component_index) {
        for (const Vec3i coordinate : components[component_index]) {
            core_lesion_by_block_.emplace(coordinate, component_ids[component_index]);
        }
    }

    assigned_lesion_by_block_ = core_lesion_by_block_;
    if (config_.halo_blocks > 0) {
        std::vector<Vec3i> non_core_coordinates;
        non_core_coordinates.reserve(blocks_.size());
        for (const auto& [coordinate, block] : blocks_) {
            if (!core_lesion_by_block_.contains(coordinate) &&
                (block.cell_count != 0 || block.occupied_voxel_count != 0)) {
                non_core_coordinates.push_back(coordinate);
            }
        }
        std::sort(non_core_coordinates.begin(), non_core_coordinates.end(), coordinate_less);
        for (const Vec3i coordinate : non_core_coordinates) {
            int best_distance = std::numeric_limits<int>::max();
            LesionId best_id = kNoLesionId;
            for (int dx = -config_.halo_blocks; dx <= config_.halo_blocks; ++dx) {
                for (int dy = -config_.halo_blocks; dy <= config_.halo_blocks; ++dy) {
                    for (int dz = -config_.halo_blocks; dz <= config_.halo_blocks; ++dz) {
                        const int distance = std::max({std::abs(dx), std::abs(dy),
                                                       std::abs(dz)});
                        if (distance == 0 || distance > best_distance) continue;
                        bool valid = false;
                        const Vec3i candidate = checked_offset(
                            coordinate, dx, dy, dz, valid);
                        if (!valid) continue;
                        const auto core = core_lesion_by_block_.find(candidate);
                        if (core == core_lesion_by_block_.end()) continue;
                        if (distance < best_distance ||
                            (distance == best_distance && core->second < best_id)) {
                            best_distance = distance;
                            best_id = core->second;
                        }
                    }
                }
            }
            if (best_id != kNoLesionId) {
                assigned_lesion_by_block_.emplace(coordinate, best_id);
            }
        }
    }

    std::map<LesionId, LesionSummary3D> summaries;
    for (std::size_t component_index = 0; component_index < components.size();
         ++component_index) {
        LesionSummary3D& summary = summaries[component_ids[component_index]];
        summary.id = component_ids[component_index];
        summary.core_blocks = components[component_index];
    }

    std::vector<std::pair<Vec3i, LesionId>> assigned;
    assigned.reserve(assigned_lesion_by_block_.size());
    for (const auto& entry : assigned_lesion_by_block_) assigned.push_back(entry);
    std::sort(assigned.begin(), assigned.end(), [](const auto& lhs, const auto& rhs) {
        return coordinate_less(lhs.first, rhs.first);
    });

    struct CentroidAccumulator {
        long double cell_x{};
        long double cell_y{};
        long double cell_z{};
        long double occupied_x{};
        long double occupied_y{};
        long double occupied_z{};
    };
    std::map<LesionId, CentroidAccumulator> centroid_sums;
    std::unordered_set<LesionId> bounds_initialized;
    for (const auto& [coordinate, id] : assigned) {
        const auto block_found = blocks_.find(coordinate);
        if (block_found == blocks_.end()) continue;
        const BlockAggregate& block = block_found->second;
        LesionSummary3D& summary = summaries.at(id);
        summary.cell_count += block.cell_count;
        summary.occupied_voxel_count += block.occupied_voxel_count;
        summary.biological_volume += block.biological_volume;
        ++summary.assigned_block_count;
        CentroidAccumulator& sums = centroid_sums[id];
        sums.cell_x += block.cell_sum_x;
        sums.cell_y += block.cell_sum_y;
        sums.cell_z += block.cell_sum_z;
        sums.occupied_x += block.occupied_sum_x;
        sums.occupied_y += block.occupied_sum_y;
        sums.occupied_z += block.occupied_sum_z;

        const Vec3i minimum = block_minimum_site(coordinate, config_.block_edge);
        const Vec3i maximum = block_maximum_site(coordinate, config_.block_edge);
        if (bounds_initialized.insert(id).second) {
            summary.minimum_site = minimum;
            summary.maximum_site = maximum;
        } else {
            summary.minimum_site.x = std::min(summary.minimum_site.x, minimum.x);
            summary.minimum_site.y = std::min(summary.minimum_site.y, minimum.y);
            summary.minimum_site.z = std::min(summary.minimum_site.z, minimum.z);
            summary.maximum_site.x = std::max(summary.maximum_site.x, maximum.x);
            summary.maximum_site.y = std::max(summary.maximum_site.y, maximum.y);
            summary.maximum_site.z = std::max(summary.maximum_site.z, maximum.z);
        }
    }

    lesions_.clear();
    lesions_.reserve(summaries.size());
    for (auto& [id, summary] : summaries) {
        const CentroidAccumulator& sums = centroid_sums[id];
        if (summary.occupied_voxel_count != 0) {
            const long double denominator = summary.occupied_voxel_count;
            summary.centroid = {static_cast<double>(sums.occupied_x / denominator),
                                static_cast<double>(sums.occupied_y / denominator),
                                static_cast<double>(sums.occupied_z / denominator)};
        } else if (summary.cell_count != 0) {
            const long double denominator = summary.cell_count;
            summary.centroid = {static_cast<double>(sums.cell_x / denominator),
                                static_cast<double>(sums.cell_y / denominator),
                                static_cast<double>(sums.cell_z / denominator)};
        }
        lesions_.push_back(std::move(summary));
    }

    previous_core_identity_ = core_lesion_by_block_;
    topology_dirty_ = false;
    statistics_dirty_ = false;
    return delta;
}

std::optional<LesionBlockStats3D> LesionIndex3D::block_stats(
    Vec3i coordinate) const {
    const auto found = blocks_.find(coordinate);
    if (found == blocks_.end()) return std::nullopt;
    const auto assigned = assigned_lesion_by_block_.find(coordinate);
    return LesionBlockStats3D{
        found->second.cell_count,
        found->second.occupied_voxel_count,
        found->second.biological_volume,
        core_lesion_by_block_.contains(coordinate),
        assigned == assigned_lesion_by_block_.end() ? kNoLesionId
                                                    : assigned->second,
    };
}

bool LesionIndex3D::is_core_block(Vec3i coordinate) const noexcept {
    return core_lesion_by_block_.contains(coordinate);
}

const LesionSummary3D* LesionIndex3D::find_lesion(LesionId id) const noexcept {
    const auto found = std::lower_bound(
        lesions_.begin(), lesions_.end(), id,
        [](const LesionSummary3D& lesion, LesionId value) { return lesion.id < value; });
    return found != lesions_.end() && found->id == id ? &*found : nullptr;
}

std::optional<LesionId> LesionIndex3D::lesion_at_site(Vec3i occupied_site) const noexcept {
    const auto found = assigned_lesion_by_block_.find(block_coordinate(occupied_site));
    if (found == assigned_lesion_by_block_.end()) return std::nullopt;
    return found->second;
}

std::optional<LesionId> LesionIndex3D::lesion_for_face(
    const ExposedFace3D& face) const noexcept {
    return lesion_at_site(face.inside);
}

std::optional<LesionId> LesionIndex3D::lesion_for_face(
    const ExposedFace3D& face,
    const CellStore3D& cells,
    const SparseChunkGrid3D& grid) const noexcept {
    const Slot owner = grid.owner(face.inside);
    if (owner == kEmptySlot || !cells.valid(owner)) return std::nullopt;
    return lesion_for_anchor(cells.anchor(owner));
}

bool LesionIndex3D::lesion_owns_site(LesionId id, Vec3i occupied_site) const noexcept {
    const std::optional<LesionId> owner = lesion_at_site(occupied_site);
    return id != kNoLesionId && owner.has_value() && *owner == id;
}

bool LesionIndex3D::lesion_owns_face(LesionId id,
                                     const ExposedFace3D& face) const noexcept {
    return lesion_owns_site(id, face.inside);
}

void LesionIndex3D::set_next_lesion_id(LesionId value) {
    LesionId maximum_id = kNoLesionId;
    for (const LesionSummary3D& lesion : lesions_) {
        maximum_id = std::max(maximum_id, lesion.id);
    }
    for (const auto& [coordinate, id] : previous_core_identity_) {
        (void)coordinate;
        maximum_id = std::max(maximum_id, id);
    }
    if (value == kNoLesionId || value <= maximum_id || value < next_lesion_id_) {
        throw std::invalid_argument(
            "next lesion ID must not decrease and must exceed all existing IDs");
    }
    next_lesion_id_ = value;
}

std::vector<LesionCoreIdentity3D> LesionIndex3D::snapshot_core_identity() const {
    std::vector<LesionCoreIdentity3D> result;
    result.reserve(core_lesion_by_block_.size());
    for (const auto& [block, lesion_id] : core_lesion_by_block_) {
        result.push_back({block, lesion_id});
    }
    std::sort(result.begin(), result.end(), [](const LesionCoreIdentity3D& lhs,
                                               const LesionCoreIdentity3D& rhs) {
        return coordinate_less(lhs.block, rhs.block);
    });
    return result;
}

void LesionIndex3D::restore_core_identity(
    std::span<const LesionCoreIdentity3D> identity,
    LesionId next_id) {
    if (observations_dirty()) {
        throw std::logic_error(
            "cannot restore lesion identity while observation blocks are dirty");
    }
    IdentityMap restored;
    restored.reserve(identity.size());
    LesionId maximum_id = kNoLesionId;
    for (const LesionCoreIdentity3D& entry : identity) {
        if (entry.lesion_id == kNoLesionId) {
            throw std::invalid_argument("checkpoint lesion identity contains zero ID");
        }
        if (!restored.emplace(entry.block, entry.lesion_id).second) {
            throw std::invalid_argument("checkpoint lesion identity contains duplicate block");
        }
        const auto observed = blocks_.find(entry.block);
        if (observed == blocks_.end() ||
            observed->second.cell_count < config_.minimum_cells_per_core_block ||
            observed->second.biological_volume <
                config_.minimum_biological_volume_per_core_block) {
            throw std::invalid_argument(
                "checkpoint lesion core block has no qualifying observation");
        }
        const double occupied_fraction =
            static_cast<double>(observed->second.occupied_voxel_count) /
            static_cast<double>(block_voxel_capacity_);
        if (occupied_fraction + 1e-15 <
            config_.core_deactivation_occupied_fraction) {
            throw std::invalid_argument(
                "checkpoint lesion core block is below deactivation threshold");
        }
        maximum_id = std::max(maximum_id, entry.lesion_id);
    }
    if (next_id == kNoLesionId || next_id <= maximum_id) {
        throw std::invalid_argument(
            "checkpoint next lesion ID must exceed every restored ID");
    }

    previous_core_identity_ = std::move(restored);
    core_lesion_by_block_.clear();
    assigned_lesion_by_block_.clear();
    lesions_.clear();
    next_lesion_id_ = next_id;
    topology_dirty_ = true;
    statistics_dirty_ = true;
    refresh_topology();

    if (core_lesion_by_block_.size() != identity.size()) {
        throw std::invalid_argument(
            "checkpoint lesion identity does not match rebuilt topology");
    }
    for (const LesionCoreIdentity3D& entry : identity) {
        const auto rebuilt = core_lesion_by_block_.find(entry.block);
        if (rebuilt == core_lesion_by_block_.end() ||
            rebuilt->second != entry.lesion_id) {
            throw std::invalid_argument(
                "checkpoint lesion ID continuation failed topology validation");
        }
    }
}

std::vector<LesionDirtyBlockState3D>
LesionIndex3D::snapshot_dirty_block_state() const {
    std::vector<Vec3i> coordinates(dirty_blocks_.begin(), dirty_blocks_.end());
    std::sort(coordinates.begin(), coordinates.end(), coordinate_less);
    std::vector<LesionDirtyBlockState3D> result;
    result.reserve(coordinates.size());
    for (const Vec3i coordinate : coordinates) {
        const auto observed = blocks_.find(coordinate);
        if (observed == blocks_.end()) {
            result.push_back({.block = coordinate, .exists = false});
            continue;
        }
        const BlockAggregate& block = observed->second;
        result.push_back({
            .block = coordinate,
            .exists = true,
            .cell_count = block.cell_count,
            .occupied_voxel_count = block.occupied_voxel_count,
            .biological_volume = block.biological_volume,
            .cell_coordinate_sum_x = block.cell_sum_x,
            .cell_coordinate_sum_y = block.cell_sum_y,
            .cell_coordinate_sum_z = block.cell_sum_z,
            .occupied_coordinate_sum_x = block.occupied_sum_x,
            .occupied_coordinate_sum_y = block.occupied_sum_y,
            .occupied_coordinate_sum_z = block.occupied_sum_z,
        });
    }
    return result;
}

void LesionIndex3D::restore_checkpoint_state(
    std::span<const LesionCoreIdentity3D> core_identity,
    LesionId next_id,
    std::span<const LesionDirtyBlockState3D> dirty_states) {
    if (observations_dirty()) {
        throw std::logic_error(
            "checkpoint lesion restore requires a clean current observation rebuild");
    }

    std::unordered_set<Vec3i, Vec3iHash> coordinates;
    coordinates.reserve(dirty_states.size());
    for (const LesionDirtyBlockState3D& state : dirty_states) {
        if (!coordinates.insert(state.block).second) {
            throw std::invalid_argument(
                "checkpoint lesion dirty state contains duplicate block");
        }
        if (!std::isfinite(state.biological_volume) ||
            state.biological_volume < 0.0) {
            throw std::invalid_argument(
                "checkpoint lesion dirty biological volume must be finite and nonnegative");
        }

        const bool payload_is_zero =
            state.cell_count == 0 && state.occupied_voxel_count == 0 &&
            state.biological_volume == 0.0 &&
            state.cell_coordinate_sum_x == 0 &&
            state.cell_coordinate_sum_y == 0 &&
            state.cell_coordinate_sum_z == 0 &&
            state.occupied_coordinate_sum_x == 0 &&
            state.occupied_coordinate_sum_y == 0 &&
            state.occupied_coordinate_sum_z == 0;
        if (!state.exists) {
            if (!payload_is_zero) {
                throw std::invalid_argument(
                    "absent checkpoint lesion dirty block must have an empty payload");
            }
            continue;
        }
        if (state.cell_count == 0 && state.occupied_voxel_count == 0) {
            throw std::invalid_argument(
                "present checkpoint lesion dirty block must not be empty");
        }
        if (state.cell_count > std::numeric_limits<Slot>::max()) {
            throw std::invalid_argument(
                "checkpoint lesion dirty cell count exceeds stable slot capacity");
        }
        if (state.occupied_voxel_count > block_voxel_capacity_) {
            throw std::invalid_argument(
                "checkpoint lesion dirty occupancy exceeds block capacity");
        }
        if ((state.cell_count == 0 && state.biological_volume != 0.0) ||
            (state.cell_count != 0 && state.biological_volume <= 0.0)) {
            throw std::invalid_argument(
                "checkpoint lesion dirty cell count and biological volume disagree");
        }

        const Vec3i minimum = block_minimum_site(state.block, config_.block_edge);
        const Vec3i maximum = block_maximum_site(state.block, config_.block_edge);
        validate_sum_within_block(state.cell_coordinate_sum_x, state.cell_count,
                                  minimum.x, maximum.x,
                                  "checkpoint lesion cell x sum");
        validate_sum_within_block(state.cell_coordinate_sum_y, state.cell_count,
                                  minimum.y, maximum.y,
                                  "checkpoint lesion cell y sum");
        validate_sum_within_block(state.cell_coordinate_sum_z, state.cell_count,
                                  minimum.z, maximum.z,
                                  "checkpoint lesion cell z sum");
        validate_sum_within_block(state.occupied_coordinate_sum_x,
                                  state.occupied_voxel_count,
                                  minimum.x, maximum.x,
                                  "checkpoint lesion occupied x sum");
        validate_sum_within_block(state.occupied_coordinate_sum_y,
                                  state.occupied_voxel_count,
                                  minimum.y, maximum.y,
                                  "checkpoint lesion occupied y sum");
        validate_sum_within_block(state.occupied_coordinate_sum_z,
                                  state.occupied_voxel_count,
                                  minimum.z, maximum.z,
                                  "checkpoint lesion occupied z sum");
    }

    for (const LesionDirtyBlockState3D& state : dirty_states) {
        if (!state.exists) {
            blocks_.erase(state.block);
            continue;
        }
        blocks_.insert_or_assign(state.block, BlockAggregate{
            .cell_count = state.cell_count,
            .occupied_voxel_count = state.occupied_voxel_count,
            .biological_volume = state.biological_volume,
            .cell_sum_x = state.cell_coordinate_sum_x,
            .cell_sum_y = state.cell_coordinate_sum_y,
            .cell_sum_z = state.cell_coordinate_sum_z,
            .occupied_sum_x = state.occupied_coordinate_sum_x,
            .occupied_sum_y = state.occupied_coordinate_sum_y,
            .occupied_sum_z = state.occupied_coordinate_sum_z,
        });
    }

    restore_core_identity(core_identity, next_id);

    dirty_blocks_.reserve(dirty_states.size());
    for (const LesionDirtyBlockState3D& state : dirty_states) {
        dirty_blocks_.insert(state.block);
    }
}

std::size_t LesionIndex3D::allocated_bytes() const noexcept {
    const auto map_bytes = [](const auto& map) {
        return map.bucket_count() * sizeof(void*) +
               map.size() * (sizeof(typename std::decay_t<decltype(map)>::value_type) +
                             2 * sizeof(void*));
    };
    std::size_t bytes = map_bytes(blocks_) + map_bytes(previous_core_identity_) +
                        map_bytes(core_lesion_by_block_) +
                        map_bytes(assigned_lesion_by_block_) +
                        lesions_.capacity() * sizeof(LesionSummary3D);
    bytes += dirty_blocks_.bucket_count() * sizeof(void*) +
             dirty_blocks_.size() * (sizeof(Vec3i) + 2 * sizeof(void*));
    for (const LesionSummary3D& lesion : lesions_) {
        bytes += lesion.core_blocks.capacity() * sizeof(Vec3i);
    }
    return bytes;
}

}  // namespace atcg3d
