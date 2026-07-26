#pragma once

#include <array>
#include <cstddef>
#include <span>
#include <vector>

#include "core/cell_store.hpp"
#include "engine/simulation.hpp"

namespace atcg3d {

struct SimulationSnapshotView3D {
    const CellStore3D& cells;
    std::span<const Slot> slots;
    const LesionIndex3D& lesion_index;
    SimulationClock3D clock;
    DisplayRadiusConfig radii{};

    std::size_t size() const noexcept { return slots.size(); }
    Slot slot(std::size_t index) const { return slots[index]; }

    std::array<float, 3> center(std::size_t index) const {
        const Slot cell_slot = slot(index);
        const Vec3i anchor = cells.anchor(cell_slot);
        const float offset = cells.stage(cell_slot) == CellStage::large ? 1.0F : 0.5F;
        return {static_cast<float>(anchor.x) + offset,
                static_cast<float>(anchor.y) + offset,
                static_cast<float>(anchor.z) + offset};
    }

    float display_radius(std::size_t index) const {
        switch (cells.stage(slot(index))) {
            case CellStage::large: return static_cast<float>(radii.large);
            case CellStage::small: return static_cast<float>(radii.small);
            case CellStage::ultrasmall: return static_cast<float>(radii.ultrasmall);
        }
        return static_cast<float>(radii.small);
    }

    LesionId lesion_id(std::size_t index) const noexcept {
        return lesion_index.lesion_for_anchor(cells.anchor(slot(index)))
            .value_or(kNoLesionId);
    }

    std::size_t total_cell_count() const noexcept { return cells.alive_count(); }
    std::size_t total_slot_count() const noexcept { return cells.slot_count(); }
    CellUid uid(std::size_t index) const { return cells.uid(slot(index)); }
    std::uint32_t clone_id(std::size_t index) const {
        return cells.clone_id(slot(index));
    }
    CellType type(std::size_t index) const { return cells.type(slot(index)); }
    CellStage stage(std::size_t index) const { return cells.stage(slot(index)); }
    std::uint8_t viability(std::size_t index) const {
        return cells.viability(slot(index));
    }
};

// Read-only view over an owning asynchronous-output copy.  It intentionally
// contains no SparseChunkGrid3D, scheduler, or density index: rendering a
// frame must not rebuild the simulation's derived runtime state.
struct FrozenCellSnapshotView3D {
    std::span<const CellInit> cells;
    std::span<const Slot> cell_slots;
    std::size_t slot_count{};
    std::span<const std::size_t> indices;
    std::span<const LesionId> lesion_ids;
    SimulationClock3D clock;
    DisplayRadiusConfig radii{};

    std::size_t size() const noexcept { return indices.size(); }
    const CellInit& cell(std::size_t index) const { return cells[indices[index]]; }
    Slot slot(std::size_t index) const { return cell_slots[indices[index]]; }

    std::array<float, 3> center(std::size_t index) const {
        const CellInit& value = cell(index);
        const float offset = value.stage == CellStage::large ? 1.0F : 0.5F;
        return {static_cast<float>(value.anchor.x) + offset,
                static_cast<float>(value.anchor.y) + offset,
                static_cast<float>(value.anchor.z) + offset};
    }

    float display_radius(std::size_t index) const {
        switch (cell(index).stage) {
            case CellStage::large: return static_cast<float>(radii.large);
            case CellStage::small: return static_cast<float>(radii.small);
            case CellStage::ultrasmall: return static_cast<float>(radii.ultrasmall);
        }
        return static_cast<float>(radii.small);
    }

    LesionId lesion_id(std::size_t index) const noexcept {
        return lesion_ids[indices[index]];
    }
    std::size_t total_cell_count() const noexcept { return cells.size(); }
    std::size_t total_slot_count() const noexcept { return slot_count; }
    CellUid uid(std::size_t index) const { return cell(index).uid; }
    std::uint32_t clone_id(std::size_t index) const { return cell(index).clone_id; }
    CellType type(std::size_t index) const { return cell(index).type; }
    CellStage stage(std::size_t index) const { return cell(index).stage; }
    std::uint8_t viability(std::size_t index) const { return cell(index).viability; }
};

}  // namespace atcg3d
