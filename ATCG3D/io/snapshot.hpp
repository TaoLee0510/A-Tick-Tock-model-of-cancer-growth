#pragma once

#include <array>
#include <span>

#include "core/cell_store.hpp"
#include "engine/simulation.hpp"

namespace atcg3d {

struct SimulationSnapshotView3D {
    const CellStore3D& cells;
    std::span<const Slot> slots;
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
};

}  // namespace atcg3d
