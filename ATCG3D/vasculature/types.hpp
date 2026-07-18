#pragma once

#include <cstdint>
#include <limits>

#include "core/types.hpp"

namespace atcg3d {

using VesselId = std::uint64_t;
using LesionId = std::uint64_t;
using VesselNodeUid = std::uint64_t;
using VesselTipUid = std::uint64_t;
using VesselNodeSlot = std::uint32_t;
using VesselTipSlot = std::uint32_t;

inline constexpr LesionId kNoLesionId = 0;

inline constexpr VesselNodeSlot kEmptyVesselNodeSlot =
    std::numeric_limits<VesselNodeSlot>::max();
inline constexpr VesselTipSlot kEmptyVesselTipSlot =
    std::numeric_limits<VesselTipSlot>::max();

enum class VesselBranchRole : std::uint8_t {
    root = 0,
    inward = 1,
    outward = 2,
};

enum class VesselTipStatus : std::uint8_t {
    dormant = 0,
    active = 1,
    blocked = 2,
    merged = 3,
    reached_target = 4,
    max_length = 5,
    boundary_stop = 6,
    complete = 7,
};

enum VesselVoxelBits : std::uint8_t {
    kVesselVoxelNone = 0,
    kVesselVoxelOccupied = 1U << 0U,
    kVesselVoxelPerfused = 1U << 1U,
    kVesselVoxelRoot = 1U << 2U,
    kVesselVoxelInward = 1U << 3U,
    kVesselVoxelOutward = 1U << 4U,
};

constexpr std::uint8_t vessel_role_bit(VesselBranchRole role) noexcept {
    switch (role) {
        case VesselBranchRole::root: return kVesselVoxelRoot;
        case VesselBranchRole::inward: return kVesselVoxelInward;
        case VesselBranchRole::outward: return kVesselVoxelOutward;
    }
    return kVesselVoxelNone;
}

constexpr bool vessel_tip_terminal(VesselTipStatus status) noexcept {
    return status != VesselTipStatus::dormant && status != VesselTipStatus::active;
}

}  // namespace atcg3d
