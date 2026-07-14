#pragma once

#include <compare>
#include <cstdint>
#include <functional>
#include <limits>

namespace atcg3d {

using CellUid = std::uint64_t;
using Slot = std::uint32_t;
using DirectionId = std::uint8_t;

inline constexpr Slot kEmptySlot = std::numeric_limits<Slot>::max();
inline constexpr DirectionId kStayDirection = 0;

struct Vec3i {
    std::int32_t x{};
    std::int32_t y{};
    std::int32_t z{};

    constexpr auto operator<=>(const Vec3i&) const = default;

    constexpr Vec3i operator+(const Vec3i& rhs) const {
        return {static_cast<std::int32_t>(x + rhs.x),
                static_cast<std::int32_t>(y + rhs.y),
                static_cast<std::int32_t>(z + rhs.z)};
    }

    constexpr Vec3i operator-(const Vec3i& rhs) const {
        return {static_cast<std::int32_t>(x - rhs.x),
                static_cast<std::int32_t>(y - rhs.y),
                static_cast<std::int32_t>(z - rhs.z)};
    }
};

struct Vec3iHash {
    std::size_t operator()(const Vec3i& value) const noexcept {
        std::uint64_t h = static_cast<std::uint32_t>(value.x);
        h ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(value.y)) << 21U;
        h ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(value.z)) << 42U;
        h ^= h >> 30U;
        h *= 0xbf58476d1ce4e5b9ULL;
        h ^= h >> 27U;
        h *= 0x94d049bb133111ebULL;
        h ^= h >> 31U;
        return static_cast<std::size_t>(h);
    }
};

enum class CellType : std::uint8_t {
    r = 1,
    K = 2,
};

enum class CellStage : std::uint8_t {
    large = 0,
    small = 1,
    ultrasmall = 2,
};

enum CellFlags : std::uint8_t {
    kNoFlags = 0,
    kMigrationActive = 1U << 0U,
    kDirtyDensity = 1U << 1U,
};

}  // namespace atcg3d
