#pragma once

#include <cstdint>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace atcg3d {

// Exact incremental counts for a family of quantized axis-aligned density
// windows. A biological cell contributes once by anchor, irrespective of its
// footprint. Updating one anchor touches only the query blocks whose window
// contains that anchor (27 blocks for the legacy 70/32 mapping).
class QuantizedBoxCountIndex3D {
public:
    QuantizedBoxCountIndex3D(int window_edge,
                             int query_block_edge,
                             bool thin_layer);

    void add(Vec3i anchor);
    void remove(Vec3i anchor);
    void move(Vec3i from, Vec3i to);

    Vec3i query_block(Vec3i anchor) const noexcept;
    std::uint64_t count(Vec3i query_block) const noexcept;
    std::uint64_t resident_count(Vec3i query_block) const noexcept;
    std::vector<Vec3i> affected_blocks(Vec3i anchor) const;
    std::vector<Vec3i> resident_blocks() const;

    template <class Visitor>
    void for_each_affected_block(Vec3i anchor, Visitor&& visitor) const {
        const std::int64_t query_edge = query_block_edge_;
        const std::int64_t lower = (window_edge_ - 1) / 2;
        const std::int64_t upper = window_edge_ - lower - 1;
        const auto range = [&](std::int32_t value) {
            const std::int64_t shifted =
                static_cast<std::int64_t>(value) - query_edge / 2;
            return std::pair<std::int64_t, std::int64_t>{
                ceil_div(shifted - upper, query_edge),
                floor_div(shifted + lower, query_edge)};
        };
        const auto [minimum_x, maximum_x] = range(anchor.x);
        const auto [minimum_y, maximum_y] = range(anchor.y);
        const auto [minimum_z, maximum_z] = thin_layer_
            ? std::pair<std::int64_t, std::int64_t>{anchor.z, anchor.z}
            : range(anchor.z);
        for (std::int64_t x = minimum_x; x <= maximum_x; ++x) {
            for (std::int64_t y = minimum_y; y <= maximum_y; ++y) {
                for (std::int64_t z = minimum_z; z <= maximum_z; ++z) {
                    if (x < std::numeric_limits<std::int32_t>::min() ||
                        x > std::numeric_limits<std::int32_t>::max() ||
                        y < std::numeric_limits<std::int32_t>::min() ||
                        y > std::numeric_limits<std::int32_t>::max() ||
                        z < std::numeric_limits<std::int32_t>::min() ||
                        z > std::numeric_limits<std::int32_t>::max()) {
                        continue;
                    }
                    visitor(Vec3i{static_cast<std::int32_t>(x),
                                  static_cast<std::int32_t>(y),
                                  static_cast<std::int32_t>(z)});
                }
            }
        }
    }

    std::size_t block_count() const noexcept { return counts_.size(); }
    std::size_t allocated_bytes() const noexcept;

private:
    static std::int64_t floor_div(std::int64_t value,
                                  std::int64_t divisor) noexcept;
    static std::int64_t ceil_div(std::int64_t value,
                                 std::int64_t divisor) noexcept;

    int window_edge_{};
    int query_block_edge_{};
    bool thin_layer_{};
    std::unordered_map<Vec3i, std::uint64_t, Vec3iHash> counts_;
    std::unordered_map<Vec3i, std::uint64_t, Vec3iHash> residents_;
};

}  // namespace atcg3d
