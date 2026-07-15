#pragma once

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace atcg3d {

struct DensityCounts3D {
    std::uint64_t r{};
    std::uint64_t K{};

    std::uint64_t total() const noexcept { return r + K; }
    DensityCounts3D& operator+=(const DensityCounts3D& other) noexcept {
        r += other.r;
        K += other.K;
        return *this;
    }
};

class BlockDensityIndex3D {
public:
    explicit BlockDensityIndex3D(int block_edge);

    int block_edge() const noexcept { return block_edge_; }
    std::size_t block_count() const noexcept { return counts_.size(); }
    std::size_t allocated_bytes() const noexcept;
    std::size_t quantized_cache_size() const noexcept { return quantized_box_cache_.size(); }
    static constexpr std::size_t quantized_cache_capacity() noexcept { return 65536; }

    void add(Vec3i anchor, CellType type, Slot slot);
    void remove(Vec3i anchor, CellType type, Slot slot);
    void move(Vec3i from, Vec3i to, CellType type, Slot slot);

    template <class Visitor>
    void for_each_slot_in_block(Vec3i coordinate, Visitor&& visitor) const {
        const auto iterator = counts_.find(coordinate);
        if (iterator == counts_.end()) {
            return;
        }
        for (const AnchorEntry& entry : iterator->second.anchors) {
            visitor(entry.slot);
        }
    }

    template <class Visitor>
    void for_each_slot_in_box(Vec3i minimum, Vec3i maximum, Visitor&& visitor) const {
        if (minimum.x > maximum.x || minimum.y > maximum.y || minimum.z > maximum.z) {
            return;
        }
        const Vec3i first = block_coordinate(minimum);
        const Vec3i last = block_coordinate(maximum);
        for (int bx = first.x; bx <= last.x; ++bx) {
            for (int by = first.y; by <= last.y; ++by) {
                for (int bz = first.z; bz <= last.z; ++bz) {
                    const auto iterator = counts_.find({bx, by, bz});
                    if (iterator == counts_.end()) {
                        continue;
                    }
                    const Vec3i block_min{bx * block_edge_, by * block_edge_, bz * block_edge_};
                    for (const AnchorEntry& entry : iterator->second.anchors) {
                        const Vec3i anchor = block_min + Vec3i{entry.x, entry.y, entry.z};
                        if (anchor.x >= minimum.x && anchor.y >= minimum.y &&
                            anchor.z >= minimum.z && anchor.x <= maximum.x &&
                            anchor.y <= maximum.y && anchor.z <= maximum.z) {
                            visitor(entry.slot);
                        }
                    }
                }
            }
        }
    }

    DensityCounts3D block_counts(Vec3i block_coordinate) const;
    DensityCounts3D estimate_box(Vec3i minimum, Vec3i maximum) const;
    DensityCounts3D estimate_quantized_box(Vec3i anchor,
                                           int window_edge,
                                           int query_block_edge,
                                           bool thin_layer = false) const;
    double estimate_directional_density(Vec3i anchor,
                                        DirectionId direction,
                                        int radius,
                                        double half_angle_degrees,
                                        bool thin_layer = false) const;

private:
    struct AnchorEntry {
        std::uint8_t x{};
        std::uint8_t y{};
        std::uint8_t z{};
        std::uint8_t type{};
        Slot slot{kEmptySlot};
    };

    struct Block {
        DensityCounts3D counts;
        std::vector<AnchorEntry> anchors;
    };

    struct CacheLayout {
        int window_edge{};
        int query_block_edge{};
        bool thin_layer{};

        bool operator==(const CacheLayout&) const = default;
    };

    struct QuantizedCacheKey {
        Vec3i coordinate{};
        int window_edge{};
        int query_block_edge{};
        bool thin_layer{};

        bool operator==(const QuantizedCacheKey&) const = default;
    };

    struct QuantizedCacheKeyHash {
        std::size_t operator()(const QuantizedCacheKey& key) const noexcept;
    };

    static int floor_div(int value, int divisor) noexcept;
    Vec3i block_coordinate(Vec3i site) const noexcept;
    AnchorEntry local_entry(Vec3i site, CellType type, Slot slot) const noexcept;
    void invalidate_quantized_cache(Vec3i anchor);
    void register_cache_layout(CacheLayout layout) const;

    int block_edge_{};
    std::unordered_map<Vec3i, Block, Vec3iHash> counts_;
    mutable std::vector<CacheLayout> cache_layouts_;
    mutable std::unordered_map<QuantizedCacheKey, DensityCounts3D,
                               QuantizedCacheKeyHash> quantized_box_cache_;
};

}  // namespace atcg3d
