#pragma once

#include <cstddef>
#include <cstdint>
#include <unordered_map>
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

    void add(Vec3i anchor, CellType type);
    void remove(Vec3i anchor, CellType type);
    void move(Vec3i from, Vec3i to, CellType type);

    DensityCounts3D block_counts(Vec3i block_coordinate) const;
    DensityCounts3D estimate_box(Vec3i minimum, Vec3i maximum) const;
    DensityCounts3D estimate_quantized_box(Vec3i anchor,
                                           int window_edge,
                                           int query_block_edge) const;
    double estimate_directional_density(Vec3i anchor,
                                        DirectionId direction,
                                        int radius,
                                        double half_angle_degrees) const;

private:
    struct AnchorEntry {
        std::uint8_t x{};
        std::uint8_t y{};
        std::uint8_t z{};
        std::uint8_t type{};
    };

    struct Block {
        DensityCounts3D counts;
        std::vector<AnchorEntry> anchors;
    };

    struct CachedBox {
        std::uint64_t generation{};
        DensityCounts3D counts;
    };

    static int floor_div(int value, int divisor) noexcept;
    Vec3i block_coordinate(Vec3i site) const noexcept;
    AnchorEntry local_entry(Vec3i site, CellType type) const noexcept;

    int block_edge_{};
    std::unordered_map<Vec3i, Block, Vec3iHash> counts_;
    std::uint64_t generation_{1};
    mutable std::unordered_map<Vec3i, CachedBox, Vec3iHash> quantized_box_cache_;
};

}  // namespace atcg3d
