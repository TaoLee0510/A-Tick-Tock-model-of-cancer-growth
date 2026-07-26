#pragma once

#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <shared_mutex>
#include <unordered_map>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace atcg3d {

class QuantizedBoxCountIndex3D;

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

// Execution-only diagnostics for the most recent anchor move while exact
// per-slot growth-window counts were enabled. The queried lattice-site
// volumes describe the disjoint box regions presented to the block index; they
// do not affect biology, checkpoint state, or deterministic simulation output.
struct LocalWindowMoveDiagnostics3D {
    std::uint64_t affected_window_query_sites{};
    std::uint64_t moving_window_query_sites{};
    std::uint64_t affected_slot_visits{};
};

class BlockDensityIndex3D {
public:
    explicit BlockDensityIndex3D(int block_edge);

    int block_edge() const noexcept { return block_edge_; }
    std::size_t block_count() const noexcept { return counts_.size(); }
    std::size_t allocated_bytes() const;
    std::size_t quantized_cache_size() const;
    static constexpr std::size_t quantized_cache_capacity() noexcept { return 65536; }

    void add(Vec3i anchor, CellType type, Slot slot);
    void remove(Vec3i anchor, CellType type, Slot slot);
    void move(Vec3i from, Vec3i to, CellType type, Slot slot);
    void attach_quantized_count_index(QuantizedBoxCountIndex3D* index) noexcept {
        incremental_quantized_index_ = index;
    }
    void configure_local_window_counts(int window_edge, bool thin_layer);
    void begin_local_window_bulk_load();
    void finish_local_window_bulk_load(std::size_t slot_count, int workers);
    bool has_local_window_counts(int window_edge, bool thin_layer) const noexcept {
        return local_window_edge_ == window_edge &&
               local_window_thin_layer_ == thin_layer;
    }
    DensityCounts3D local_window_counts(Slot slot) const;
    const LocalWindowMoveDiagnostics3D&
    last_local_window_move_diagnostics() const noexcept {
        return last_local_window_move_diagnostics_;
    }
    void reset_quantized_cache();

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
    // Evaluates all 26 lattice cones while visiting nearby anchors only once.
    // Element zero is unused and remains zero.
    std::array<double, 27> estimate_all_directional_densities(
        Vec3i anchor,
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
    void register_cache_layout_locked(CacheLayout layout) const;
    void adjust_local_window_counts(Vec3i changed_anchor,
                                    CellType changed_type,
                                    bool add);
    void set_local_window_counts(Slot slot, DensityCounts3D counts);

    int block_edge_{};
    std::unordered_map<Vec3i, Block, Vec3iHash> counts_;
    mutable std::vector<CacheLayout> cache_layouts_;
    mutable std::unordered_map<QuantizedCacheKey, DensityCounts3D,
                               QuantizedCacheKeyHash> quantized_box_cache_;
    mutable std::shared_mutex quantized_cache_mutex_;
    mutable std::atomic<bool> has_quantized_cache_layouts_{false};
    QuantizedBoxCountIndex3D* incremental_quantized_index_{};
    struct LocalCounts {
        std::uint32_t r{};
        std::uint32_t K{};
    };
    std::vector<LocalCounts> local_window_counts_;
    int local_window_edge_{};
    bool local_window_thin_layer_{};
    bool local_window_bulk_loading_{};
    LocalWindowMoveDiagnostics3D last_local_window_move_diagnostics_;
};

}  // namespace atcg3d
