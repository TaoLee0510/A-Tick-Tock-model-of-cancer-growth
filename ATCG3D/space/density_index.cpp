#include "space/density_index.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <vector>

#include "geometry/footprint.hpp"
#include "space/quantized_box_count_index.hpp"

namespace atcg3d {
namespace {

std::int64_t floor_div64(std::int64_t value, std::int64_t divisor) noexcept {
    std::int64_t quotient = value / divisor;
    const std::int64_t remainder = value % divisor;
    if (remainder != 0 && ((remainder < 0) != (divisor < 0))) {
        --quotient;
    }
    return quotient;
}

std::int64_t ceil_div64(std::int64_t value, std::int64_t divisor) noexcept {
    return -floor_div64(-value, divisor);
}

struct ConeVolumeCacheEntry {
    DirectionId direction{};
    int radius{};
    double half_angle_degrees{};
    bool thin_layer{};
    std::size_t site_count{};
};

std::size_t cached_cone_site_count(DirectionId direction,
                                   int radius,
                                   double half_angle_degrees,
                                   bool thin_layer) {
    // Directional density can be evaluated by proposal workers in parallel.
    // A small thread-local numeric cache avoids locks and avoids regenerating
    // the full discrete cone on every cell query.
    thread_local std::vector<ConeVolumeCacheEntry> cache;
    const auto found = std::find_if(
        cache.begin(), cache.end(), [&](const ConeVolumeCacheEntry& entry) {
            return entry.direction == direction && entry.radius == radius &&
                   entry.half_angle_degrees == half_angle_degrees &&
                   entry.thin_layer == thin_layer;
        });
    if (found != cache.end()) return found->site_count;
    const std::vector<Vec3i> offsets =
        directional_cone_offsets(direction, radius, half_angle_degrees);
    const std::size_t count = thin_layer
        ? static_cast<std::size_t>(std::count_if(
              offsets.begin(), offsets.end(),
              [](Vec3i offset) { return offset.z == 0; }))
        : offsets.size();
    if (cache.size() >= 128) cache.clear();
    cache.push_back({direction, radius, half_angle_degrees, thin_layer, count});
    return count;
}

bool inside_directional_cone(Vec3i offset,
                             Vec3i forward,
                             int radius,
                             double minimum_cosine,
                             bool thin_layer) noexcept {
    if (thin_layer && offset.z != 0) return false;
    const int distance = std::max(
        {std::abs(offset.x), std::abs(offset.y), std::abs(offset.z)});
    if (distance == 0 || distance > radius) return false;
    const double offset_length =
        std::sqrt(static_cast<double>(squared_length(offset)));
    const double forward_length =
        std::sqrt(static_cast<double>(squared_length(forward)));
    const double cosine = static_cast<double>(dot(offset, forward)) /
                          (offset_length * forward_length);
    return cosine + 1e-12 >= minimum_cosine;
}

struct InclusiveBox3D {
    Vec3i minimum{};
    Vec3i maximum{};
};

struct BoxDifference3D {
    std::array<InclusiveBox3D, 6> boxes{};
    std::size_t size{};
    std::uint64_t lattice_sites{};
};

bool valid_box(const InclusiveBox3D& box) noexcept {
    return box.minimum.x <= box.maximum.x &&
           box.minimum.y <= box.maximum.y &&
           box.minimum.z <= box.maximum.z;
}

bool contains(const InclusiveBox3D& box, Vec3i site) noexcept {
    return site.x >= box.minimum.x && site.x <= box.maximum.x &&
           site.y >= box.minimum.y && site.y <= box.maximum.y &&
           site.z >= box.minimum.z && site.z <= box.maximum.z;
}

std::uint64_t lattice_site_count(const InclusiveBox3D& box) noexcept {
    if (!valid_box(box)) return 0;
    const auto extent = [](std::int32_t minimum, std::int32_t maximum) {
        return static_cast<std::uint64_t>(
            static_cast<std::int64_t>(maximum) -
            static_cast<std::int64_t>(minimum) + 1);
    };
    return extent(box.minimum.x, box.maximum.x) *
           extent(box.minimum.y, box.maximum.y) *
           extent(box.minimum.z, box.maximum.z);
}

BoxDifference3D subtract_box(const InclusiveBox3D& lhs,
                             const InclusiveBox3D& rhs) {
    BoxDifference3D result;
    const InclusiveBox3D overlap{
        {std::max(lhs.minimum.x, rhs.minimum.x),
         std::max(lhs.minimum.y, rhs.minimum.y),
         std::max(lhs.minimum.z, rhs.minimum.z)},
        {std::min(lhs.maximum.x, rhs.maximum.x),
         std::min(lhs.maximum.y, rhs.maximum.y),
         std::min(lhs.maximum.z, rhs.maximum.z)}};
    const auto append = [&](InclusiveBox3D box) {
        if (!valid_box(box)) return;
        result.boxes[result.size++] = box;
        result.lattice_sites += lattice_site_count(box);
    };
    if (!valid_box(overlap)) {
        append(lhs);
        return result;
    }

    // Six disjoint slabs exactly cover lhs - overlap. Splitting x first, then
    // y inside the overlapping x range, and finally z inside overlapping x/y
    // avoids scanning any lattice coordinate twice.
    if (lhs.minimum.x < overlap.minimum.x) {
        append({lhs.minimum,
                {static_cast<std::int32_t>(overlap.minimum.x - 1),
                 lhs.maximum.y, lhs.maximum.z}});
    }
    if (overlap.maximum.x < lhs.maximum.x) {
        append({{static_cast<std::int32_t>(overlap.maximum.x + 1),
                 lhs.minimum.y, lhs.minimum.z},
                lhs.maximum});
    }
    const std::int32_t middle_min_x = overlap.minimum.x;
    const std::int32_t middle_max_x = overlap.maximum.x;
    if (lhs.minimum.y < overlap.minimum.y) {
        append({{middle_min_x, lhs.minimum.y, lhs.minimum.z},
                {middle_max_x,
                 static_cast<std::int32_t>(overlap.minimum.y - 1),
                 lhs.maximum.z}});
    }
    if (overlap.maximum.y < lhs.maximum.y) {
        append({{middle_min_x,
                 static_cast<std::int32_t>(overlap.maximum.y + 1),
                 lhs.minimum.z},
                {middle_max_x, lhs.maximum.y, lhs.maximum.z}});
    }
    const std::int32_t middle_min_y = overlap.minimum.y;
    const std::int32_t middle_max_y = overlap.maximum.y;
    if (lhs.minimum.z < overlap.minimum.z) {
        append({{middle_min_x, middle_min_y, lhs.minimum.z},
                {middle_max_x, middle_max_y,
                 static_cast<std::int32_t>(overlap.minimum.z - 1)}});
    }
    if (overlap.maximum.z < lhs.maximum.z) {
        append({{middle_min_x, middle_min_y,
                 static_cast<std::int32_t>(overlap.maximum.z + 1)},
                {middle_max_x, middle_max_y, lhs.maximum.z}});
    }
    return result;
}

}  // namespace

BlockDensityIndex3D::BlockDensityIndex3D(int block_edge) : block_edge_(block_edge) {
    if (block_edge_ <= 0) {
        throw std::invalid_argument("density block edge must be positive");
    }
}

void BlockDensityIndex3D::configure_local_window_counts(
    int window_edge, bool thin_layer) {
    if (window_edge <= 0) {
        throw std::invalid_argument("local density window edge must be positive");
    }
    if (!counts_.empty()) {
        throw std::logic_error(
            "local density counts must be configured before inserting anchors");
    }
    local_window_edge_ = window_edge;
    local_window_thin_layer_ = thin_layer;
    local_window_counts_.clear();
}

void BlockDensityIndex3D::begin_local_window_bulk_load() {
    if (local_window_edge_ <= 0) {
        throw std::logic_error("local density counts are not configured");
    }
    if (!counts_.empty()) {
        throw std::logic_error("local density bulk load requires an empty index");
    }
    local_window_counts_.clear();
    local_window_bulk_loading_ = true;
}

void BlockDensityIndex3D::finish_local_window_bulk_load(
    std::size_t slot_count, int workers) {
    if (local_window_edge_ <= 0) return;
    if (!local_window_bulk_loading_) {
        throw std::logic_error("local density bulk load is not active");
    }
    struct Item {
        Slot slot{};
        Vec3i anchor{};
    };
    std::vector<Item> items;
    for (const auto& [block_coordinate_value, block] : counts_) {
        const Vec3i block_min{
            block_coordinate_value.x * block_edge_,
            block_coordinate_value.y * block_edge_,
            block_coordinate_value.z * block_edge_};
        for (const AnchorEntry& entry : block.anchors) {
            if (entry.slot >= slot_count) {
                throw std::logic_error("local density bulk slot exceeds layout");
            }
            items.push_back({
                entry.slot,
                block_min + Vec3i{entry.x, entry.y, entry.z}});
        }
    }
    local_window_counts_.assign(slot_count, {});
    const int lower = (local_window_edge_ - 1) / 2;
    const int upper = local_window_edge_ - lower - 1;
    const auto rebuild = [&](std::size_t index) {
        const Item item = items[index];
        const DensityCounts3D counts = estimate_box(
            item.anchor - Vec3i{lower, lower,
                                local_window_thin_layer_ ? 0 : lower},
            item.anchor + Vec3i{upper, upper,
                                local_window_thin_layer_ ? 0 : upper});
        local_window_counts_[item.slot] = {
            static_cast<std::uint32_t>(counts.r),
            static_cast<std::uint32_t>(counts.K)};
    };
    if (workers <= 1 || items.size() <= 1) {
        for (std::size_t index = 0; index < items.size(); ++index) rebuild(index);
    } else {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(workers)
        for (std::int64_t index = 0;
             index < static_cast<std::int64_t>(items.size()); ++index) {
            rebuild(static_cast<std::size_t>(index));
        }
#else
        for (std::size_t index = 0; index < items.size(); ++index) rebuild(index);
#endif
    }
    local_window_bulk_loading_ = false;
}

DensityCounts3D BlockDensityIndex3D::local_window_counts(Slot slot) const {
    if (local_window_edge_ <= 0 || slot >= local_window_counts_.size()) {
        throw std::out_of_range("local density count slot is unavailable");
    }
    return {local_window_counts_[slot].r, local_window_counts_[slot].K};
}

void BlockDensityIndex3D::set_local_window_counts(
    Slot slot, DensityCounts3D counts) {
    if (counts.r > std::numeric_limits<std::uint32_t>::max() ||
        counts.K > std::numeric_limits<std::uint32_t>::max()) {
        throw std::overflow_error("local density count exceeds uint32 capacity");
    }
    if (slot >= local_window_counts_.size()) {
        local_window_counts_.resize(static_cast<std::size_t>(slot) + 1);
    }
    local_window_counts_[slot] = {
        static_cast<std::uint32_t>(counts.r),
        static_cast<std::uint32_t>(counts.K)};
}

void BlockDensityIndex3D::adjust_local_window_counts(
    Vec3i changed_anchor, CellType changed_type, bool add_value) {
    if (local_window_edge_ <= 0 || local_window_bulk_loading_) return;
    const int lower = (local_window_edge_ - 1) / 2;
    const int upper = local_window_edge_ - lower - 1;
    const Vec3i minimum = changed_anchor -
        Vec3i{upper, upper, local_window_thin_layer_ ? 0 : upper};
    const Vec3i maximum = changed_anchor +
        Vec3i{lower, lower, local_window_thin_layer_ ? 0 : lower};
    for_each_slot_in_box(minimum, maximum, [&](Slot slot) {
        if (slot >= local_window_counts_.size()) {
            throw std::logic_error("local density count slot is missing");
        }
        std::uint32_t& count = changed_type == CellType::r
            ? local_window_counts_[slot].r : local_window_counts_[slot].K;
        if (add_value) {
            if (count == std::numeric_limits<std::uint32_t>::max()) {
                throw std::overflow_error("local density count overflow");
            }
            ++count;
        } else {
            if (count == 0) {
                throw std::logic_error("local density count underflow");
            }
            --count;
        }
    });
}

int BlockDensityIndex3D::floor_div(int value, int divisor) noexcept {
    int quotient = value / divisor;
    const int remainder = value % divisor;
    if (remainder != 0 && ((remainder < 0) != (divisor < 0))) {
        --quotient;
    }
    return quotient;
}

Vec3i BlockDensityIndex3D::block_coordinate(Vec3i site) const noexcept {
    return {floor_div(site.x, block_edge_), floor_div(site.y, block_edge_),
            floor_div(site.z, block_edge_)};
}

BlockDensityIndex3D::AnchorEntry BlockDensityIndex3D::local_entry(
    Vec3i site, CellType type, Slot slot) const noexcept {
    const Vec3i block = block_coordinate(site);
    return {static_cast<std::uint8_t>(site.x - block.x * block_edge_),
            static_cast<std::uint8_t>(site.y - block.y * block_edge_),
            static_cast<std::uint8_t>(site.z - block.z * block_edge_),
            static_cast<std::uint8_t>(type), slot};
}

std::size_t BlockDensityIndex3D::QuantizedCacheKeyHash::operator()(
    const QuantizedCacheKey& key) const noexcept {
    std::uint64_t value = static_cast<std::uint64_t>(Vec3iHash{}(key.coordinate));
    value ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(key.window_edge)) << 17U;
    value ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(key.query_block_edge)) << 43U;
    value ^= key.thin_layer ? 0xd6e8feb86659fd93ULL : 0ULL;
    value ^= value >> 30U;
    value *= 0xbf58476d1ce4e5b9ULL;
    value ^= value >> 27U;
    value *= 0x94d049bb133111ebULL;
    return static_cast<std::size_t>(value ^ (value >> 31U));
}

void BlockDensityIndex3D::register_cache_layout_locked(CacheLayout layout) const {
    if (std::find(cache_layouts_.begin(), cache_layouts_.end(), layout) ==
        cache_layouts_.end()) {
        cache_layouts_.push_back(layout);
        has_quantized_cache_layouts_.store(true, std::memory_order_release);
    }
}

void BlockDensityIndex3D::invalidate_quantized_cache(Vec3i anchor) {
    if (!has_quantized_cache_layouts_.load(std::memory_order_acquire)) return;
    const std::unique_lock lock(quantized_cache_mutex_);
    for (const CacheLayout layout : cache_layouts_) {
        const std::int64_t query_edge = layout.query_block_edge;
        const std::int64_t lower = (layout.window_edge - 1) / 2;
        const std::int64_t upper = layout.window_edge - lower - 1;
        const auto coordinate_range = [&](std::int32_t value) {
            const std::int64_t shifted = static_cast<std::int64_t>(value) - query_edge / 2;
            return std::pair<std::int64_t, std::int64_t>{
                ceil_div64(shifted - upper, query_edge),
                floor_div64(shifted + lower, query_edge)};
        };
        const auto [minimum_x, maximum_x] = coordinate_range(anchor.x);
        const auto [minimum_y, maximum_y] = coordinate_range(anchor.y);
        const auto [minimum_z, maximum_z] = layout.thin_layer
            ? std::pair<std::int64_t, std::int64_t>{
                  anchor.z, anchor.z}
            : coordinate_range(anchor.z);
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
                    quantized_box_cache_.erase({
                        {static_cast<std::int32_t>(x), static_cast<std::int32_t>(y),
                         static_cast<std::int32_t>(z)},
                        layout.window_edge, layout.query_block_edge,
                        layout.thin_layer});
                }
            }
        }
    }
}

void BlockDensityIndex3D::add(Vec3i anchor, CellType type, Slot slot) {
    if (slot == kEmptySlot) {
        throw std::invalid_argument("density index slot must not be empty");
    }
    adjust_local_window_counts(anchor, type, true);
    invalidate_quantized_cache(anchor);
    Block& block = counts_[block_coordinate(anchor)];
    DensityCounts3D& counts = block.counts;
    if (type == CellType::r) {
        ++counts.r;
    } else {
        ++counts.K;
    }
    block.anchors.push_back(local_entry(anchor, type, slot));
    if (incremental_quantized_index_ != nullptr) {
        incremental_quantized_index_->add(anchor);
    }
    if (local_window_edge_ > 0 && !local_window_bulk_loading_) {
        const int lower = (local_window_edge_ - 1) / 2;
        const int upper = local_window_edge_ - lower - 1;
        set_local_window_counts(
            slot,
            estimate_box(
                anchor - Vec3i{lower, lower,
                               local_window_thin_layer_ ? 0 : lower},
                anchor + Vec3i{upper, upper,
                               local_window_thin_layer_ ? 0 : upper}));
    }
}

void BlockDensityIndex3D::remove(Vec3i anchor, CellType type, Slot slot) {
    const Vec3i coordinate = block_coordinate(anchor);
    auto iterator = counts_.find(coordinate);
    if (iterator == counts_.end()) {
        throw std::logic_error("density index removal from absent block");
    }
    const AnchorEntry target = local_entry(anchor, type, slot);
    auto& anchors = iterator->second.anchors;
    const auto found = std::find_if(anchors.begin(), anchors.end(), [&](const AnchorEntry& entry) {
        return entry.x == target.x && entry.y == target.y && entry.z == target.z &&
               entry.type == target.type && entry.slot == target.slot;
    });
    if (found == anchors.end()) {
        throw std::logic_error("density index anchor entry is missing");
    }
    std::uint64_t& value = type == CellType::r
        ? iterator->second.counts.r : iterator->second.counts.K;
    if (value == 0) {
        throw std::logic_error("density index count underflow");
    }
    adjust_local_window_counts(anchor, type, false);
    invalidate_quantized_cache(anchor);
    --value;
    *found = anchors.back();
    anchors.pop_back();
    if (incremental_quantized_index_ != nullptr) {
        incremental_quantized_index_->remove(anchor);
    }
    if (iterator->second.counts.total() == 0) {
        counts_.erase(iterator);
    }
    if (local_window_edge_ > 0 && !local_window_bulk_loading_ &&
        slot < local_window_counts_.size()) {
        local_window_counts_[slot] = {};
    }
}

void BlockDensityIndex3D::move(Vec3i from, Vec3i to, CellType type, Slot slot) {
    last_local_window_move_diagnostics_ = {};
    if (from == to) {
        return;
    }
    if (local_window_edge_ > 0 && !local_window_bulk_loading_) {
        const Vec3i source_coordinate = block_coordinate(from);
        auto source_block = counts_.find(source_coordinate);
        if (source_block == counts_.end()) {
            throw std::logic_error("density index move from absent block");
        }
        const AnchorEntry source = local_entry(from, type, slot);
        auto source_entry = std::find_if(
            source_block->second.anchors.begin(),
            source_block->second.anchors.end(),
            [&](const AnchorEntry& entry) {
                return entry.x == source.x && entry.y == source.y &&
                       entry.z == source.z && entry.type == source.type &&
                       entry.slot == source.slot;
            });
        if (source_entry == source_block->second.anchors.end()) {
            throw std::logic_error("density index move anchor entry is missing");
        }
        if (slot >= local_window_counts_.size()) {
            throw std::logic_error("local density count slot is missing");
        }

        const int lower = (local_window_edge_ - 1) / 2;
        const int upper = local_window_edge_ - lower - 1;
        const int z_lower = local_window_thin_layer_ ? 0 : lower;
        const int z_upper = local_window_thin_layer_ ? 0 : upper;

        // These boxes contain cell anchors whose cached observation windows
        // contain the changed anchor. Anchors in their overlap continue to see
        // this cell after the move and therefore need no counter mutation.
        const InclusiveBox3D old_affected{
            from - Vec3i{upper, upper, z_upper},
            from + Vec3i{lower, lower, z_lower}};
        const InclusiveBox3D new_affected{
            to - Vec3i{upper, upper, z_upper},
            to + Vec3i{lower, lower, z_lower}};
        const BoxDifference3D no_longer_affected =
            subtract_box(old_affected, new_affected);
        const BoxDifference3D newly_affected =
            subtract_box(new_affected, old_affected);
        last_local_window_move_diagnostics_.affected_window_query_sites =
            no_longer_affected.lattice_sites + newly_affected.lattice_sites;

        const auto adjust_cached_slot = [&](Slot affected_slot, bool add_value) {
            ++last_local_window_move_diagnostics_.affected_slot_visits;
            if (affected_slot == slot) return;
            if (affected_slot >= local_window_counts_.size()) {
                throw std::logic_error("local density count slot is missing");
            }
            std::uint32_t& count = type == CellType::r
                ? local_window_counts_[affected_slot].r
                : local_window_counts_[affected_slot].K;
            if (add_value) {
                if (count == std::numeric_limits<std::uint32_t>::max()) {
                    throw std::overflow_error("local density count overflow");
                }
                ++count;
            } else {
                if (count == 0) {
                    throw std::logic_error("local density count underflow");
                }
                --count;
            }
        };
        for (std::size_t index = 0; index < no_longer_affected.size; ++index) {
            const InclusiveBox3D& box = no_longer_affected.boxes[index];
            for_each_slot_in_box(
                box.minimum, box.maximum,
                [&](Slot affected_slot) {
                    adjust_cached_slot(affected_slot, false);
                });
        }
        for (std::size_t index = 0; index < newly_affected.size; ++index) {
            const InclusiveBox3D& box = newly_affected.boxes[index];
            for_each_slot_in_box(
                box.minimum, box.maximum,
                [&](Slot affected_slot) {
                    adjust_cached_slot(affected_slot, true);
                });
        }

        // The moving cell's own observation window is translated too. Update
        // its cached count from the old value using only the two window
        // differences. Its self contribution is present before and after the
        // move, so exclude the old anchor when it lies in a removed slab.
        DensityCounts3D moving_counts{
            local_window_counts_[slot].r,
            local_window_counts_[slot].K};
        const InclusiveBox3D old_window{
            from - Vec3i{lower, lower, z_lower},
            from + Vec3i{upper, upper, z_upper}};
        const InclusiveBox3D new_window{
            to - Vec3i{lower, lower, z_lower},
            to + Vec3i{upper, upper, z_upper}};
        const BoxDifference3D leaving_window =
            subtract_box(old_window, new_window);
        const BoxDifference3D entering_window =
            subtract_box(new_window, old_window);
        last_local_window_move_diagnostics_.moving_window_query_sites =
            leaving_window.lattice_sites + entering_window.lattice_sites;

        const auto counts_without_moving_cell =
            [&](const InclusiveBox3D& box) {
            DensityCounts3D result = estimate_box(box.minimum, box.maximum);
            if (contains(box, from)) {
                std::uint64_t& own_type =
                    type == CellType::r ? result.r : result.K;
                if (own_type == 0) {
                    throw std::logic_error(
                        "moving cell is absent from its density window");
                }
                --own_type;
            }
            return result;
        };
        for (std::size_t index = 0; index < leaving_window.size; ++index) {
            const DensityCounts3D removed =
                counts_without_moving_cell(leaving_window.boxes[index]);
            if (removed.r > moving_counts.r || removed.K > moving_counts.K) {
                throw std::logic_error("moving local density count underflow");
            }
            moving_counts.r -= removed.r;
            moving_counts.K -= removed.K;
        }
        for (std::size_t index = 0; index < entering_window.size; ++index) {
            const DensityCounts3D added =
                counts_without_moving_cell(entering_window.boxes[index]);
            if (added.r > std::numeric_limits<std::uint64_t>::max() -
                              moving_counts.r ||
                added.K > std::numeric_limits<std::uint64_t>::max() -
                              moving_counts.K) {
                throw std::overflow_error("moving local density count overflow");
            }
            moving_counts += added;
        }

        invalidate_quantized_cache(from);
        invalidate_quantized_cache(to);
        const Vec3i target_coordinate = block_coordinate(to);
        const AnchorEntry target = local_entry(to, type, slot);
        if (source_coordinate == target_coordinate) {
            *source_entry = target;
        } else {
            std::uint64_t& source_count = type == CellType::r
                ? source_block->second.counts.r
                : source_block->second.counts.K;
            if (source_count == 0) {
                throw std::logic_error("density index count underflow");
            }
            --source_count;
            *source_entry = source_block->second.anchors.back();
            source_block->second.anchors.pop_back();
            if (source_block->second.counts.total() == 0) {
                counts_.erase(source_block);
            }

            Block& target_block = counts_[target_coordinate];
            if (type == CellType::r) {
                ++target_block.counts.r;
            } else {
                ++target_block.counts.K;
            }
            target_block.anchors.push_back(target);
        }
        if (incremental_quantized_index_ != nullptr) {
            incremental_quantized_index_->move(from, to);
        }
        set_local_window_counts(slot, moving_counts);
        return;
    }
    if (block_coordinate(from) == block_coordinate(to)) {
        invalidate_quantized_cache(from);
        invalidate_quantized_cache(to);
        auto iterator = counts_.find(block_coordinate(from));
        if (iterator == counts_.end()) {
            throw std::logic_error("density index move from absent block");
        }
        const AnchorEntry source = local_entry(from, type, slot);
        const AnchorEntry target = local_entry(to, type, slot);
        const auto found = std::find_if(
            iterator->second.anchors.begin(), iterator->second.anchors.end(),
            [&](const AnchorEntry& entry) {
                return entry.x == source.x && entry.y == source.y && entry.z == source.z &&
                       entry.type == source.type && entry.slot == source.slot;
            });
        if (found == iterator->second.anchors.end()) {
            throw std::logic_error("density index move anchor entry is missing");
        }
        *found = target;
        if (incremental_quantized_index_ != nullptr) {
            incremental_quantized_index_->move(from, to);
        }
        return;
    }
    remove(from, type, slot);
    add(to, type, slot);
}

void BlockDensityIndex3D::reset_quantized_cache() {
    const std::unique_lock lock(quantized_cache_mutex_);
    cache_layouts_.clear();
    quantized_box_cache_.clear();
    has_quantized_cache_layouts_.store(false, std::memory_order_release);
}

DensityCounts3D BlockDensityIndex3D::block_counts(Vec3i block_coordinate_value) const {
    const auto iterator = counts_.find(block_coordinate_value);
    return iterator == counts_.end() ? DensityCounts3D{} : iterator->second.counts;
}

std::size_t BlockDensityIndex3D::allocated_bytes() const {
    std::size_t bytes = counts_.bucket_count() * sizeof(void*);
    for (const auto& [coordinate, block] : counts_) {
        (void)coordinate;
        bytes += sizeof(Vec3i) + sizeof(Block) + 2 * sizeof(void*) +
                 block.anchors.capacity() * sizeof(AnchorEntry);
    }
    bytes += local_window_counts_.capacity() * sizeof(LocalCounts);
    {
        const std::shared_lock lock(quantized_cache_mutex_);
        bytes += quantized_box_cache_.bucket_count() * sizeof(void*) +
                 quantized_box_cache_.size() *
                     (sizeof(QuantizedCacheKey) + sizeof(DensityCounts3D) +
                      2 * sizeof(void*));
        bytes += cache_layouts_.capacity() * sizeof(CacheLayout);
    }
    return bytes;
}

std::size_t BlockDensityIndex3D::quantized_cache_size() const {
    const std::shared_lock lock(quantized_cache_mutex_);
    return quantized_box_cache_.size();
}

DensityCounts3D BlockDensityIndex3D::estimate_box(Vec3i minimum, Vec3i maximum) const {
    const Vec3i first = block_coordinate(minimum);
    const Vec3i last = block_coordinate(maximum);
    DensityCounts3D result;
    for (int bx = first.x; bx <= last.x; ++bx) {
        for (int by = first.y; by <= last.y; ++by) {
            for (int bz = first.z; bz <= last.z; ++bz) {
                const Vec3i coordinate{bx, by, bz};
                const auto iterator = counts_.find(coordinate);
                if (iterator == counts_.end()) continue;
                const Vec3i block_min{bx * block_edge_, by * block_edge_, bz * block_edge_};
                const Vec3i block_max = block_min + Vec3i{block_edge_ - 1,
                                                          block_edge_ - 1,
                                                          block_edge_ - 1};
                const bool complete = minimum.x <= block_min.x && minimum.y <= block_min.y &&
                    minimum.z <= block_min.z && maximum.x >= block_max.x &&
                    maximum.y >= block_max.y && maximum.z >= block_max.z;
                if (complete) {
                    result += iterator->second.counts;
                    continue;
                }
                for (const AnchorEntry& entry : iterator->second.anchors) {
                    const Vec3i anchor = block_min + Vec3i{entry.x, entry.y, entry.z};
                    if (anchor.x < minimum.x || anchor.y < minimum.y || anchor.z < minimum.z ||
                        anchor.x > maximum.x || anchor.y > maximum.y || anchor.z > maximum.z) {
                        continue;
                    }
                    if (entry.type == static_cast<std::uint8_t>(CellType::r)) ++result.r;
                    else ++result.K;
                }
            }
        }
    }
    return result;
}

DensityCounts3D BlockDensityIndex3D::estimate_quantized_box(
    Vec3i anchor, int window_edge, int query_block_edge, bool thin_layer) const {
    if (window_edge <= 0 || query_block_edge <= 0) {
        throw std::invalid_argument("quantized density query edges must be positive");
    }
    const QuantizedCacheKey key{{floor_div(anchor.x, query_block_edge),
                                 floor_div(anchor.y, query_block_edge),
                                 thin_layer ? anchor.z
                                            : floor_div(anchor.z, query_block_edge)},
                                window_edge, query_block_edge, thin_layer};
    {
        const std::shared_lock lock(quantized_cache_mutex_);
        const auto cached = quantized_box_cache_.find(key);
        if (cached != quantized_box_cache_.end()) {
            return cached->second;
        }
    }
    const Vec3i center{key.coordinate.x * query_block_edge + query_block_edge / 2,
                       key.coordinate.y * query_block_edge + query_block_edge / 2,
                       thin_layer
                           ? key.coordinate.z
                           : key.coordinate.z * query_block_edge +
                                 query_block_edge / 2};
    const int lower = (window_edge - 1) / 2;
    const int upper = window_edge - lower - 1;
    const DensityCounts3D result = estimate_box(
        center - Vec3i{lower, lower, thin_layer ? 0 : lower},
        center + Vec3i{upper, upper, thin_layer ? 0 : upper});
    const std::unique_lock lock(quantized_cache_mutex_);
    register_cache_layout_locked({window_edge, query_block_edge, thin_layer});
    const auto cached = quantized_box_cache_.find(key);
    if (cached != quantized_box_cache_.end()) {
        return cached->second;
    }
    if (quantized_box_cache_.size() >= quantized_cache_capacity()) {
        quantized_box_cache_.clear();
    }
    quantized_box_cache_.emplace(key, result);
    return result;
}

double BlockDensityIndex3D::estimate_directional_density(Vec3i anchor,
                                                          DirectionId direction,
                                                          int radius,
                                                          double half_angle_degrees,
                                                          bool thin_layer) const {
    const std::size_t cone_site_count =
        cached_cone_site_count(direction, radius, half_angle_degrees, thin_layer);
    if (cone_site_count == 0) {
        return 1.0;
    }
    const Vec3i forward = direction_vector(direction);
    const double minimum_cosine =
        std::cos(half_angle_degrees * std::acos(-1.0) / 180.0);
    const Vec3i minimum = anchor - Vec3i{radius, radius, thin_layer ? 0 : radius};
    const Vec3i maximum = anchor + Vec3i{radius, radius, thin_layer ? 0 : radius};
    const Vec3i first = block_coordinate(minimum);
    const Vec3i last = block_coordinate(maximum);
    std::uint64_t count = 0;
    for (int bx = first.x; bx <= last.x; ++bx) {
        for (int by = first.y; by <= last.y; ++by) {
            for (int bz = first.z; bz <= last.z; ++bz) {
                const auto iterator = counts_.find({bx, by, bz});
                if (iterator == counts_.end()) continue;
                const Vec3i block_min{
                    bx * block_edge_, by * block_edge_, bz * block_edge_};
                for (const AnchorEntry& entry : iterator->second.anchors) {
                    const Vec3i site =
                        block_min + Vec3i{entry.x, entry.y, entry.z};
                    if (inside_directional_cone(
                            site - anchor, forward, radius, minimum_cosine,
                            thin_layer)) {
                        ++count;
                    }
                }
            }
        }
    }
    return static_cast<double>(count) / static_cast<double>(cone_site_count);
}

std::array<double, 27>
BlockDensityIndex3D::estimate_all_directional_densities(
    Vec3i anchor,
    int radius,
    double half_angle_degrees,
    bool thin_layer) const {
    std::array<double, 27> result{};
    std::array<std::uint64_t, 27> counts{};
    std::array<double, 27> forward_lengths{};
    std::array<std::size_t, 27> cone_site_counts{};
    for (DirectionId direction = 1; direction <= 26; ++direction) {
        forward_lengths[direction] = std::sqrt(
            static_cast<double>(squared_length(direction_vector(direction))));
        cone_site_counts[direction] = cached_cone_site_count(
            direction, radius, half_angle_degrees, thin_layer);
    }

    const double minimum_cosine =
        std::cos(half_angle_degrees * std::acos(-1.0) / 180.0);
    const Vec3i minimum =
        anchor - Vec3i{radius, radius, thin_layer ? 0 : radius};
    const Vec3i maximum =
        anchor + Vec3i{radius, radius, thin_layer ? 0 : radius};
    const Vec3i first = block_coordinate(minimum);
    const Vec3i last = block_coordinate(maximum);
    for (int bx = first.x; bx <= last.x; ++bx) {
        for (int by = first.y; by <= last.y; ++by) {
            for (int bz = first.z; bz <= last.z; ++bz) {
                const auto iterator = counts_.find({bx, by, bz});
                if (iterator == counts_.end()) continue;
                const Vec3i block_min{
                    bx * block_edge_, by * block_edge_, bz * block_edge_};
                for (const AnchorEntry& entry : iterator->second.anchors) {
                    const Vec3i offset =
                        block_min + Vec3i{entry.x, entry.y, entry.z} - anchor;
                    if (thin_layer && offset.z != 0) continue;
                    const int distance = std::max(
                        {std::abs(offset.x), std::abs(offset.y),
                         std::abs(offset.z)});
                    if (distance == 0 || distance > radius) continue;
                    const double offset_length = std::sqrt(
                        static_cast<double>(squared_length(offset)));
                    for (DirectionId direction = 1; direction <= 26;
                         ++direction) {
                        const double cosine =
                            static_cast<double>(
                                dot(offset, direction_vector(direction))) /
                            (offset_length * forward_lengths[direction]);
                        if (cosine + 1e-12 >= minimum_cosine) {
                            ++counts[direction];
                        }
                    }
                }
            }
        }
    }
    for (DirectionId direction = 1; direction <= 26; ++direction) {
        result[direction] = cone_site_counts[direction] == 0
            ? 1.0
            : static_cast<double>(counts[direction]) /
                  static_cast<double>(cone_site_counts[direction]);
    }
    return result;
}

}  // namespace atcg3d
