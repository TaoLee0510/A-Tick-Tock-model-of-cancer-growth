#include "space/density_index.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "geometry/footprint.hpp"

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

}  // namespace

BlockDensityIndex3D::BlockDensityIndex3D(int block_edge) : block_edge_(block_edge) {
    if (block_edge_ <= 0) {
        throw std::invalid_argument("density block edge must be positive");
    }
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

void BlockDensityIndex3D::register_cache_layout(CacheLayout layout) const {
    if (std::find(cache_layouts_.begin(), cache_layouts_.end(), layout) ==
        cache_layouts_.end()) {
        cache_layouts_.push_back(layout);
    }
}

void BlockDensityIndex3D::invalidate_quantized_cache(Vec3i anchor) {
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
    invalidate_quantized_cache(anchor);
    Block& block = counts_[block_coordinate(anchor)];
    DensityCounts3D& counts = block.counts;
    if (type == CellType::r) {
        ++counts.r;
    } else {
        ++counts.K;
    }
    block.anchors.push_back(local_entry(anchor, type, slot));
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
    invalidate_quantized_cache(anchor);
    --value;
    *found = anchors.back();
    anchors.pop_back();
    if (iterator->second.counts.total() == 0) {
        counts_.erase(iterator);
    }
}

void BlockDensityIndex3D::move(Vec3i from, Vec3i to, CellType type, Slot slot) {
    if (from == to) {
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
        return;
    }
    remove(from, type, slot);
    add(to, type, slot);
}

DensityCounts3D BlockDensityIndex3D::block_counts(Vec3i block_coordinate_value) const {
    const auto iterator = counts_.find(block_coordinate_value);
    return iterator == counts_.end() ? DensityCounts3D{} : iterator->second.counts;
}

std::size_t BlockDensityIndex3D::allocated_bytes() const noexcept {
    std::size_t bytes = counts_.bucket_count() * sizeof(void*);
    for (const auto& [coordinate, block] : counts_) {
        (void)coordinate;
        bytes += sizeof(Vec3i) + sizeof(Block) + 2 * sizeof(void*) +
                 block.anchors.capacity() * sizeof(AnchorEntry);
    }
    bytes += quantized_box_cache_.bucket_count() * sizeof(void*) +
             quantized_box_cache_.size() *
                 (sizeof(QuantizedCacheKey) + sizeof(DensityCounts3D) + 2 * sizeof(void*));
    bytes += cache_layouts_.capacity() * sizeof(CacheLayout);
    return bytes;
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
    register_cache_layout({window_edge, query_block_edge, thin_layer});
    const auto cached = quantized_box_cache_.find(key);
    if (cached != quantized_box_cache_.end()) {
        return cached->second;
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

}  // namespace atcg3d
