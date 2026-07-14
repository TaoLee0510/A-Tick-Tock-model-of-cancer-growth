#include "space/density_index.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_set>

#include "geometry/footprint.hpp"

namespace atcg3d {

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
    Vec3i site, CellType type) const noexcept {
    const Vec3i block = block_coordinate(site);
    return {static_cast<std::uint8_t>(site.x - block.x * block_edge_),
            static_cast<std::uint8_t>(site.y - block.y * block_edge_),
            static_cast<std::uint8_t>(site.z - block.z * block_edge_),
            static_cast<std::uint8_t>(type)};
}

void BlockDensityIndex3D::add(Vec3i anchor, CellType type) {
    ++generation_;
    Block& block = counts_[block_coordinate(anchor)];
    DensityCounts3D& counts = block.counts;
    if (type == CellType::r) {
        ++counts.r;
    } else {
        ++counts.K;
    }
    block.anchors.push_back(local_entry(anchor, type));
}

void BlockDensityIndex3D::remove(Vec3i anchor, CellType type) {
    ++generation_;
    const Vec3i coordinate = block_coordinate(anchor);
    auto iterator = counts_.find(coordinate);
    if (iterator == counts_.end()) {
        throw std::logic_error("density index removal from absent block");
    }
    std::uint64_t& value = type == CellType::r
        ? iterator->second.counts.r : iterator->second.counts.K;
    if (value == 0) {
        throw std::logic_error("density index count underflow");
    }
    --value;
    const AnchorEntry target = local_entry(anchor, type);
    auto& anchors = iterator->second.anchors;
    const auto found = std::find_if(anchors.begin(), anchors.end(), [&](const AnchorEntry& entry) {
        return entry.x == target.x && entry.y == target.y && entry.z == target.z &&
               entry.type == target.type;
    });
    if (found == anchors.end()) {
        throw std::logic_error("density index anchor entry is missing");
    }
    *found = anchors.back();
    anchors.pop_back();
    if (iterator->second.counts.total() == 0) {
        counts_.erase(iterator);
    }
}

void BlockDensityIndex3D::move(Vec3i from, Vec3i to, CellType type) {
    if (block_coordinate(from) == block_coordinate(to)) {
        ++generation_;
        auto iterator = counts_.find(block_coordinate(from));
        if (iterator == counts_.end()) {
            throw std::logic_error("density index move from absent block");
        }
        const AnchorEntry source = local_entry(from, type);
        const AnchorEntry target = local_entry(to, type);
        const auto found = std::find_if(
            iterator->second.anchors.begin(), iterator->second.anchors.end(),
            [&](const AnchorEntry& entry) {
                return entry.x == source.x && entry.y == source.y && entry.z == source.z &&
                       entry.type == source.type;
            });
        if (found == iterator->second.anchors.end()) {
            throw std::logic_error("density index move anchor entry is missing");
        }
        *found = target;
        return;
    }
    remove(from, type);
    add(to, type);
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
                 (sizeof(Vec3i) + sizeof(CachedBox) + 2 * sizeof(void*));
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
    Vec3i anchor, int window_edge, int query_block_edge) const {
    if (window_edge <= 0 || query_block_edge <= 0) {
        throw std::invalid_argument("quantized density query edges must be positive");
    }
    const Vec3i key{floor_div(anchor.x, query_block_edge),
                    floor_div(anchor.y, query_block_edge),
                    floor_div(anchor.z, query_block_edge)};
    const auto cached = quantized_box_cache_.find(key);
    if (cached != quantized_box_cache_.end() && cached->second.generation == generation_) {
        return cached->second.counts;
    }
    const Vec3i center{key.x * query_block_edge + query_block_edge / 2,
                       key.y * query_block_edge + query_block_edge / 2,
                       key.z * query_block_edge + query_block_edge / 2};
    const int lower = (window_edge - 1) / 2;
    const int upper = window_edge - lower - 1;
    const DensityCounts3D result = estimate_box(
        center - Vec3i{lower, lower, lower}, center + Vec3i{upper, upper, upper});
    quantized_box_cache_[key] = {generation_, result};
    return result;
}

double BlockDensityIndex3D::estimate_directional_density(Vec3i anchor,
                                                          DirectionId direction,
                                                          int radius,
                                                          double half_angle_degrees) const {
    const auto offsets = directional_cone_offsets(direction, radius, half_angle_degrees);
    if (offsets.empty()) {
        return 1.0;
    }
    std::unordered_set<Vec3i, Vec3iHash> sampled_blocks;
    sampled_blocks.reserve(offsets.size());
    for (const Vec3i offset : offsets) {
        sampled_blocks.insert(block_coordinate(anchor + offset));
    }
    DensityCounts3D count;
    for (const Vec3i coordinate : sampled_blocks) {
        count += block_counts(coordinate);
    }
    const double block_volume = static_cast<double>(block_edge_) * block_edge_ * block_edge_;
    return static_cast<double>(count.total()) /
           (static_cast<double>(sampled_blocks.size()) * block_volume);
}

}  // namespace atcg3d
