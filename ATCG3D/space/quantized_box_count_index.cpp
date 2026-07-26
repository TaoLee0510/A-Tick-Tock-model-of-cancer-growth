#include "space/quantized_box_count_index.hpp"

#include <algorithm>
#include <stdexcept>

namespace atcg3d {

QuantizedBoxCountIndex3D::QuantizedBoxCountIndex3D(
    int window_edge, int query_block_edge, bool thin_layer)
    : window_edge_(window_edge),
      query_block_edge_(query_block_edge),
      thin_layer_(thin_layer) {
    if (window_edge_ <= 0 || query_block_edge_ <= 0) {
        throw std::invalid_argument(
            "quantized box count index edges must be positive");
    }
}

std::int64_t QuantizedBoxCountIndex3D::floor_div(
    std::int64_t value, std::int64_t divisor) noexcept {
    std::int64_t quotient = value / divisor;
    if (value % divisor < 0) --quotient;
    return quotient;
}

std::int64_t QuantizedBoxCountIndex3D::ceil_div(
    std::int64_t value, std::int64_t divisor) noexcept {
    return -floor_div(-value, divisor);
}

Vec3i QuantizedBoxCountIndex3D::query_block(Vec3i anchor) const noexcept {
    return {static_cast<std::int32_t>(floor_div(anchor.x, query_block_edge_)),
            static_cast<std::int32_t>(floor_div(anchor.y, query_block_edge_)),
            thin_layer_
                ? anchor.z
                : static_cast<std::int32_t>(
                      floor_div(anchor.z, query_block_edge_))};
}

std::vector<Vec3i> QuantizedBoxCountIndex3D::affected_blocks(
    Vec3i anchor) const {
    std::vector<Vec3i> result;
    for_each_affected_block(
        anchor, [&](Vec3i block) { result.push_back(block); });
    return result;
}

void QuantizedBoxCountIndex3D::add(Vec3i anchor) {
    for_each_affected_block(anchor, [&](Vec3i block) { ++counts_[block]; });
    ++residents_[query_block(anchor)];
}

void QuantizedBoxCountIndex3D::remove(Vec3i anchor) {
    for_each_affected_block(anchor, [&](Vec3i block) {
        const auto found = counts_.find(block);
        if (found == counts_.end() || found->second == 0) {
            throw std::logic_error("quantized box count index underflow");
        }
        if (--found->second == 0) counts_.erase(found);
    });
    const Vec3i resident_block = query_block(anchor);
    const auto resident = residents_.find(resident_block);
    if (resident == residents_.end() || resident->second == 0) {
        throw std::logic_error("quantized box resident count underflow");
    }
    if (--resident->second == 0) residents_.erase(resident);
}

void QuantizedBoxCountIndex3D::move(Vec3i from, Vec3i to) {
    if (from == to) return;
    remove(from);
    add(to);
}

std::uint64_t QuantizedBoxCountIndex3D::count(Vec3i query_block_value) const
    noexcept {
    const auto found = counts_.find(query_block_value);
    return found == counts_.end() ? 0 : found->second;
}

std::uint64_t QuantizedBoxCountIndex3D::resident_count(
    Vec3i query_block_value) const noexcept {
    const auto found = residents_.find(query_block_value);
    return found == residents_.end() ? 0 : found->second;
}

std::vector<Vec3i> QuantizedBoxCountIndex3D::resident_blocks() const {
    std::vector<Vec3i> result;
    result.reserve(residents_.size());
    for (const auto& [block, count_value] : residents_) {
        if (count_value != 0) result.push_back(block);
    }
    std::sort(result.begin(), result.end());
    return result;
}

std::size_t QuantizedBoxCountIndex3D::allocated_bytes() const noexcept {
    return counts_.bucket_count() * sizeof(void*) +
           counts_.size() *
               (sizeof(Vec3i) + sizeof(std::uint64_t) + 2 * sizeof(void*)) +
           residents_.bucket_count() * sizeof(void*) +
           residents_.size() *
               (sizeof(Vec3i) + sizeof(std::uint64_t) + 2 * sizeof(void*));
}

}  // namespace atcg3d
