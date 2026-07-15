#include "space/chunk_grid.hpp"

#include <algorithm>
#include <stdexcept>

#include "vasculature/vessel_grid.hpp"

namespace atcg3d {

SparseChunkGrid3D::SparseChunkGrid3D(int chunk_edge, DomainPolicy domain)
    : chunk_edge_(chunk_edge),
      chunk_voxels_(static_cast<std::size_t>(chunk_edge) * static_cast<std::size_t>(chunk_edge) *
                    static_cast<std::size_t>(chunk_edge)),
      domain_(domain) {
    if (chunk_edge_ <= 0) {
        throw std::invalid_argument("chunk edge must be positive");
    }
}

int SparseChunkGrid3D::floor_div(int value, int divisor) noexcept {
    int quotient = value / divisor;
    const int remainder = value % divisor;
    if (remainder != 0 && ((remainder < 0) != (divisor < 0))) {
        --quotient;
    }
    return quotient;
}

SparseChunkGrid3D::Address SparseChunkGrid3D::address(Vec3i site) const {
    const Vec3i chunk{floor_div(site.x, chunk_edge_), floor_div(site.y, chunk_edge_),
                      floor_div(site.z, chunk_edge_)};
    const int local_x = site.x - chunk.x * chunk_edge_;
    const int local_y = site.y - chunk.y * chunk_edge_;
    const int local_z = site.z - chunk.z * chunk_edge_;
    const auto index = static_cast<std::uint32_t>((local_z * chunk_edge_ + local_y) * chunk_edge_ + local_x);
    return {chunk, index};
}

SparseChunkGrid3D::Chunk* SparseChunkGrid3D::find_chunk(Vec3i coordinate) {
    const auto iterator = chunks_.find(coordinate);
    return iterator == chunks_.end() ? nullptr : iterator->second.get();
}

const SparseChunkGrid3D::Chunk* SparseChunkGrid3D::find_chunk(Vec3i coordinate) const {
    const auto iterator = chunks_.find(coordinate);
    return iterator == chunks_.end() ? nullptr : iterator->second.get();
}

SparseChunkGrid3D::Chunk& SparseChunkGrid3D::ensure_chunk(Vec3i coordinate) {
    auto [iterator, inserted] = chunks_.try_emplace(coordinate);
    if (inserted) {
        iterator->second = std::make_unique<Chunk>(chunk_voxels_);
    }
    return *iterator->second;
}

Slot SparseChunkGrid3D::owner(Vec3i site) const {
    if (!domain_.contains(site)) {
        return kEmptySlot;
    }
    const Address location = address(site);
    const Chunk* chunk = find_chunk(location.chunk);
    return chunk == nullptr ? kEmptySlot : chunk->owner[location.local_index];
}

std::vector<Slot> SparseChunkGrid3D::occupants(Vec3i site) const {
    std::vector<Slot> result;
    if (!domain_.contains(site)) {
        return result;
    }
    const Address location = address(site);
    const Chunk* chunk = find_chunk(location.chunk);
    if (chunk == nullptr || chunk->owner[location.local_index] == kEmptySlot) {
        return result;
    }
    result.push_back(chunk->owner[location.local_index]);
    const auto overflow = chunk->overflow.find(location.local_index);
    if (overflow != chunk->overflow.end()) {
        result.insert(result.end(), overflow->second.begin(), overflow->second.end());
    }
    return result;
}

bool SparseChunkGrid3D::empty(Vec3i site) const {
    return domain_.contains(site) && owner(site) == kEmptySlot;
}

bool SparseChunkGrid3D::available(Vec3i site) const {
    return domain_.contains(site) && empty(site) && !blocked_by_vessel(site);
}

bool SparseChunkGrid3D::blocked_by_vessel(Vec3i site) const {
    return vessels_ != nullptr && vessels_->occupied(site);
}

bool SparseChunkGrid3D::place_single(Vec3i site, Slot slot) {
    if (!available(site) || slot == kEmptySlot) {
        return false;
    }
    const Address location = address(site);
    ensure_chunk(location.chunk).owner[location.local_index] = slot;
    return true;
}

bool SparseChunkGrid3D::add_colocated(Vec3i site, Slot slot) {
    if (!domain_.contains(site) || blocked_by_vessel(site) || slot == kEmptySlot) {
        return false;
    }
    const Address location = address(site);
    Chunk& chunk = ensure_chunk(location.chunk);
    Slot& primary = chunk.owner[location.local_index];
    if (primary == kEmptySlot) {
        primary = slot;
        return true;
    }
    if (primary == slot) {
        return false;
    }
    std::vector<Slot>& overflow = chunk.overflow[location.local_index];
    if (std::find(overflow.begin(), overflow.end(), slot) != overflow.end()) {
        return false;
    }
    overflow.push_back(slot);
    return true;
}

bool SparseChunkGrid3D::remove(Vec3i site, Slot slot) {
    if (!domain_.contains(site)) {
        return false;
    }
    const Address location = address(site);
    Chunk* chunk = find_chunk(location.chunk);
    if (chunk == nullptr) {
        return false;
    }
    Slot& primary = chunk->owner[location.local_index];
    auto overflow_iterator = chunk->overflow.find(location.local_index);
    if (primary == slot) {
        if (overflow_iterator != chunk->overflow.end() && !overflow_iterator->second.empty()) {
            primary = overflow_iterator->second.back();
            overflow_iterator->second.pop_back();
            if (overflow_iterator->second.empty()) {
                chunk->overflow.erase(overflow_iterator);
            }
        } else {
            primary = kEmptySlot;
        }
        return true;
    }
    if (overflow_iterator == chunk->overflow.end()) {
        return false;
    }
    auto& overflow = overflow_iterator->second;
    const auto found = std::find(overflow.begin(), overflow.end(), slot);
    if (found == overflow.end()) {
        return false;
    }
    *found = overflow.back();
    overflow.pop_back();
    if (overflow.empty()) {
        chunk->overflow.erase(overflow_iterator);
    }
    return true;
}

bool SparseChunkGrid3D::can_place_large(Vec3i anchor) const {
    const auto sites = large_footprint(anchor);
    return std::all_of(sites.begin(), sites.end(), [this](Vec3i site) { return available(site); });
}

bool SparseChunkGrid3D::place_large(Vec3i anchor, Slot slot) {
    if (!can_place_large(anchor)) {
        return false;
    }
    for (const Vec3i site : large_footprint(anchor)) {
        const Address location = address(site);
        ensure_chunk(location.chunk).owner[location.local_index] = slot;
    }
    return true;
}

void SparseChunkGrid3D::remove_large(Vec3i anchor, Slot slot) {
    for (const Vec3i site : large_footprint(anchor)) {
        remove(site, slot);
    }
}

bool SparseChunkGrid3D::can_move_large(Vec3i anchor, Vec3i displacement) const {
    const auto entering = entering_voxels(anchor, displacement);
    return std::all_of(entering.begin(), entering.end(), [this](Vec3i site) { return available(site); });
}

std::size_t SparseChunkGrid3D::allocated_bytes() const noexcept {
    std::size_t bytes = chunks_.bucket_count() * sizeof(void*);
    for (const auto& [coordinate, chunk] : chunks_) {
        (void)coordinate;
        bytes += sizeof(Chunk) + chunk->owner.capacity() * sizeof(Slot);
        for (const auto& [index, overflow] : chunk->overflow) {
            (void)index;
            bytes += sizeof(std::uint32_t) + overflow.capacity() * sizeof(Slot);
        }
    }
    return bytes;
}

}  // namespace atcg3d
