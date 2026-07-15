#include "vasculature/vessel_grid.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace atcg3d {

SparseVesselGrid3D::SparseVesselGrid3D(int chunk_edge, DomainPolicy domain)
    : chunk_edge_(chunk_edge),
      chunk_voxels_(static_cast<std::size_t>(chunk_edge) *
                    static_cast<std::size_t>(chunk_edge) *
                    static_cast<std::size_t>(chunk_edge)),
      domain_(std::move(domain)) {
    if (chunk_edge_ <= 0) {
        throw std::invalid_argument("vessel chunk edge must be positive");
    }
}

int SparseVesselGrid3D::floor_div(int value, int divisor) noexcept {
    int quotient = value / divisor;
    const int remainder = value % divisor;
    if (remainder != 0 && ((remainder < 0) != (divisor < 0))) {
        --quotient;
    }
    return quotient;
}

SparseVesselGrid3D::Address SparseVesselGrid3D::address(Vec3i site) const {
    const Vec3i chunk{floor_div(site.x, chunk_edge_), floor_div(site.y, chunk_edge_),
                      floor_div(site.z, chunk_edge_)};
    const int local_x = site.x - chunk.x * chunk_edge_;
    const int local_y = site.y - chunk.y * chunk_edge_;
    const int local_z = site.z - chunk.z * chunk_edge_;
    const auto index = static_cast<std::uint32_t>(
        (local_z * chunk_edge_ + local_y) * chunk_edge_ + local_x);
    return {chunk, index};
}

SparseVesselGrid3D::Chunk* SparseVesselGrid3D::find_chunk(Vec3i coordinate) {
    const auto iterator = chunks_.find(coordinate);
    return iterator == chunks_.end() ? nullptr : iterator->second.get();
}

const SparseVesselGrid3D::Chunk* SparseVesselGrid3D::find_chunk(Vec3i coordinate) const {
    const auto iterator = chunks_.find(coordinate);
    return iterator == chunks_.end() ? nullptr : iterator->second.get();
}

SparseVesselGrid3D::Chunk& SparseVesselGrid3D::ensure_chunk(Vec3i coordinate) {
    auto [iterator, inserted] = chunks_.try_emplace(coordinate);
    if (inserted) {
        iterator->second = std::make_unique<Chunk>(chunk_voxels_);
    }
    return *iterator->second;
}

std::uint8_t& SparseVesselGrid3D::mask_ref(Vec3i site) {
    const Address location = address(site);
    return ensure_chunk(location.chunk).masks[location.local_index];
}

std::uint8_t SparseVesselGrid3D::mask(Vec3i site) const {
    if (!domain_.contains(site)) return kVesselVoxelNone;
    const Address location = address(site);
    const Chunk* chunk = find_chunk(location.chunk);
    return chunk == nullptr ? kVesselVoxelNone : chunk->masks[location.local_index];
}

VesselId SparseVesselGrid3D::vessel_id(Vec3i site) const {
    if (!domain_.contains(site)) return 0;
    const Address location = address(site);
    const Chunk* chunk = find_chunk(location.chunk);
    return chunk == nullptr ? 0 : chunk->vessel_ids[location.local_index];
}

bool SparseVesselGrid3D::all_in_domain(std::span<const Vec3i> sites) const {
    return std::all_of(sites.begin(), sites.end(),
                       [this](Vec3i site) { return domain_.contains(site); });
}

bool SparseVesselGrid3D::any_occupied(std::span<const Vec3i> sites) const {
    return std::any_of(sites.begin(), sites.end(),
                       [this](Vec3i site) { return occupied(site); });
}

VesselPlacementResult3D SparseVesselGrid3D::add(Vec3i site,
                                                 VesselBranchRole role,
                                                 bool is_perfused,
                                                 VesselId owner_vessel_id) {
    return add_sites(std::span<const Vec3i>(&site, 1), role, is_perfused,
                     owner_vessel_id);
}

VesselPlacementResult3D SparseVesselGrid3D::add_sites(std::span<const Vec3i> sites,
                                                       VesselBranchRole role,
                                                       bool is_perfused,
                                                       VesselId owner_vessel_id) {
    if (!all_in_domain(sites)) return {false, 0};
    const std::uint8_t additions = static_cast<std::uint8_t>(
        kVesselVoxelOccupied | vessel_role_bit(role) |
        (is_perfused ? kVesselVoxelPerfused : kVesselVoxelNone));
    std::size_t newly_occupied = 0;
    for (const Vec3i site : sites) {
        const Address location = address(site);
        Chunk& chunk = ensure_chunk(location.chunk);
        std::uint8_t& value = chunk.masks[location.local_index];
        if ((value & kVesselVoxelOccupied) == 0) {
            ++newly_occupied;
            ++occupied_voxel_count_;
            chunk.vessel_ids[location.local_index] = owner_vessel_id;
        } else if (chunk.vessel_ids[location.local_index] == 0 &&
                   owner_vessel_id != 0) {
            // Upgrade legacy/test occupancy without replacing a known owner.
            chunk.vessel_ids[location.local_index] = owner_vessel_id;
        }
        value = static_cast<std::uint8_t>(value | additions);
    }
    return {true, newly_occupied};
}

bool SparseVesselGrid3D::mark_perfused(Vec3i site) {
    if (!occupied(site)) return false;
    std::uint8_t& value = mask_ref(site);
    const bool changed = (value & kVesselVoxelPerfused) == 0;
    value = static_cast<std::uint8_t>(value | kVesselVoxelPerfused);
    return changed;
}

std::size_t SparseVesselGrid3D::mark_perfused(std::span<const Vec3i> sites) {
    std::size_t changed = 0;
    for (const Vec3i site : sites) {
        if (mark_perfused(site)) ++changed;
    }
    return changed;
}

std::vector<Vec3i> SparseVesselGrid3D::occupied_sites() const {
    std::vector<Vec3i> result;
    result.reserve(occupied_voxel_count_);
    for (const auto& [coordinate, chunk] : chunks_) {
        for (std::uint32_t index = 0; index < chunk->masks.size(); ++index) {
            if ((chunk->masks[index] & kVesselVoxelOccupied) == 0) continue;
            const int local_x = static_cast<int>(index % static_cast<std::uint32_t>(chunk_edge_));
            const auto yz = index / static_cast<std::uint32_t>(chunk_edge_);
            const int local_y = static_cast<int>(yz % static_cast<std::uint32_t>(chunk_edge_));
            const int local_z = static_cast<int>(yz / static_cast<std::uint32_t>(chunk_edge_));
            result.push_back({coordinate.x * chunk_edge_ + local_x,
                              coordinate.y * chunk_edge_ + local_y,
                              coordinate.z * chunk_edge_ + local_z});
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

void SparseVesselGrid3D::clear() noexcept {
    chunks_.clear();
    occupied_voxel_count_ = 0;
}

std::size_t SparseVesselGrid3D::allocated_bytes() const noexcept {
    std::size_t bytes = chunks_.bucket_count() * sizeof(void*);
    for (const auto& [coordinate, chunk] : chunks_) {
        (void)coordinate;
        bytes += sizeof(Chunk) + chunk->masks.capacity() * sizeof(std::uint8_t) +
                 chunk->vessel_ids.capacity() * sizeof(VesselId);
    }
    return bytes;
}

}  // namespace atcg3d
