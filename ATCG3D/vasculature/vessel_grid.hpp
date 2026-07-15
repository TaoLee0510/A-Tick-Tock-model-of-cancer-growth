#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <unordered_map>
#include <vector>

#include "space/domain.hpp"
#include "vasculature/types.hpp"

namespace atcg3d {

struct VesselPlacementResult3D {
    bool placed{};
    std::size_t newly_occupied_voxels{};
};

class SparseVesselGrid3D {
public:
    SparseVesselGrid3D(int chunk_edge, DomainPolicy domain);

    int chunk_edge() const noexcept { return chunk_edge_; }
    std::size_t chunk_count() const noexcept { return chunks_.size(); }
    std::size_t occupied_voxel_count() const noexcept { return occupied_voxel_count_; }
    std::size_t allocated_bytes() const noexcept;

    bool in_domain(Vec3i site) const noexcept { return domain_.contains(site); }
    std::uint8_t mask(Vec3i site) const;
    bool occupied(Vec3i site) const {
        return (mask(site) & kVesselVoxelOccupied) != 0;
    }
    bool perfused(Vec3i site) const {
        return (mask(site) & kVesselVoxelPerfused) != 0;
    }
    VesselId vessel_id(Vec3i site) const;
    bool has_any(Vec3i site, std::uint8_t bits) const {
        return (mask(site) & bits) != 0;
    }

    bool all_in_domain(std::span<const Vec3i> sites) const;
    bool any_occupied(std::span<const Vec3i> sites) const;
    VesselPlacementResult3D add(Vec3i site,
                                VesselBranchRole role,
                                bool perfused = false,
                                VesselId vessel_id = 0);
    VesselPlacementResult3D add_sites(std::span<const Vec3i> sites,
                                      VesselBranchRole role,
                                      bool perfused = false,
                                      VesselId vessel_id = 0);
    bool mark_perfused(Vec3i site);
    std::size_t mark_perfused(std::span<const Vec3i> sites);

    std::vector<Vec3i> occupied_sites() const;
    void clear() noexcept;

private:
    struct Chunk {
        explicit Chunk(std::size_t voxel_count)
            : masks(voxel_count, kVesselVoxelNone), vessel_ids(voxel_count, 0) {}
        std::vector<std::uint8_t> masks;
        // Primary owner is used only for collision/anastomosis decisions. Role
        // bits still retain every overlapping branch classification.
        std::vector<VesselId> vessel_ids;
    };

    struct Address {
        Vec3i chunk;
        std::uint32_t local_index{};
    };

    static int floor_div(int value, int divisor) noexcept;
    Address address(Vec3i site) const;
    Chunk* find_chunk(Vec3i coordinate);
    const Chunk* find_chunk(Vec3i coordinate) const;
    Chunk& ensure_chunk(Vec3i coordinate);
    std::uint8_t& mask_ref(Vec3i site);

    int chunk_edge_{};
    std::size_t chunk_voxels_{};
    DomainPolicy domain_;
    std::unordered_map<Vec3i, std::unique_ptr<Chunk>, Vec3iHash> chunks_;
    std::size_t occupied_voxel_count_{};
};

}  // namespace atcg3d
