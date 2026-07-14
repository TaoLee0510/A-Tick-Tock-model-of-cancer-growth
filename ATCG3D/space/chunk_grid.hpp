#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

#include "core/types.hpp"
#include "geometry/footprint.hpp"
#include "space/domain.hpp"

namespace atcg3d {

class SparseChunkGrid3D {
public:
    SparseChunkGrid3D(int chunk_edge, DomainPolicy domain);

    int chunk_edge() const noexcept { return chunk_edge_; }
    std::size_t chunk_count() const noexcept { return chunks_.size(); }
    std::size_t allocated_bytes() const noexcept;

    Slot owner(Vec3i site) const;
    std::vector<Slot> occupants(Vec3i site) const;
    bool empty(Vec3i site) const;
    bool available(Vec3i site) const;

    bool place_single(Vec3i site, Slot slot);
    bool add_colocated(Vec3i site, Slot slot);
    bool remove(Vec3i site, Slot slot);

    bool can_place_large(Vec3i anchor) const;
    bool place_large(Vec3i anchor, Slot slot);
    void remove_large(Vec3i anchor, Slot slot);
    bool can_move_large(Vec3i anchor, Vec3i displacement) const;

private:
    struct Chunk {
        explicit Chunk(std::size_t voxel_count) : owner(voxel_count, kEmptySlot) {}
        std::vector<Slot> owner;
        std::unordered_map<std::uint32_t, std::vector<Slot>> overflow;
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

    int chunk_edge_{};
    std::size_t chunk_voxels_{};
    DomainPolicy domain_;
    std::unordered_map<Vec3i, std::unique_ptr<Chunk>, Vec3iHash> chunks_;
};

}  // namespace atcg3d
