#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <unordered_map>
#include <vector>

#include "core/types.hpp"
#include "engine/environment.hpp"

namespace atcg3d {

enum class VascularInfluenceProfile3D : std::uint8_t {
    linear_cutoff = 0,
    exponential = 1,
};

// Stores a monotonically increasing relief fraction in sparse dense blocks.
// Overlapping perfused vessel sources combine by max rather than addition.
class VascularInfluenceField3D : public LocalDensityModifier3D {
public:
    VascularInfluenceField3D(int block_edge,
                             float cutoff_radius_voxels,
                             float maximum_relief_fraction,
                             VascularInfluenceProfile3D profile =
                                 VascularInfluenceProfile3D::linear_cutoff,
                             float decay_length_voxels = 1.0F);

    int block_edge() const noexcept { return block_edge_; }
    float cutoff_radius_voxels() const noexcept { return cutoff_radius_voxels_; }
    float maximum_relief_fraction() const noexcept { return maximum_relief_fraction_; }
    VascularInfluenceProfile3D profile() const noexcept { return profile_; }
    float decay_length_voxels() const noexcept { return decay_length_voxels_; }
    std::size_t block_count() const noexcept { return blocks_.size(); }
    std::size_t influenced_voxel_count() const noexcept { return influenced_voxel_count_; }
    std::size_t allocated_bytes() const noexcept;

    float relief(Vec3i site) const;
    double retained_density(Vec3i site) const noexcept override;
    double effective_density(double raw_density, Vec3i site) const;
    double local_capacity_multiplier(Vec3i site) const;

    std::size_t add_source(Vec3i vessel_voxel);
    std::size_t add_sources(std::span<const Vec3i> vessel_voxels);
    std::size_t add_capsule(Vec3i start, Vec3i end, float diameter_voxels);
    void clear() noexcept;

private:
    struct Block {
        explicit Block(std::size_t voxel_count) : relief(voxel_count, 0.0F) {}
        std::vector<float> relief;
    };

    struct Address {
        Vec3i block;
        std::uint32_t local_index{};
    };

    static int floor_div(int value, int divisor) noexcept;
    Address address(Vec3i site) const;
    Block* find_block(Vec3i coordinate);
    const Block* find_block(Vec3i coordinate) const;
    Block& ensure_block(Vec3i coordinate);
    bool raise(Vec3i site, float value);

    int block_edge_{};
    std::size_t block_voxels_{};
    float cutoff_radius_voxels_{};
    float maximum_relief_fraction_{};
    VascularInfluenceProfile3D profile_{VascularInfluenceProfile3D::linear_cutoff};
    float decay_length_voxels_{1.0F};
    std::unordered_map<Vec3i, std::unique_ptr<Block>, Vec3iHash> blocks_;
    std::size_t influenced_voxel_count_{};
};

}  // namespace atcg3d
