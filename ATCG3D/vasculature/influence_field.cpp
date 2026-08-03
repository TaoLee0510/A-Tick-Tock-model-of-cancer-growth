#include "vasculature/influence_field.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "vasculature/geometry.hpp"

namespace atcg3d {

VascularInfluenceField3D::VascularInfluenceField3D(
    int block_edge,
    float cutoff_radius_voxels,
    float maximum_relief_fraction,
    VascularInfluenceProfile3D profile,
    float decay_length_voxels)
    : block_edge_(block_edge),
      block_voxels_(static_cast<std::size_t>(block_edge) *
                    static_cast<std::size_t>(block_edge) *
                    static_cast<std::size_t>(block_edge)),
      cutoff_radius_voxels_(cutoff_radius_voxels),
      maximum_relief_fraction_(maximum_relief_fraction),
      profile_(profile),
      decay_length_voxels_(decay_length_voxels) {
    if (block_edge_ <= 0) {
        throw std::invalid_argument("vascular influence block edge must be positive");
    }
    if (!(cutoff_radius_voxels_ > 0.0F) || !std::isfinite(cutoff_radius_voxels_)) {
        throw std::invalid_argument("vascular influence cutoff must be finite and positive");
    }
    if (maximum_relief_fraction_ < 0.0F || maximum_relief_fraction_ >= 1.0F ||
        !std::isfinite(maximum_relief_fraction_)) {
        throw std::invalid_argument("vascular maximum relief must be finite in [0,1)");
    }
    if (!(decay_length_voxels_ > 0.0F) || !std::isfinite(decay_length_voxels_)) {
        throw std::invalid_argument("vascular influence decay length must be finite and positive");
    }
}

int VascularInfluenceField3D::floor_div(int value, int divisor) noexcept {
    int quotient = value / divisor;
    const int remainder = value % divisor;
    if (remainder != 0 && ((remainder < 0) != (divisor < 0))) {
        --quotient;
    }
    return quotient;
}

VascularInfluenceField3D::Address VascularInfluenceField3D::address(Vec3i site) const {
    const Vec3i block{floor_div(site.x, block_edge_), floor_div(site.y, block_edge_),
                      floor_div(site.z, block_edge_)};
    const int local_x = site.x - block.x * block_edge_;
    const int local_y = site.y - block.y * block_edge_;
    const int local_z = site.z - block.z * block_edge_;
    const auto index = static_cast<std::uint32_t>(
        (local_z * block_edge_ + local_y) * block_edge_ + local_x);
    return {block, index};
}

VascularInfluenceField3D::Block* VascularInfluenceField3D::find_block(Vec3i coordinate) {
    const auto iterator = blocks_.find(coordinate);
    return iterator == blocks_.end() ? nullptr : iterator->second.get();
}

const VascularInfluenceField3D::Block* VascularInfluenceField3D::find_block(
    Vec3i coordinate) const {
    const auto iterator = blocks_.find(coordinate);
    return iterator == blocks_.end() ? nullptr : iterator->second.get();
}

VascularInfluenceField3D::Block& VascularInfluenceField3D::ensure_block(Vec3i coordinate) {
    auto [iterator, inserted] = blocks_.try_emplace(coordinate);
    if (inserted) iterator->second = std::make_unique<Block>(block_voxels_);
    return *iterator->second;
}

float VascularInfluenceField3D::relief(Vec3i site) const {
    const Address location = address(site);
    const Block* block = find_block(location.block);
    return block == nullptr ? 0.0F : block->relief[location.local_index];
}

double VascularInfluenceField3D::retained_density(Vec3i site) const noexcept {
    return 1.0 - static_cast<double>(relief(site));
}

double VascularInfluenceField3D::effective_density(double raw_density, Vec3i site) const {
    if (!std::isfinite(raw_density) || raw_density < 0.0) {
        throw std::invalid_argument("raw density must be finite and nonnegative");
    }
    return raw_density * retained_density(site);
}

double VascularInfluenceField3D::local_capacity_multiplier(Vec3i site) const {
    return 1.0 / (1.0 - static_cast<double>(relief(site)));
}

bool VascularInfluenceField3D::raise(Vec3i site, float value) {
    if (!(value > 0.0F)) return false;
    const Address location = address(site);
    float& current = ensure_block(location.block).relief[location.local_index];
    if (value <= current) return false;
    if (current == 0.0F) ++influenced_voxel_count_;
    current = value;
    return true;
}

std::size_t VascularInfluenceField3D::add_source(Vec3i vessel_voxel) {
    if (maximum_relief_fraction_ == 0.0F) return 0;
    const int radius = static_cast<int>(std::ceil(cutoff_radius_voxels_));
    const double cutoff = cutoff_radius_voxels_;
    std::size_t changed = 0;
    for (int dx = -radius; dx <= radius; ++dx) {
        for (int dy = -radius; dy <= radius; ++dy) {
            for (int dz = -radius; dz <= radius; ++dz) {
                const double dx_value = dx;
                const double dy_value = dy;
                const double dz_value = dz;
                const double distance = std::sqrt(dx_value * dx_value + dy_value * dy_value +
                                                  dz_value * dz_value);
                if (distance >= cutoff) continue;
                const float value = profile_ == VascularInfluenceProfile3D::linear_cutoff
                    ? static_cast<float>(maximum_relief_fraction_ *
                                         (1.0 - distance / cutoff))
                    : static_cast<float>(maximum_relief_fraction_ *
                                         std::exp(-distance / decay_length_voxels_));
                if (raise(vessel_voxel + Vec3i{dx, dy, dz}, value)) ++changed;
            }
        }
    }
    return changed;
}

std::size_t VascularInfluenceField3D::add_sources(
    std::span<const Vec3i> vessel_voxels) {
    std::size_t changed = 0;
    for (const Vec3i site : vessel_voxels) changed += add_source(site);
    return changed;
}

std::size_t VascularInfluenceField3D::add_capsule(Vec3i start,
                                                   Vec3i end,
                                                   float diameter_voxels) {
    const std::vector<Vec3i> voxels = rasterize_capsule(start, end, diameter_voxels);
    return add_sources(voxels);
}

void VascularInfluenceField3D::clear() noexcept {
    blocks_.clear();
    influenced_voxel_count_ = 0;
}

std::size_t VascularInfluenceField3D::allocated_bytes() const noexcept {
    std::size_t bytes = blocks_.bucket_count() * sizeof(void*);
    for (const auto& [coordinate, block] : blocks_) {
        (void)coordinate;
        bytes += sizeof(Block) + block->relief.capacity() * sizeof(float);
    }
    return bytes;
}

}  // namespace atcg3d
