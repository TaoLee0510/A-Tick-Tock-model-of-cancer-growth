#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <unordered_set>
#include <vector>

#include "core/types.hpp"

namespace atcg3d {

inline constexpr std::array<Vec3i, 6> kAxisFaceNormals3D{{
    {-1, 0, 0}, {1, 0, 0}, {0, -1, 0},
    {0, 1, 0}, {0, 0, -1}, {0, 0, 1},
}};

struct ExposedFace3D {
    Vec3i inside{};
    Vec3i outward_normal{};

    Vec3i outside() const { return inside + outward_normal; }
    constexpr auto operator<=>(const ExposedFace3D&) const = default;
};

struct ExposedFace3DHash {
    std::size_t operator()(const ExposedFace3D& face) const noexcept {
        const std::size_t inside = Vec3iHash{}(face.inside);
        const std::size_t normal = Vec3iHash{}(face.outward_normal);
        return inside ^ (normal + static_cast<std::size_t>(0x9e3779b9U) +
                         (inside << 6U) + (inside >> 2U));
    }
};

class TumorSurfaceIndex3D {
public:
    std::size_t size() const noexcept { return faces_.size(); }
    bool empty() const noexcept { return faces_.empty(); }
    bool contains(const ExposedFace3D& face) const { return faces_.contains(face); }
    void clear() noexcept {
        faces_.clear();
        centroid_sum_x_ = 0.0L;
        centroid_sum_y_ = 0.0L;
        centroid_sum_z_ = 0.0L;
    }

    template <class IsOccupied>
    void rebuild(std::span<const Vec3i> occupied_sites, IsOccupied&& is_occupied) {
        clear();
        for (const Vec3i inside : occupied_sites) {
            if (!is_occupied(inside)) continue;
            for (const Vec3i normal : kAxisFaceNormals3D) {
                if (!is_occupied(inside + normal)) {
                    insert_face({inside, normal});
                }
            }
        }
    }

    template <class VisitOccupied, class IsOccupied>
    void rebuild_from_visitor(VisitOccupied&& visit_occupied, IsOccupied&& is_occupied) {
        clear();
        visit_occupied([&](Vec3i inside) {
            if (!is_occupied(inside)) return;
            for (const Vec3i normal : kAxisFaceNormals3D) {
                if (!is_occupied(inside + normal)) {
                    insert_face({inside, normal});
                }
            }
        });
    }

    // changed_sites are sites whose occupancy may have changed. Each changed
    // site and its six possible inside neighbors are recomputed locally.
    template <class IsOccupied>
    void refresh(std::span<const Vec3i> changed_sites, IsOccupied&& is_occupied) {
        std::unordered_set<Vec3i, Vec3iHash> candidates;
        candidates.reserve(changed_sites.size() * 7U);
        for (const Vec3i site : changed_sites) {
            candidates.insert(site);
            for (const Vec3i normal : kAxisFaceNormals3D) {
                candidates.insert(site + normal);
            }
        }
        for (const Vec3i inside : candidates) {
            for (const Vec3i normal : kAxisFaceNormals3D) {
                erase_face({inside, normal});
            }
            if (!is_occupied(inside)) continue;
            for (const Vec3i normal : kAxisFaceNormals3D) {
                if (!is_occupied(inside + normal)) {
                    insert_face({inside, normal});
                }
            }
        }
    }

    std::vector<ExposedFace3D> faces() const;
    Vec3i approximate_centroid() const noexcept;
    std::vector<ExposedFace3D> sample_without_replacement(
        std::size_t count,
        double minimum_separation_voxels,
        std::uint64_t seed,
        std::uint64_t event_sequence = 0) const;

    // Samples only faces whose outward-normal ray is unobstructed to infinity.
    // For each transverse lattice column, this is the greatest axial + face
    // and the least axial - face. Closed-cavity walls are therefore excluded.
    std::vector<ExposedFace3D> sample_external_without_replacement(
        std::size_t count,
        double minimum_separation_voxels,
        std::uint64_t seed,
        std::uint64_t event_sequence = 0) const;

private:
    void insert_face(const ExposedFace3D& face);
    void erase_face(const ExposedFace3D& face) noexcept;

    std::unordered_set<ExposedFace3D, ExposedFace3DHash> faces_;
    // These sums make the seed-orientation centroid query O(1). Long double
    // avoids overflowing fixed-width integer accumulators for sparse domains
    // with large signed coordinates.
    long double centroid_sum_x_{0.0L};
    long double centroid_sum_y_{0.0L};
    long double centroid_sum_z_{0.0L};
};

}  // namespace atcg3d
