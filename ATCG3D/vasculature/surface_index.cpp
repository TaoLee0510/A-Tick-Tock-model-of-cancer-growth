#include "vasculature/surface_index.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

#include "core/stateless_rng.hpp"

namespace atcg3d {
namespace {

constexpr std::uint64_t kSurfaceSamplingEvent = 0x7661736373757266ULL;

struct SurfaceColumnKey3D {
    std::uint8_t axis{};
    std::int32_t transverse_a{};
    std::int32_t transverse_b{};

    bool operator==(const SurfaceColumnKey3D&) const = default;
};

struct SurfaceColumnKey3DHash {
    std::size_t operator()(const SurfaceColumnKey3D& key) const noexcept {
        std::uint64_t value = splitmix64(key.axis);
        value ^= splitmix64(static_cast<std::uint32_t>(key.transverse_a) +
                            0x9e3779b9ULL);
        value ^= splitmix64(static_cast<std::uint32_t>(key.transverse_b) +
                            0x85ebca6bULL);
        return static_cast<std::size_t>(splitmix64(value));
    }
};

struct SurfaceColumnCoordinate3D {
    SurfaceColumnKey3D key;
    std::int32_t axial{};
    std::int8_t sign{};
};

struct SurfaceColumnExtrema3D {
    std::int32_t least_negative{std::numeric_limits<std::int32_t>::max()};
    std::int32_t greatest_positive{std::numeric_limits<std::int32_t>::min()};
    bool has_negative{false};
    bool has_positive{false};
};

SurfaceColumnCoordinate3D surface_column(const ExposedFace3D& face) noexcept {
    if (face.outward_normal.x != 0) {
        return {{0U, face.inside.y, face.inside.z}, face.inside.x,
                static_cast<std::int8_t>(face.outward_normal.x)};
    }
    if (face.outward_normal.y != 0) {
        return {{1U, face.inside.x, face.inside.z}, face.inside.y,
                static_cast<std::int8_t>(face.outward_normal.y)};
    }
    return {{2U, face.inside.x, face.inside.y}, face.inside.z,
            static_cast<std::int8_t>(face.outward_normal.z)};
}

std::uint64_t face_key(const ExposedFace3D& face) noexcept {
    std::uint64_t value = splitmix64(static_cast<std::uint32_t>(face.inside.x));
    value ^= splitmix64(static_cast<std::uint32_t>(face.inside.y) + 0x9e3779b9ULL);
    value ^= splitmix64(static_cast<std::uint32_t>(face.inside.z) + 0x85ebca6bULL);
    value ^= splitmix64(static_cast<std::uint32_t>(face.outward_normal.x + 1) |
                        (static_cast<std::uint64_t>(face.outward_normal.y + 1) << 2U) |
                        (static_cast<std::uint64_t>(face.outward_normal.z + 1) << 4U));
    return splitmix64(value);
}

bool separated_from(const ExposedFace3D& candidate,
                    const std::vector<ExposedFace3D>& selected,
                    double minimum_separation_voxels) {
    if (minimum_separation_voxels <= 0.0) return true;
    const Vec3i candidate_normal = candidate.outward_normal;
    const long double candidate_x = 2.0L * candidate.inside.x + candidate_normal.x;
    const long double candidate_y = 2.0L * candidate.inside.y + candidate_normal.y;
    const long double candidate_z = 2.0L * candidate.inside.z + candidate_normal.z;
    const long double minimum_squared =
        4.0L * minimum_separation_voxels * minimum_separation_voxels;
    for (const ExposedFace3D& existing : selected) {
        const long double dx = candidate_x -
            (2.0L * existing.inside.x + existing.outward_normal.x);
        const long double dy = candidate_y -
            (2.0L * existing.inside.y + existing.outward_normal.y);
        const long double dz = candidate_z -
            (2.0L * existing.inside.z + existing.outward_normal.z);
        if (dx * dx + dy * dy + dz * dz < minimum_squared) return false;
    }
    return true;
}

template <class IncludeFace>
std::vector<ExposedFace3D> sample_stable_hash_top_k(
    const std::unordered_set<ExposedFace3D, ExposedFace3DHash>& faces,
    std::size_t count,
    double minimum_separation_voxels,
    std::uint64_t seed,
    std::uint64_t event_sequence,
    IncludeFace&& include_face) {
    struct RankedFace {
        std::uint64_t rank{};
        ExposedFace3D face;
    };

    const auto ranked_less = [](const RankedFace& lhs, const RankedFace& rhs) {
        if (lhs.rank != rhs.rank) return lhs.rank < rhs.rank;
        return lhs.face < rhs.face;
    };
    const std::size_t candidate_count = std::min(count, faces.size());
    if (candidate_count == 0U) return {};

    // front() is the worst retained candidate. Each better face replaces it,
    // retaining the exact global stable-hash top K in one pass.
    std::vector<RankedFace> ranked;
    ranked.reserve(candidate_count);
    for (const ExposedFace3D& face : faces) {
        if (!include_face(face)) continue;
        const std::uint64_t key = face_key(face);
        RankedFace candidate{
            rng_word(seed, key, kSurfaceSamplingEvent, event_sequence), face};
        if (ranked.size() < candidate_count) {
            ranked.push_back(candidate);
            std::push_heap(ranked.begin(), ranked.end(), ranked_less);
        } else if (ranked_less(candidate, ranked.front())) {
            std::pop_heap(ranked.begin(), ranked.end(), ranked_less);
            ranked.back() = candidate;
            std::push_heap(ranked.begin(), ranked.end(), ranked_less);
        }
    }
    std::sort_heap(ranked.begin(), ranked.end(), ranked_less);

    std::vector<ExposedFace3D> result;
    result.reserve(ranked.size());
    for (const RankedFace& candidate : ranked) {
        if (separated_from(candidate.face, result, minimum_separation_voxels)) {
            result.push_back(candidate.face);
            if (result.size() == count) break;
        }
    }
    return result;
}

}  // namespace

void TumorSurfaceIndex3D::insert_face(const ExposedFace3D& face) {
    const auto [unused, inserted] = faces_.insert(face);
    (void)unused;
    if (!inserted) return;
    centroid_sum_x_ += face.inside.x;
    centroid_sum_y_ += face.inside.y;
    centroid_sum_z_ += face.inside.z;
}

void TumorSurfaceIndex3D::erase_face(const ExposedFace3D& face) noexcept {
    const auto found = faces_.find(face);
    if (found == faces_.end()) return;
    centroid_sum_x_ -= found->inside.x;
    centroid_sum_y_ -= found->inside.y;
    centroid_sum_z_ -= found->inside.z;
    faces_.erase(found);
}

std::vector<ExposedFace3D> TumorSurfaceIndex3D::faces() const {
    std::vector<ExposedFace3D> result(faces_.begin(), faces_.end());
    std::sort(result.begin(), result.end());
    return result;
}

Vec3i TumorSurfaceIndex3D::approximate_centroid() const noexcept {
    if (faces_.empty()) return {};
    const long double count = static_cast<long double>(faces_.size());
    return {static_cast<std::int32_t>(std::llround(centroid_sum_x_ / count)),
            static_cast<std::int32_t>(std::llround(centroid_sum_y_ / count)),
            static_cast<std::int32_t>(std::llround(centroid_sum_z_ / count))};
}

std::vector<ExposedFace3D> TumorSurfaceIndex3D::sample_without_replacement(
    std::size_t count,
    double minimum_separation_voxels,
    std::uint64_t seed,
    std::uint64_t event_sequence) const {
    if (minimum_separation_voxels < 0.0 || !std::isfinite(minimum_separation_voxels)) {
        throw std::invalid_argument("surface seed separation must be finite and nonnegative");
    }
    return sample_stable_hash_top_k(
        faces_, count, minimum_separation_voxels, seed, event_sequence,
        [](const ExposedFace3D&) { return true; });
}

std::vector<ExposedFace3D> TumorSurfaceIndex3D::sample_external_without_replacement(
    std::size_t count,
    double minimum_separation_voxels,
    std::uint64_t seed,
    std::uint64_t event_sequence) const {
    return sample_external_subset_without_replacement(
        count, minimum_separation_voxels, seed, event_sequence,
        [](const ExposedFace3D&) { return true; });
}

std::vector<ExposedFace3D>
TumorSurfaceIndex3D::sample_external_subset_without_replacement(
    std::size_t count,
    double minimum_separation_voxels,
    std::uint64_t seed,
    std::uint64_t event_sequence,
    const std::function<bool(const ExposedFace3D&)>& include_face) const {
    if (minimum_separation_voxels < 0.0 || !std::isfinite(minimum_separation_voxels)) {
        throw std::invalid_argument("surface seed separation must be finite and nonnegative");
    }
    if (!include_face) {
        throw std::invalid_argument("surface subset predicate must be callable");
    }
    if (count == 0U || faces_.empty()) return {};

    // First pass: determine the ray-visible extreme face independently for
    // the positive and negative direction of every transverse lattice column.
    std::unordered_map<SurfaceColumnKey3D,
                       SurfaceColumnExtrema3D,
                       SurfaceColumnKey3DHash> extrema;
    extrema.reserve(faces_.size());
    for (const ExposedFace3D& face : faces_) {
        if (!include_face(face)) continue;
        const SurfaceColumnCoordinate3D column = surface_column(face);
        SurfaceColumnExtrema3D& limits = extrema[column.key];
        if (column.sign > 0) {
            limits.greatest_positive =
                std::max(limits.greatest_positive, column.axial);
            limits.has_positive = true;
        } else {
            limits.least_negative = std::min(limits.least_negative, column.axial);
            limits.has_negative = true;
        }
    }

    // Second pass is performed by the stable-hash sampler. Only O(K) ranked
    // candidates are retained; the surface itself is never copied or sorted.
    const auto is_external = [&extrema, &include_face](const ExposedFace3D& face) {
        if (!include_face(face)) return false;
        const SurfaceColumnCoordinate3D column = surface_column(face);
        const auto found = extrema.find(column.key);
        if (found == extrema.end()) return false;
        if (column.sign > 0) {
            return found->second.has_positive &&
                   column.axial == found->second.greatest_positive;
        }
        return found->second.has_negative &&
               column.axial == found->second.least_negative;
    };
    return sample_stable_hash_top_k(
        faces_, count, minimum_separation_voxels, seed, event_sequence, is_external);
}

}  // namespace atcg3d
