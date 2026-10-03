#pragma once

#include <array>
#include <bit>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/types.hpp"

namespace atcg3d {

// Shared finite-grid source geometry for the versioned ABM/PDE contract.
// Integer ABM sites are mapped to the same voxel centres used by the PDE.
struct StaticVascularGeometry3D {
    std::string model{"disabled"};
    std::array<int, 3> shape{1, 1, 1};
    std::array<double, 3> origin{0.0, 0.0, -0.5};
    double spacing_voxels{1.0};
    std::string source_mode{"abm_perfusion"};
    std::string synthetic_axis{"y"};
    std::array<double, 3> synthetic_center{};
    double synthetic_radius_voxels{1.5};
    std::vector<Vec3i> static_sources;
    bool thin_layer{true};

    bool enabled() const noexcept { return model != "disabled"; }
    bool contains(Vec3i site) const noexcept {
        if (!enabled()) return true;
        const std::array<double, 3> point{
            static_cast<double>(site.x), static_cast<double>(site.y),
            static_cast<double>(thin_layer ? 0 : site.z)};
        for (int axis = 0; axis < 3; ++axis) {
            const double coordinate = (point[axis] - origin[axis]) / spacing_voxels;
            if (coordinate < 0.0 || coordinate >= shape[axis]) return false;
        }
        return true;
    }
    bool source(Vec3i site) const noexcept {
        if (!enabled() || !contains(site)) return false;
        const std::array<double, 3> point{
            static_cast<double>(site.x), static_cast<double>(site.y),
            static_cast<double>(thin_layer ? 0 : site.z)};
        std::array<int, 3> voxel{};
        for (int axis = 0; axis < 3; ++axis) {
            voxel[axis] = static_cast<int>(std::floor(
                (point[axis] - origin[axis]) / spacing_voxels));
        }
        return source_voxel(voxel);
    }
    bool source_voxel(const std::array<int, 3>& voxel) const noexcept {
        if (!enabled()) return false;
        for (const Vec3i configured : static_sources) {
            const std::array<double, 3> p{static_cast<double>(configured.x),
                static_cast<double>(configured.y), static_cast<double>(configured.z)};
            bool match = true;
            for (int axis = 0; axis < 3; ++axis) {
                match = match && voxel[axis] == static_cast<int>(std::floor(
                    (p[axis] - origin[axis]) / spacing_voxels));
            }
            if (match) return true;
        }
        if (source_mode != "synthetic_central_line" &&
            source_mode != "abm_plus_synthetic_line") return false;
        const int line_axis = synthetic_axis == "x" ? 0 :
            (synthetic_axis == "y" ? 1 : 2);
        double distance_squared = 0.0;
        for (int axis = 0; axis < 3; ++axis) {
            if (axis == line_axis || (thin_layer && axis == 2)) continue;
            const double centre = origin[axis] + (voxel[axis] + 0.5) * spacing_voxels;
            const double offset = centre - synthetic_center[axis];
            distance_squared += offset * offset;
        }
        return distance_squared <= synthetic_radius_voxels * synthetic_radius_voxels;
    }
    void validate() const {
        if (!enabled()) return;
        if (model != "shared_static_v3" || !(spacing_voxels > 0.0) ||
            !std::isfinite(spacing_voxels) || !(synthetic_radius_voxels > 0.0) ||
            !std::isfinite(synthetic_radius_voxels) ||
            (synthetic_axis != "x" && synthetic_axis != "y" && synthetic_axis != "z") ||
            (source_mode != "abm_perfusion" && source_mode != "static_voxels" &&
             source_mode != "synthetic_central_line" &&
             source_mode != "abm_plus_synthetic_line")) {
            throw std::invalid_argument("invalid shared static vascular geometry");
        }
        std::uint64_t voxel_count = 1;
        for (int axis = 0; axis < 3; ++axis) {
            if (shape[axis] <= 0 || !std::isfinite(origin[axis]) ||
                !std::isfinite(synthetic_center[axis]) || shape[axis] > 32768 ||
                voxel_count > 500000000ULL / static_cast<std::uint64_t>(shape[axis])) {
                throw std::invalid_argument("invalid shared vascular grid");
            }
            voxel_count *= static_cast<std::uint64_t>(shape[axis]);
        }
        if (thin_layer != (shape[2] == 1)) {
            throw std::invalid_argument("shared vascular z extent disagrees with dimensionality");
        }
        for (const Vec3i site : static_sources) {
            if (!contains(site)) throw std::invalid_argument("static source outside vascular grid");
        }
    }
    std::string to_json() const {
        std::ostringstream out;
        out << std::setprecision(17) << "{\"model\":\"" << model
            << "\",\"shape\":[" << shape[0] << ',' << shape[1] << ',' << shape[2]
            << "],\"origin\":[" << origin[0] << ',' << origin[1] << ',' << origin[2]
            << "],\"spacing_voxels\":" << spacing_voxels
            << ",\"source_mode\":\"" << source_mode << "\",\"synthetic_axis\":\""
            << synthetic_axis << "\",\"synthetic_center\":[" << synthetic_center[0]
            << ',' << synthetic_center[1] << ',' << synthetic_center[2]
            << "],\"synthetic_radius_voxels\":" << synthetic_radius_voxels
            << ",\"thin_layer\":" << (thin_layer ? "true" : "false")
            << ",\"static_sources\":[";
        for (std::size_t i = 0; i < static_sources.size(); ++i) {
            if (i != 0) out << ',';
            out << '[' << static_sources[i].x << ',' << static_sources[i].y << ','
                << static_sources[i].z << ']';
        }
        out << "]}";
        return out.str();
    }
    std::uint64_t fingerprint() const noexcept {
        std::uint64_t state = 0x5354415449435633ULL;
        const auto add = [&state](std::uint64_t value) {
            value += 0x9e3779b97f4a7c15ULL;
            value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
            value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
            value ^= value >> 31U;
            state ^= value + (state << 6U) + (state >> 2U);
        };
        for (const auto* text : {&model, &source_mode, &synthetic_axis}) {
            add(text->size());
            for (const unsigned char c : *text) add(c);
        }
        for (int axis = 0; axis < 3; ++axis) {
            add(shape[axis]);
            add(std::bit_cast<std::uint64_t>(origin[axis]));
            add(std::bit_cast<std::uint64_t>(synthetic_center[axis]));
        }
        add(std::bit_cast<std::uint64_t>(spacing_voxels));
        add(std::bit_cast<std::uint64_t>(synthetic_radius_voxels));
        add(thin_layer);
        add(static_sources.size());
        for (const Vec3i site : static_sources) { add(site.x); add(site.y); add(site.z); }
        return state;
    }
};

}  // namespace atcg3d
