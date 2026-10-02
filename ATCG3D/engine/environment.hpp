#pragma once

#include <cstddef>
#include <cstdint>

#include "core/types.hpp"

namespace atcg3d {

class CellStore3D;
class SparseVesselGrid3D;

// Supplies a local multiplier for density-dependent biology. A retained
// density of 1 leaves raw counts unchanged; values in (0,1) represent a
// locally increased effective carrying capacity.
class LocalDensityModifier3D {
public:
    virtual ~LocalDensityModifier3D() = default;
    virtual double retained_density(Vec3i site) const noexcept = 0;

    // Normalized effective resource supply used only by explicitly
    // resource-guided migration profiles.  Density-only environments and the
    // uncoupled ABM expose a spatially uniform value, preserving their legacy
    // direction probabilities.
    virtual double normalized_resource(Vec3i site) const noexcept {
        (void)site;
        return 1.0;
    }
    virtual bool contains_resource_site(Vec3i site) const noexcept {
        (void)site;
        return true;
    }
    virtual double shared_density_limit() const noexcept { return 0.0; }
    virtual double shared_carrying_capacity() const noexcept { return 0.0; }
    virtual double growth_resource_scale(Vec3i) const noexcept { return 1.0; }
    virtual bool pure_nutrient_guidance() const noexcept { return false; }
    virtual double nutrient_direction_weight(Vec3i, DirectionId) const { return 1.0; }
};

struct EnvironmentInitializationResult3D {
    bool refresh_cell_rates{};
};

// Optional deterministic environment coupled to the event-driven cell and
// vessel system. The implementation owns its continuous state and refresh
// clock; Simulation3D owns biological rescheduling.
class EnvironmentCoupling3D : public LocalDensityModifier3D {
public:
    ~EnvironmentCoupling3D() override = default;

    virtual EnvironmentInitializationResult3D initialize(
        double now_hours,
        const CellStore3D& cells,
        const SparseVesselGrid3D& vessels) = 0;
    virtual void refresh(double now_hours,
                         const CellStore3D& cells,
                         const SparseVesselGrid3D& vessels) = 0;

    virtual double next_refresh_time_hours() const noexcept = 0;
    virtual std::uint32_t schedule_generation() const noexcept = 0;
    virtual std::uint64_t refresh_count() const noexcept = 0;
    virtual std::size_t allocated_bytes() const noexcept = 0;
    virtual std::uint64_t field_checksum() const noexcept = 0;
    virtual bool individual_refractory() const noexcept { return false; }
    virtual bool activation_ready(CellUid, double, double) { return true; }
    virtual void activation_expired(CellUid, double) {}
};

}  // namespace atcg3d
