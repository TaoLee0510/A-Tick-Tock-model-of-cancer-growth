#pragma once

#include <array>
#include <filesystem>
#include <functional>
#include <map>
#include <memory>
#include <vector>

#include "engine/environment.hpp"
#include "model/moving_tumor_front.hpp"
#include "config/structured_config.hpp"

namespace atcg3d::hybrid { class HybridModel3D; }

namespace atcg3d::shared_rules {

Model3DConfig abm_config(const structured_pde::StructuredPdeConfig3D& config);

class SharedResourceEnvironment3D final : public EnvironmentCoupling3D {
    friend class atcg3d::hybrid::HybridModel3D;
public:
    explicit SharedResourceEnvironment3D(structured_pde::StructuredPdeConfig3D config);
    EnvironmentInitializationResult3D initialize(double now, const CellStore3D& cells,
                                                  const SparseVesselGrid3D& vessels) override;
    void refresh(double now, const CellStore3D& cells, const SparseVesselGrid3D& vessels) override;
    double retained_density(Vec3i) const noexcept override { return 1.0; }
    double shared_density_limit() const noexcept override;
    double shared_carrying_capacity() const noexcept override;
    double growth_resource_scale(Vec3i site) const noexcept override;
    double normalized_resource(Vec3i site) const noexcept override;
    bool contains_resource_site(Vec3i site) const noexcept override;
    bool pure_nutrient_guidance() const noexcept override { return true; }
    double nutrient_direction_weight(Vec3i site, DirectionId direction) const override;
    std::array<double,2> external_growth_counts(Vec3i site) const override { return external_counts_ ? external_counts_(site) : std::array<double,2>{}; }
    double external_activation_density(Vec3i site, CellStage stage) const override { return external_activation_ ? external_activation_(site,stage) : 0.0; }
    bool destination_available(Vec3i site) const noexcept override { return !external_destination_ || external_destination_(site); }
    bool individual_refractory() const noexcept override { return true; }
    bool activation_ready(CellUid uid, double now, double density) override;
    void activation_expired(CellUid uid, double now) override;
    bool migration_operator_enabled() const noexcept override {
        return config_.migration_operator_enabled;
    }
    bool division_operator_enabled() const noexcept override {
        return config_.division_operator_enabled;
    }
    double next_refresh_time_hours() const noexcept override { return next_refresh_; }
    std::uint32_t schedule_generation() const noexcept override { return generation_; }
    std::uint64_t refresh_count() const noexcept override { return refreshes_; }
    std::size_t allocated_bytes() const noexcept override;
    std::uint64_t field_checksum() const noexcept override;
    const std::vector<double>& nutrient() const noexcept { return nutrient_; }
    const std::vector<double>& vessel_fraction() const noexcept { return vessels_; }
    const std::vector<std::uint8_t>& tumour_mask() const noexcept { return tumour_mask_; }
    void save_checkpoint(const std::filesystem::path& path, std::uint64_t abm_checksum) const;
    void load_checkpoint(const std::filesystem::path& path, std::uint64_t abm_checksum);

private:
    bool externally_driven_{};
    std::function<std::array<double,2>(Vec3i)> external_counts_;
    std::function<double(Vec3i,CellStage)> external_activation_;
    std::function<bool(Vec3i)> external_destination_;
    struct Refractory { double until{}; bool armed{true}; };
    struct RowSpan { int dy{}, dx0{}, dx1{}; };
    std::unique_ptr<continuum::AngiogenesisField3D> angiogenesis_;
    std::size_t index(int x, int y, int z) const noexcept;
    std::size_t location(Vec3i site) const noexcept;
    void assemble(const CellStore3D& cells, const SparseVesselGrid3D& vessels);
    bool source(int x, int y, int z, std::size_t here) const noexcept;
    void rebuild_prefix();
    void rebuild_sources();
    structured_pde::StructuredPdeConfig3D config_;
    std::uint64_t fingerprint_{};
    StaticVascularGeometry3D geometry_;
    std::vector<double> nutrient_, next_, consumers_, occupied_, vessels_, row_prefix_;
    std::vector<std::uint8_t> tumour_mask_;
    continuum::MovingTumorFrontWorkspace2D front_workspace_;
    std::array<std::vector<Vec3i>, 27> sector_offsets_;
    std::array<std::vector<RowSpan>, 27> sector_rows_;
    std::map<CellUid, Refractory> refractory_;
    double last_refresh_{}, next_refresh_{};
    std::uint32_t generation_{1};
    std::uint64_t refreshes_{};
    bool restored_{};
};

}  // namespace atcg3d::shared_rules
