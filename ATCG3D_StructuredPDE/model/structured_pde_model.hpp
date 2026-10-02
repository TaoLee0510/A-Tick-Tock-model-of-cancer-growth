#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <memory>
#include <vector>

#include "config/structured_config.hpp"
#include "model/paged_field.hpp"
#include "model/division_renewal.hpp"
#include "model/moving_tumor_front.hpp"

namespace atcg3d {
class Simulation3D;
namespace hybrid { class HybridModel3D; }
}

namespace atcg3d::structured_pde {

enum class StructuredStage3D : std::size_t {
    small = 0,
    large = 1,
};

inline constexpr std::size_t kStructuredStageCount3D = 2;

struct StructuredInitialFields3D {
    std::array<std::vector<double>, kStructuredStageCount3D> r_normal;
    std::array<std::vector<double>, kStructuredStageCount3D> r_active;
    std::array<std::vector<double>, kStructuredStageCount3D> K;
    std::array<std::vector<double>, kStructuredStageCount3D>
        active_remaining_hours;
    std::vector<double> vessel_fraction;
    // Optional v7 ordinary-r subsets. Empty vectors mean no refractory mass.
    std::array<std::vector<double>, kStructuredStageCount3D> r_refractory;
    std::array<std::vector<double>, kStructuredStageCount3D> refractory_remaining_hours;
};

struct StructuredPdeDiagnostics3D {
    std::array<double, kStructuredStageCount3D> r_normal_mass{};
    std::array<double, kStructuredStageCount3D> r_active_mass{};
    std::array<double, kStructuredStageCount3D> K_mass{};
    double r_total{};
    double r_active_total{};
    double active_fraction{};
    double K_total{};
    double assembled_r_consumption_per_hour{};
    double assembled_K_consumption_per_hour{};
    double occupied_volume{};
    double maximum_occupied_fraction{};
    double mean_nutrient{};
    double maximum_nutrient{};
    double r_mean_nutrient{};
    double K_mean_nutrient{};
    double r_mean_radius{};
    double K_mean_radius{};
    double r_radius_50{};
    double r_radius_90{};
    double r_radius_99{};
    double vessel_volume{};
    double tumour_volume{};
    double tumour_front_volume{};
    double tumour_mean_nutrient{};
    double tumour_front_mean_nutrient{};
};

struct StructuredActiveBounds3D {
    int x0{};
    int y0{};
    int z0{};
    int x1{};
    int y1{};
    int z1{};
    bool valid{};
};

struct StructuredDirectionRowSpan3D {
    int dy{};
    int dx0{};
    int dx1{};
};

class StructuredPdeModel3D {
    friend class atcg3d::hybrid::HybridModel3D;
public:
    explicit StructuredPdeModel3D(StructuredPdeConfig3D config);

    void initialize_from_abm(const Simulation3D& simulation);
    void initialize_from_arrays(StructuredInitialFields3D fields,
                                double time_hours = 0.0);
    bool step();

    const StructuredPdeConfig3D& config() const noexcept { return config_; }
    double time_hours() const noexcept { return time_hours_; }
    std::uint64_t step_count() const noexcept { return step_count_; }
    std::uint64_t nutrient_solve_count() const noexcept {
        return nutrient_solve_count_;
    }
    const DivisionRenewal3D* division_renewal() const noexcept { return renewal_.get(); }
    const DivisionRenewal3D* activation_distribution() const noexcept { return duration_.get(); }
    const DivisionRenewal3D* activation_rates() const noexcept { return velocity_.get(); }
    std::size_t voxel_count() const noexcept { return voxel_count_; }
    std::size_t moving_direction_count() const noexcept {
        return direction_ids_.size();
    }
    double voxel_measure() const noexcept { return voxel_measure_; }
    double large_cell_volume() const noexcept { return large_cell_volume_; }
    const continuum::AngiogenesisField3D* angiogenesis() const noexcept { return angiogenesis_.get(); }

    double r_normal(StructuredStage3D stage, std::size_t location) const noexcept;
    double r_active(StructuredStage3D stage, std::size_t location) const noexcept;
    double K(StructuredStage3D stage, std::size_t location) const noexcept;
    double activation_density(StructuredStage3D stage,
                              std::size_t location) const noexcept;
    double occupied_fraction(std::size_t location) const noexcept;
    std::array<double, 2> growth_counts_at(Vec3i site) const;
    double refractory_mass(StructuredStage3D stage, std::size_t location) const noexcept;
    double refractory_mean_hours(StructuredStage3D stage, std::size_t location) const noexcept;
    std::array<double, 3> coordinate(std::size_t location) const noexcept;
    const std::vector<double>& nutrient() const noexcept { return nutrient_; }
    const std::vector<double>& vessel_fraction() const noexcept { return vessel_; }
    const std::vector<std::uint8_t>& tumour_mask() const noexcept {
        return tumour_mask_;
    }

    StructuredPdeDiagnostics3D diagnostics() const;
    std::uint64_t state_checksum() const;
    void save_checkpoint(const std::filesystem::path& path) const;
    void load_checkpoint(const std::filesystem::path& path);

private:
    std::unique_ptr<DivisionRenewal3D> renewal_;
    std::unique_ptr<DivisionRenewal3D> duration_;
    std::unique_ptr<DivisionRenewal3D> velocity_;
    void finish_duration_transport();
    void finish_velocity_transport();
    double division_channel_mass(std::size_t location, std::size_t channel) const;
    void finish_division_transport();
    std::size_t index(int x, int y, int z) const noexcept;
    bool grid_coordinate(Vec3i site, int& x, int& y, int& z) const noexcept;
    void add_synthetic_vessel();
    void advance_angiogenesis(double dt);
    void clear_cells_from_vessels();
    void solve_nutrient();
    void advance_transient_nutrient(double dt);
    void rebuild_moving_tumour_front();
    bool nutrient_source(int x, int y, int z,
                         std::size_t location) const noexcept;
    void build_activation_density();
    void refresh_activation(double dt);
    void migrate_normal_and_K(double dt);
    void migrate_active(double dt);
    void exchange_active_r_with_K(double dt);
    void expire_active(std::size_t stage, double dt);
    void react(double dt);
    void build_local_counts(std::vector<double>& r_counts,
                            std::vector<double>& K_counts) const;
    std::vector<std::size_t> eligible_initial_directions(
        std::size_t location) const;
    std::vector<double> guided_direction_weights(
        std::size_t location,
        bool allow_crowded_sectors = false) const;
    void build_guidance_prefix(double dt);
    bool vessel_blocks_cells(std::size_t location) const noexcept;
    double capacity_multiplier(double nutrient) const noexcept;
    double mean_growth_rate(CellType type) const noexcept;
    double normal_diffusion(StructuredStage3D stage, CellType type) const noexcept;
    double active_rate(StructuredStage3D stage) const noexcept;
    void include_active_location(std::size_t stage,
                                 int x,
                                 int y,
                                 int z) noexcept;
    void include_population_location(int x, int y, int z) noexcept;
    StructuredActiveBounds3D expanded_bounds(
        const StructuredActiveBounds3D& bounds) const noexcept;
    void clear_active_work(const StructuredActiveBounds3D& bounds);
    void shrink_active_bounds(std::size_t stage);
    void shrink_population_bounds();
    void validate_resources() const;
    void validate_state() const;

    StructuredPdeConfig3D config_;
    std::unique_ptr<continuum::AngiogenesisField3D> angiogenesis_;
    // Derived external-agent fields are populated only by hybrid v1. They
    // participate in local rates/resources, never in PDE transported mass.
    std::vector<double> external_r_,external_K_,external_occupied_;
    std::function<void()> external_after_vascular_advance_;
    std::size_t voxel_count_{};
    double voxel_measure_{};
    double large_cell_volume_{};
    std::vector<DirectionId> direction_ids_;
    std::vector<std::vector<std::size_t>> turn_buckets_;
    std::vector<std::vector<StructuredDirectionRowSpan3D>>
        direction_row_spans_;
    std::vector<std::size_t> direction_sector_site_counts_;
    std::vector<double> guidance_density_row_prefix_;
    std::vector<double> guidance_resource_row_prefix_;
    // V4 evaluates guidance from a prefix field that is fixed during one PDE
    // step. Cache each visited location once instead of recomputing eight
    // 70-voxel sectors for every high-rate transport substep.
    std::vector<PagedField<float>> guidance_weight_cache_;
    PagedField<std::uint32_t> guidance_weight_cache_stamp_;
    std::uint32_t guidance_weight_cache_generation_{};
    // V4 transports only locations carrying active-r mass during the many
    // high-rate microsteps. A generation stamp deduplicates target locations
    // without scanning the (potentially very large and mostly empty) bounding
    // rectangle after every microstep.
    PagedField<std::uint32_t> active_location_stamp_;
    std::uint32_t active_location_generation_{};
    StructuredActiveBounds3D guidance_prefix_bounds_;
    int guidance_prefix_pitch_{};
    // Bucket zero is the ABM "no persistent direction yet" state. Buckets
    // 1..N correspond to direction_ids_[bucket-1].
    std::array<std::vector<PagedField<float>>, kStructuredStageCount3D>
        active_direction_;
    std::array<std::vector<PagedField<float>>, kStructuredStageCount3D>
        active_clock_;
    std::array<PagedField<float>, kStructuredStageCount3D> active_total_;
    std::array<StructuredActiveBounds3D, kStructuredStageCount3D>
        active_bounds_;
    StructuredActiveBounds3D work_dirty_bounds_;
    StructuredActiveBounds3D population_bounds_;
    StructuredActiveBounds3D normal_work_dirty_bounds_;
    std::array<PagedField<double>, kStructuredStageCount3D> r_normal_;
    std::array<PagedField<double>, kStructuredStageCount3D> K_;
    std::array<PagedField<double>, kStructuredStageCount3D> r_normal_work_;
    // V7 transports an ordinary-r refractory subset and its time mass with
    // the same flux as normal r. The ratio is a mass-weighted mean clock.
    std::array<PagedField<double>, kStructuredStageCount3D> r_refractory_;
    std::array<PagedField<double>, kStructuredStageCount3D> refractory_clock_;
    std::array<PagedField<double>, kStructuredStageCount3D> refractory_work_;
    std::array<PagedField<double>, kStructuredStageCount3D> refractory_clock_work_;
    std::array<PagedField<double>, kStructuredStageCount3D> K_work_;
    std::vector<PagedField<float>> active_work_;
    std::vector<PagedField<float>> clock_work_;
    std::array<PagedField<double>, kStructuredStageCount3D>
        activation_density_;
    // V5 uses an explicit local refractory/hysteresis closure. An active
    // cohort disarms each location it occupies; the location can only re-arm
    // after its cooldown expires and the 70-window density drops below the
    // configured off threshold. This prevents immediate clock-expiry
    // reactivation in the same crowded neighbourhood.
    std::array<PagedField<float>, kStructuredStageCount3D>
        activation_cooldown_;
    std::array<PagedField<std::uint8_t>, kStructuredStageCount3D>
        activation_armed_;
    StructuredActiveBounds3D activation_density_bounds_;
    std::vector<double> nutrient_;
    std::vector<double> nutrient_next_;
    std::vector<double> vessel_;
    std::vector<std::uint8_t> tumour_mask_;
    std::vector<double> tumour_occupancy_work_;
    std::vector<std::uint8_t> tumour_local_mask_work_;
    continuum::MovingTumorFrontWorkspace2D tumour_front_workspace_;
    StructuredActiveBounds3D tumour_mask_bounds_;
    // Union of the previous and current mask boxes. Source voxels abandoned
    // by a shrinking/re-shaped mask must be refreshed in both nutrient
    // buffers once before they can leave the active solve region.
    StructuredActiveBounds3D nutrient_update_bounds_;
    std::size_t tumour_voxel_count_{};
    std::size_t tumour_front_voxel_count_{};
    double time_hours_{};
    double next_nutrient_refresh_hours_{};
    std::uint64_t step_count_{};
    std::uint64_t nutrient_solve_count_{};
    bool initialized_{};
};

}  // namespace atcg3d::structured_pde
