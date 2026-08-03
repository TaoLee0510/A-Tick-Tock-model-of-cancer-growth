#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <vector>

#include "config/continuum_config.hpp"

namespace atcg3d {
class Simulation3D;
}

namespace atcg3d::continuum {

enum class PopulationField3D : std::size_t {
    r_small = 0,
    r_large = 1,
    K_small = 2,
    K_large = 3,
};

inline constexpr std::size_t kPopulationFieldCount3D = 4;

struct ContinuumDiagnostics3D {
    std::array<double, kPopulationFieldCount3D> population_mass{};
    std::array<double, 2> type_mass{};
    std::array<double, 2> mean_nutrient_by_type{};
    std::array<double, 2> mean_radius_by_type{};
    double occupied_volume{};
    double maximum_occupied_fraction{};
    double mean_nutrient{};
    double maximum_nutrient{};
    double vessel_volume{};
};

class ContinuumModel3D {
public:
    explicit ContinuumModel3D(ContinuumModelConfig3D config);

    void initialize_from_abm(const Simulation3D& simulation);
    void initialize_from_arrays(
        std::array<std::vector<double>, kPopulationFieldCount3D> populations,
        std::vector<double> vessel_fraction,
        double time_hours = 0.0);
    bool step();

    const ContinuumModelConfig3D& config() const noexcept { return config_; }
    double time_hours() const noexcept { return time_hours_; }
    std::uint64_t step_count() const noexcept { return step_count_; }
    std::uint64_t nutrient_solve_count() const noexcept {
        return nutrient_solve_count_;
    }
    std::size_t voxel_count() const noexcept { return voxel_count_; }
    double voxel_measure() const noexcept { return voxel_measure_; }
    double large_cell_volume() const noexcept { return large_cell_volume_; }
    const std::vector<double>& population(PopulationField3D field) const noexcept {
        return populations_[static_cast<std::size_t>(field)];
    }
    const std::vector<double>& nutrient() const noexcept { return nutrient_; }
    const std::vector<double>& vessel_fraction() const noexcept { return vessel_; }
    double occupied_fraction(std::size_t index) const noexcept;
    std::array<double, 3> coordinate(std::size_t index) const noexcept;
    ContinuumDiagnostics3D diagnostics() const;
    std::uint64_t state_checksum() const;

    void save_checkpoint(const std::filesystem::path& path) const;
    void load_checkpoint(const std::filesystem::path& path);

private:
    std::size_t index(int x, int y, int z) const noexcept;
    bool grid_coordinate(Vec3i site, int& x, int& y, int& z) const noexcept;
    void add_synthetic_vessel();
    void solve_nutrient();
    void migrate(double dt);
    void react(double dt);
    void build_local_counts(std::vector<double>& r_counts,
                            std::vector<double>& K_counts) const;
    double capacity_multiplier(double nutrient) const noexcept;
    double mean_growth_rate(CellType type) const noexcept;
    double base_diffusion(PopulationField3D field,
                          double occupied_fraction) const noexcept;
    void validate_state() const;

    ContinuumModelConfig3D config_;
    std::size_t voxel_count_{};
    double voxel_measure_{};
    double large_cell_volume_{};
    std::array<std::vector<double>, kPopulationFieldCount3D> populations_;
    std::array<std::vector<double>, kPopulationFieldCount3D> work_;
    std::vector<double> nutrient_;
    std::vector<double> nutrient_next_;
    std::vector<double> vessel_;
    double time_hours_{};
    double next_nutrient_refresh_hours_{};
    std::uint64_t step_count_{};
    std::uint64_t nutrient_solve_count_{};
    bool initialized_{};
};

}  // namespace atcg3d::continuum
