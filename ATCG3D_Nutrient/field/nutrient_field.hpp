#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <map>
#include <vector>

#include "config/nutrient_config.hpp"
#include "engine/environment.hpp"

namespace atcg3d {
class CellStore3D;
class SparseVesselGrid3D;
}

namespace atcg3d::nutrient {

struct NutrientFieldDiagnostics3D {
    std::size_t block_count{};
    std::size_t active_voxel_count{};
    std::size_t perfused_source_voxels{};
    std::size_t consuming_voxels{};
    double minimum{};
    double maximum{};
    double mean{};
    double last_max_update{};
};

struct NutrientVoxelSample3D {
    Vec3i site{};
    float value{};
};

class NutrientEnvironment3D final : public EnvironmentCoupling3D {
public:
    NutrientEnvironment3D(NutrientFieldConfig3D config, bool thin_layer);

    EnvironmentInitializationResult3D initialize(
        double now_hours,
        const CellStore3D& cells,
        const SparseVesselGrid3D& vessels) override;
    void refresh(double now_hours,
                 const CellStore3D& cells,
                 const SparseVesselGrid3D& vessels) override;

    double retained_density(Vec3i site) const noexcept override;
    double next_refresh_time_hours() const noexcept override {
        return next_refresh_time_hours_;
    }
    std::uint32_t schedule_generation() const noexcept override {
        return schedule_generation_;
    }
    std::uint64_t refresh_count() const noexcept override {
        return refresh_count_;
    }
    std::size_t allocated_bytes() const noexcept override;
    std::uint64_t field_checksum() const noexcept override;

    float value(Vec3i site) const noexcept;
    double capacity_multiplier(Vec3i site) const noexcept;
    double last_refresh_time_hours() const noexcept {
        return last_refresh_time_hours_;
    }
    const NutrientFieldConfig3D& config() const noexcept { return config_; }
    const NutrientFieldDiagnostics3D& diagnostics() const noexcept {
        return diagnostics_;
    }
    std::vector<NutrientVoxelSample3D> nonzero_voxels() const;

    void save_checkpoint(const std::filesystem::path& path,
                         std::uint64_t base_state_checksum,
                         double time_hours,
                         std::uint64_t completed_events) const;
    void load_checkpoint(const std::filesystem::path& path,
                         std::uint64_t expected_base_state_checksum,
                         double expected_time_hours,
                         std::uint64_t expected_completed_events);

private:
    struct Block {
        explicit Block(std::size_t count = 0)
            : value(count, 0.0F), next(count, 0.0F),
              r_consumption(count, 0.0F), K_consumption(count, 0.0F),
              perfused(count, 0U) {}
        std::vector<float> value;
        std::vector<float> next;
        std::vector<float> r_consumption;
        std::vector<float> K_consumption;
        std::vector<std::uint8_t> perfused;
    };

    struct Address {
        Vec3i block{};
        std::uint32_t index{};
    };

    static int floor_div(int value, int divisor) noexcept;
    Address address(Vec3i site) const noexcept;
    Vec3i site_from(Vec3i block, std::uint32_t index) const noexcept;
    Block* find_block(Vec3i coordinate) noexcept;
    const Block* find_block(Vec3i coordinate) const noexcept;
    bool active_site(Vec3i site) const noexcept;
    void rebuild_sources_and_sinks(const CellStore3D& cells,
                                   const SparseVesselGrid3D& vessels);
    void solve();
    void update_diagnostics();

    NutrientFieldConfig3D config_;
    bool thin_layer_{};
    std::size_t block_voxels_{};
    std::map<Vec3i, Block> blocks_;
    double last_refresh_time_hours_{};
    double next_refresh_time_hours_{};
    std::uint32_t schedule_generation_{};
    std::uint64_t refresh_count_{};
    NutrientFieldDiagnostics3D diagnostics_;
    bool loaded_checkpoint_{};
    double loaded_checkpoint_time_{};
};

}  // namespace atcg3d::nutrient
