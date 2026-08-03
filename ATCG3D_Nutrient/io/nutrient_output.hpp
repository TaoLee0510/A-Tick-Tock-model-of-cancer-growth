#pragma once

#include <filesystem>
#include <fstream>

#include "config/nutrient_config.hpp"
#include "field/nutrient_field.hpp"

namespace atcg3d {
class Simulation3D;
struct SimulationClock3D;
}

namespace atcg3d::nutrient {

class NutrientOutput3D {
public:
    NutrientOutput3D(const NutrientModelConfig3D& config,
                     NutrientEnvironment3D& environment,
                     double resume_time_hours = -1.0);

    void observe(const Simulation3D& simulation);
    void checkpoint_if_due(const Simulation3D& simulation);
    void checkpoint_now(const Simulation3D& simulation);

    static std::filesystem::path sidecar_path(
        const std::filesystem::path& base_checkpoint);
    static std::filesystem::path base_checkpoint_path(
        const std::filesystem::path& output_directory,
        const SimulationClock3D& clock);

private:
    static bool same_time(double lhs, double rhs) noexcept;
    static bool due(double now, double next, double interval) noexcept;
    static double advance(double next, double now, double interval) noexcept;
    void write_metrics(const Simulation3D& simulation);
    void write_field_snapshot(const Simulation3D& simulation);
    void write_sidecar(const Simulation3D& simulation);

    NutrientModelConfig3D config_;
    NutrientEnvironment3D& environment_;
    std::filesystem::path directory_;
    std::ofstream metrics_;
    double next_metrics_{};
    double next_field_snapshot_{};
    double next_checkpoint_{};
    double last_sidecar_time_{-1.0};
    std::size_t field_snapshot_index_{};
};

}  // namespace atcg3d::nutrient
