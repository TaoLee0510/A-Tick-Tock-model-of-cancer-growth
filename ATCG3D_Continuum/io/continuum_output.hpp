#pragma once

#include <cstddef>
#include <filesystem>
#include <fstream>

#include "model/continuum_model.hpp"

namespace atcg3d::continuum {

class ContinuumOutput3D {
public:
    ContinuumOutput3D(const ContinuumModelConfig3D& config,
                      double initial_time_hours);

    void observe(const ContinuumModel3D& model);
    void finalize(const ContinuumModel3D& model);
    void checkpoint_now(const ContinuumModel3D& model);

private:
    static bool same_time(double lhs, double rhs) noexcept;
    static bool due(double now, double next, double interval) noexcept;
    static double advance(double next, double now, double interval) noexcept;
    static std::filesystem::path checkpoint_path(
        const std::filesystem::path& directory,
        const ContinuumModel3D& model);
    void write_metrics(const ContinuumModel3D& model);
    void write_field(const ContinuumModel3D& model);
    void write_radial_profile(const ContinuumModel3D& model);
    void write_checkpoint(const ContinuumModel3D& model);

    ContinuumModelConfig3D config_;
    std::filesystem::path directory_;
    std::ofstream metrics_;
    double next_metrics_{};
    double next_field_{};
    double next_profile_{};
    double next_checkpoint_{};
    double last_checkpoint_time_{-1.0};
    std::size_t field_index_{};
    std::size_t profile_index_{};
    bool finalized_{};
};

}  // namespace atcg3d::continuum
