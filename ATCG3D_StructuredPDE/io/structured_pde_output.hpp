#pragma once

#include <filesystem>
#include <fstream>

#include "model/structured_pde_model.hpp"

namespace atcg3d::structured_pde {

class StructuredPdeOutput3D {
public:
    StructuredPdeOutput3D(const StructuredPdeConfig3D& config,
                          double initial_time_hours);
    void observe(const StructuredPdeModel3D& model);
    void checkpoint_now(const StructuredPdeModel3D& model);
    void finalize(const StructuredPdeModel3D& model);

private:
    static bool due(double now, double next, double interval) noexcept;
    static double advance(double next, double now, double interval) noexcept;
    void write_metrics(const StructuredPdeModel3D& model);
    void write_field(const StructuredPdeModel3D& model);
    void write_checkpoint(const StructuredPdeModel3D& model);

    StructuredPdeConfig3D config_;
    std::filesystem::path directory_;
    std::ofstream metrics_;
    double next_metrics_{};
    double next_field_{};
    double next_checkpoint_{};
    double last_checkpoint_time_{-1.0};
    std::size_t field_index_{};
    bool finalized_{};
};

}  // namespace atcg3d::structured_pde
