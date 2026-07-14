#pragma once

#include <cstddef>
#include <filesystem>
#include <vector>

#include "config/model_config.hpp"
#include "engine/simulation.hpp"
#include "io/run_manifest.hpp"

namespace atcg3d {

class OutputManager3D {
public:
    explicit OutputManager3D(const Model3DConfig& config);

    void observe(const Simulation3D& simulation);
    void finalize(const Simulation3D& simulation);

    const std::vector<SeriesEntry3D>& preview_entries() const noexcept { return preview_; }
    const std::vector<SeriesEntry3D>& full_entries() const noexcept { return full_; }

private:
    bool due(double now, double next, double interval) const noexcept;
    static double advance(double next, double now, double interval) noexcept;
    void write_preview(const Simulation3D& simulation);
    void write_full(const Simulation3D& simulation);
    void write_checkpoint(const Simulation3D& simulation);
    void append_lineage(const Simulation3D& simulation);
    void update_metadata();
    void write_metrics(const Simulation3D& simulation);

    Model3DConfig config_;
    std::filesystem::path run_directory_;
    std::vector<SeriesEntry3D> preview_;
    std::vector<SeriesEntry3D> full_;
    double next_preview_{};
    double next_full_{};
    double next_checkpoint_{};
    double last_preview_time_{-1.0};
    double last_full_time_{-1.0};
    std::size_t lineage_written_{};
    bool initialized_{};
    bool finalized_{};
};

}  // namespace atcg3d
