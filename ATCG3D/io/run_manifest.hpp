#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "config/model_config.hpp"

namespace atcg3d {

struct SeriesEntry3D {
    std::string name;
    double time_hours{};
};

struct ExistingRunOutput3D {
    std::vector<SeriesEntry3D> preview;
    std::vector<SeriesEntry3D> full;
    std::vector<SeriesEntry3D> vessels;
};

void write_series_atomic(const std::filesystem::path& path,
                         const std::vector<SeriesEntry3D>& entries);
void write_run_manifest_atomic(const std::filesystem::path& run_directory,
                               const Model3DConfig& config,
                               const std::vector<SeriesEntry3D>& preview,
                               const std::vector<SeriesEntry3D>& full,
                               const std::vector<SeriesEntry3D>& vessels,
                               bool checkpoint_enabled);

// Loads a run created by this schema for append-only resume. The loader is
// deliberately strict: catalogs must be ordered, use the canonical relative
// frame names, reference existing files, and describe the same dynamics as
// the checkpoint-validated resume configuration.
ExistingRunOutput3D load_existing_run_output(
    const std::filesystem::path& run_directory,
    const Model3DConfig& resume_config);

}  // namespace atcg3d
