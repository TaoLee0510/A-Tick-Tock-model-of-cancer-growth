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

void write_series_atomic(const std::filesystem::path& path,
                         const std::vector<SeriesEntry3D>& entries);
void write_run_manifest_atomic(const std::filesystem::path& run_directory,
                               const Model3DConfig& config,
                               const std::vector<SeriesEntry3D>& preview,
                               const std::vector<SeriesEntry3D>& full,
                               bool checkpoint_enabled);

}  // namespace atcg3d
