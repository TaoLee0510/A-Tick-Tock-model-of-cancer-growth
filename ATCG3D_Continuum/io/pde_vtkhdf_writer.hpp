#pragma once

#include "config/continuum_config.hpp"
#include <filesystem>
#include <functional>
#include <string>
#include <vector>

namespace atcg3d::structured_pde {
class StructuredPdeModel3D;
}
namespace atcg3d::continuum {
class ContinuumModel3D;
struct GridFieldView3D {
    std::string name;
    std::function<double(std::size_t)> value;
};
bool pde_vtkhdf_available() noexcept;
void write_grid_vtkhdf(const std::filesystem::path &path,
                       const ContinuumGridConfig3D &grid, double time,
                       const std::vector<GridFieldView3D> &fields);
void append_pde_series(const std::filesystem::path &directory,
                       const std::filesystem::path &path, double time);
void write_pde_vtkhdf(const std::filesystem::path &path,
                      const ContinuumModel3D &model);
void write_pde_vtkhdf(const std::filesystem::path &path,
                      const structured_pde::StructuredPdeModel3D &model);
} // namespace atcg3d::continuum
