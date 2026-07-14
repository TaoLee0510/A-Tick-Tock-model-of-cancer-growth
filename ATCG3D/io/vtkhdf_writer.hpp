#pragma once

#include <cstddef>
#include <filesystem>

#include "io/snapshot.hpp"

namespace atcg3d {

bool vtkhdf_output_available() noexcept;

// Writes one biological cell as one vtkPolyData point. The final file appears
// only after vtkHDFWriter has closed the temporary file successfully.
void write_vtkhdf_points_atomic(const std::filesystem::path& path,
                                const SimulationSnapshotView3D& snapshot);
std::size_t read_vtkhdf_point_count(const std::filesystem::path& path);

}  // namespace atcg3d
