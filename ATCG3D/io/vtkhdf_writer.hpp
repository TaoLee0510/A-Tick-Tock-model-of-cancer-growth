#pragma once

#include <cstddef>
#include <filesystem>

#include "io/snapshot.hpp"
#include "vasculature/vessel_store.hpp"

namespace atcg3d {

bool vtkhdf_output_available() noexcept;

// Writes one biological cell as one vtkPolyData point. The final file appears
// only after vtkHDFWriter has closed the temporary file successfully.
void write_vtkhdf_points_atomic(const std::filesystem::path& path,
                                const SimulationSnapshotView3D& snapshot,
                                int compression_level = 1);
void write_vtkhdf_points_atomic(const std::filesystem::path& path,
                                const FrozenCellSnapshotView3D& snapshot,
                                int compression_level = 1);

// Writes the vascular centerline as a separate vtkPolyData. Every live vessel
// node becomes one point and every non-root node contributes one parent-child
// line. Vessel geometry is deliberately not mixed into the biological-cell
// points-only files.
void write_vtkhdf_vessels_atomic(const std::filesystem::path& path,
                                 const VesselNodeStore3D& nodes,
                                 int compression_level = 1);
void write_vtkhdf_vessels_atomic(
    const std::filesystem::path& path,
    std::span<const VesselNodeInit3D> nodes,
    int compression_level = 1);

std::size_t read_vtkhdf_point_count(const std::filesystem::path& path);

}  // namespace atcg3d
