#include "io/vtkhdf_writer.hpp"

#include <stdexcept>

namespace atcg3d {

bool vtkhdf_output_available() noexcept {
    return false;
}

void write_vtkhdf_points_atomic(const std::filesystem::path&,
                                const SimulationSnapshotView3D&) {
    throw std::runtime_error(
        "VTK-HDF output is unavailable: configure with "
        "-DATCG3D_ENABLE_VTKHDF=ON and install VTK with IOHDF");
}

std::size_t read_vtkhdf_point_count(const std::filesystem::path&) {
    throw std::runtime_error("VTK-HDF reader is unavailable in this build");
}

}  // namespace atcg3d
