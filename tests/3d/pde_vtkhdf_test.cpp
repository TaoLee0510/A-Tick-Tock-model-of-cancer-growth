#include "engine/simulation.hpp"
#include "io/pde_vtkhdf_writer.hpp"
#include "model/structured_pde_model.hpp"
#include <H5Cpp.h>
#include <cassert>
#include <filesystem>

int main() {
    using namespace atcg3d;
    using namespace atcg3d::continuum;
    auto c = structured_pde::StructuredPdeConfig3D::load(
        std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_SharedRules/config/angiogenesis_v8.yaml");
    c.continuum.base.angiogenesis.enabled = false;
    Simulation3D source(c.continuum.base);
    source.initialize();
    structured_pde::StructuredPdeModel3D model(c);
    model.initialize_from_abm(source);
    model.step();
    const auto directory = std::filesystem::current_path() / "pde-vtkhdf-test";
    std::filesystem::remove_all(directory);
    const auto path = directory / "fields/field.vtkhdf";
    write_pde_vtkhdf(path, model);
    append_pde_series(directory, path, model.time_hours());
    H5::H5File file(path.string(), H5F_ACC_RDONLY);
    auto group = file.openGroup("VTKHDF");
    int extent[6];
    group.openAttribute("WholeExtent").read(H5::PredType::NATIVE_INT, extent);
    assert(extent[1] == 47 && extent[3] == 47 && extent[5] == 0);
    for (auto name :
         {"r_total", "r_normal_small", "r_active_large", "r_refractory_small",
          "K_total", "nutrient", "vessel_fraction", "occupied_fraction", "VEGF",
          "vessel_tips"})
        assert(H5Lexists(group.getId(),
                         ("PointData/" + std::string(name)).c_str(),
                         H5P_DEFAULT) > 0);
    auto data = file.openDataSet("VTKHDF/PointData/nutrient");
    hsize_t shape[3];
    assert(data.getSpace().getSimpleExtentNdims() == 2);
    data.getSpace().getSimpleExtentDims(shape);
    assert(shape[0] == 48 && shape[1] == 48);
    std::vector<double> values(2304);
    data.read(values.data(), H5::PredType::NATIVE_DOUBLE);
    assert(values == model.nutrient());
    double origin[3];
    group.openAttribute("Origin").read(H5::PredType::NATIVE_DOUBLE, origin);
    assert(origin[0] == -23.5 && origin[2] == 0);
    assert(std::filesystem::exists(directory / "fields.vtkhdf.series"));
    assert(!std::filesystem::exists(path.string() + ".tmp"));
    ContinuumGridConfig3D volume;
    volume.shape = {2, 3, 4};
    volume.origin = {0, 0, 0};
    volume.spacing_voxels = 1;
    write_grid_vtkhdf(directory / "volume.vtkhdf", volume, 0,
                      {{"r_total", [](auto i) { return double(i); }}});
    // Keep this tiny fixture in the build directory for external VTK reader QA.
}
