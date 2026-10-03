#include "io/pde_vtkhdf_writer.hpp"
#include "model/continuum_model.hpp"
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <yaml-cpp/yaml.h>
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
#include <H5Cpp.h>
#endif

namespace atcg3d::continuum {
bool pde_vtkhdf_available() noexcept {
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
    return true;
#else
    return false;
#endif
}
void write_grid_vtkhdf(const std::filesystem::path &path,
                       const ContinuumGridConfig3D &grid, double time,
                       const std::vector<GridFieldView3D> &fields) {
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
    if (std::filesystem::exists(path))
        throw std::runtime_error("refusing to overwrite PDE VTK-HDF");
    if (!path.parent_path().empty())
        std::filesystem::create_directories(path.parent_path());
    const auto temporary = path.string() + ".tmp";
    try {
        {
            H5::H5File file(temporary, H5F_ACC_TRUNC);
            auto group = file.createGroup("VTKHDF");
            auto attr = [&](const char *name, const H5::DataType &type,
                            hsize_t n, const void *value) {
                H5::DataSpace space(1, &n);
                auto a = group.createAttribute(name, type, space);
                a.write(type, value);
            };
            H5::StrType string_type(H5::PredType::C_S1, 9);
            H5::DataSpace scalar(H5S_SCALAR);
            auto type = group.createAttribute("Type", string_type, scalar);
            type.write(string_type, "ImageData");
            const int version[2]{2, 0};
            attr("Version", H5::PredType::NATIVE_INT, 2, version);
            const int extent[6]{0, grid.shape[0] - 1, 0, grid.shape[1] - 1,
                                0, grid.shape[2] - 1};
            attr("WholeExtent", H5::PredType::NATIVE_INT, 6, extent);
            double origin[3],
                spacing[3]{grid.spacing_voxels, grid.spacing_voxels,
                           grid.spacing_voxels};
            for (int i = 0; i < 3; ++i)
                origin[i] = grid.origin[i] + grid.spacing_voxels / 2;
            attr("Origin", H5::PredType::NATIVE_DOUBLE, 3, origin);
            attr("Spacing", H5::PredType::NATIVE_DOUBLE, 3, spacing);
            const double direction[9]{1, 0, 0, 0, 1, 0, 0, 0, 1};
            attr("Direction", H5::PredType::NATIVE_DOUBLE, 9, direction);
            auto point = group.createGroup("PointData");
            group.createGroup("CellData");
            auto metadata = group.createGroup("FieldData");
            hsize_t one = 1;
            H5::DataSpace time_space(1, &one);
            auto time_data = metadata.createDataSet(
                "time_hours", H5::PredType::NATIVE_DOUBLE, time_space);
            time_data.write(&time, H5::PredType::NATIVE_DOUBLE);
            H5::StrType scalar_name(H5::PredType::C_S1, 7);
            auto scalars =
                point.createAttribute("Scalars", scalar_name, scalar);
            scalars.write(scalar_name, "r_total");
            // VTK omits singleton image axes from scalar dataset rank.
            std::vector<int> axes;
            std::vector<hsize_t> shape, chunk, row_shape;
            for (int axis = 2; axis >= 0; --axis)
                if (grid.shape[axis] > 1) {
                    axes.push_back(axis);
                    shape.push_back(grid.shape[axis]);
                    chunk.push_back(axis == 2 ? 1
                                              : std::min(64, grid.shape[axis]));
                    row_shape.push_back(axis == 0 ? grid.shape[0] : 1);
                }
            if (axes.empty()) {
                axes.push_back(0);
                shape.push_back(1);
                chunk.push_back(1);
                row_shape.push_back(1);
            }
            H5::DataSpace space(int(shape.size()), shape.data());
            H5::DSetCreatPropList properties;
            properties.setChunk(int(chunk.size()), chunk.data());
            properties.setDeflate(1);
            H5::DataSpace memory(int(row_shape.size()), row_shape.data());
            std::vector<double> row(grid.shape[0]);
            for (const auto &field : fields) {
                auto data = point.createDataSet(
                    field.name, H5::PredType::NATIVE_DOUBLE, space, properties);
                for (int z = 0; z < grid.shape[2]; ++z)
                    for (int y = 0; y < grid.shape[1]; ++y) {
                        for (int x = 0; x < grid.shape[0]; ++x)
                            row[x] = field.value(
                                (std::size_t(z) * grid.shape[1] + y) *
                                    grid.shape[0] +
                                x);
                        std::vector<hsize_t> start;
                        for (int axis : axes)
                            start.push_back(axis == 2 ? z : axis == 1 ? y : 0);
                        auto selection = data.getSpace();
                        selection.selectHyperslab(
                            H5S_SELECT_SET, row_shape.data(), start.data());
                        data.write(row.data(), H5::PredType::NATIVE_DOUBLE,
                                   memory, selection);
                    }
            }

            file.flush(H5F_SCOPE_GLOBAL);
        }
        std::filesystem::rename(temporary, path);
    } catch (...) {
        std::filesystem::remove(temporary);
        throw;
    }
#else
    (void)path;
    (void)grid;
    (void)time;
    (void)fields;
    throw std::runtime_error(
        "PDE VTK-HDF requires ATCG3D_ENABLE_HDF5_CHECKPOINT=ON");
#endif
}
void append_pde_series(const std::filesystem::path &directory,
                       const std::filesystem::path &path, double time) {
    const auto series = directory / "fields.vtkhdf.series";
    YAML::Node files;
    if (std::filesystem::exists(series))
        files = YAML::LoadFile(series.string())["files"];
    const auto temporary = series.string() + ".tmp";
    std::ofstream out(temporary);
    out << std::setprecision(17)
        << "{\"file-series-version\":\"1.0\",\"files\":[";
    bool first = true;
    const auto entry = [&](std::string name, double t) {
        if (!first)
            out << ',';
        first = false;
        out << "{\"name\":\"" << name << "\",\"time\":" << t << '}';
    };
    for (const auto &file : files)
        entry(file["name"].as<std::string>(), file["time"].as<double>());
    entry(std::filesystem::relative(path, directory).generic_string(), time);
    out << "]}\n";
    out.close();
    if (!out)
        throw std::runtime_error("PDE series write failed");
    std::filesystem::rename(temporary, series);
}
void write_pde_vtkhdf(const std::filesystem::path &path,
                      const ContinuumModel3D &model) {
    std::vector<GridFieldView3D> fields;
    for (auto [name, id] :
         std::vector<std::pair<std::string, PopulationField3D>>{
             {"r_small", PopulationField3D::r_small},
             {"r_large", PopulationField3D::r_large},
             {"K_small", PopulationField3D::K_small},
             {"K_large", PopulationField3D::K_large}})
        fields.push_back(
            {name, [&model, id](auto i) { return model.population(id)[i]; }});
    fields.push_back(
        {"r_total", [&](auto i) {
             return model.population(PopulationField3D::r_small)[i] +
                    model.population(PopulationField3D::r_large)[i];
         }});
    fields.push_back(
        {"K_total", [&](auto i) {
             return model.population(PopulationField3D::K_small)[i] +
                    model.population(PopulationField3D::K_large)[i];
         }});
    fields.push_back({"nutrient", [&](auto i) { return model.nutrient()[i]; }});
    fields.push_back({"vessel_fraction",
                      [&](auto i) { return model.vessel_fraction()[i]; }});
    fields.push_back({"occupied_fraction",
                      [&](auto i) { return model.occupied_fraction(i); }});
    if (model.angiogenesis()) {
        fields.push_back(
            {"VEGF", [&](auto i) { return model.angiogenesis()->taf()[i]; }});
        fields.push_back({"vessel_tips", [&](auto i) {
                              return model.angiogenesis()->tips()[i];
                          }});
    }
    write_grid_vtkhdf(path, model.config().grid, model.time_hours(), fields);
}
} // namespace atcg3d::continuum
