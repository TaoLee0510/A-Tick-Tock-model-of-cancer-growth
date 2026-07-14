#include "io/vtkhdf_writer.hpp"

#include <filesystem>
#include <limits>
#include <stdexcept>

#include <vtkFloatArray.h>
#include <vtkHDFReader.h>
#include <vtkHDFWriter.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkUnsignedLongLongArray.h>

namespace atcg3d {
namespace {

template <class Array>
vtkSmartPointer<Array> make_array(const char* name, vtkIdType count) {
    vtkSmartPointer<Array> result = vtkSmartPointer<Array>::New();
    result->SetName(name);
    result->SetNumberOfComponents(1);
    result->SetNumberOfTuples(count);
    return result;
}

}  // namespace

bool vtkhdf_output_available() noexcept {
    return true;
}

void write_vtkhdf_points_atomic(const std::filesystem::path& path,
                                const SimulationSnapshotView3D& snapshot) {
    if (snapshot.size() > static_cast<std::size_t>(std::numeric_limits<vtkIdType>::max())) {
        throw std::overflow_error("snapshot has more points than this VTK build can index");
    }
    std::filesystem::create_directories(path.parent_path());
    const std::filesystem::path temporary = path.string() + ".tmp";
    std::filesystem::remove(temporary);

    try {
        const vtkIdType count = static_cast<vtkIdType>(snapshot.size());
        vtkNew<vtkPoints> points;
        points->SetDataTypeToFloat();
        points->SetNumberOfPoints(count);
        auto cell_id = make_array<vtkUnsignedLongLongArray>("cell_id", count);
        auto clone_id = make_array<vtkUnsignedIntArray>("clone_id", count);
        auto cell_type = make_array<vtkUnsignedCharArray>("cell_type", count);
        auto stage = make_array<vtkUnsignedCharArray>("stage", count);
        auto viability = make_array<vtkUnsignedCharArray>("viability", count);
        auto display_radius = make_array<vtkFloatArray>("display_radius", count);

        for (vtkIdType index = 0; index < count; ++index) {
            const std::size_t source = static_cast<std::size_t>(index);
            const Slot slot = snapshot.slot(source);
            const auto center = snapshot.center(source);
            points->SetPoint(index, center[0], center[1], center[2]);
            cell_id->SetValue(index, snapshot.cells.uid(slot));
            clone_id->SetValue(index, snapshot.cells.clone_id(slot));
            cell_type->SetValue(index, static_cast<unsigned char>(snapshot.cells.type(slot)));
            stage->SetValue(index, static_cast<unsigned char>(snapshot.cells.stage(slot)));
            viability->SetValue(index, snapshot.cells.viability(slot));
            display_radius->SetValue(index, snapshot.display_radius(source));
        }

        vtkNew<vtkPolyData> data;
        data->SetPoints(points);
        data->GetPointData()->AddArray(cell_id);
        data->GetPointData()->AddArray(clone_id);
        data->GetPointData()->AddArray(cell_type);
        data->GetPointData()->AddArray(stage);
        data->GetPointData()->AddArray(viability);
        data->GetPointData()->AddArray(display_radius);

        {
            vtkNew<vtkHDFWriter> writer;
            writer->SetFileName(temporary.string().c_str());
            writer->SetInputData(data);
            if (writer->Write() != 1) {
                throw std::runtime_error("vtkHDFWriter failed for " + temporary.string());
            }
        }
        if (!std::filesystem::is_regular_file(temporary) ||
            std::filesystem::file_size(temporary) == 0) {
            throw std::runtime_error("vtkHDFWriter produced no data for " + temporary.string());
        }
        std::filesystem::rename(temporary, path);
    } catch (...) {
        std::filesystem::remove(temporary);
        throw;
    }
}

std::size_t read_vtkhdf_point_count(const std::filesystem::path& path) {
    vtkNew<vtkHDFReader> reader;
    reader->SetFileName(path.string().c_str());
    reader->Update();
    vtkPolyData* data = vtkPolyData::SafeDownCast(reader->GetOutputDataObject(0));
    if (data == nullptr) {
        throw std::runtime_error("VTK-HDF file is not point vtkPolyData: " + path.string());
    }
    if (data->GetNumberOfCells() != 0) {
        throw std::runtime_error("VTK-HDF points-only file unexpectedly contains cells");
    }
    return static_cast<std::size_t>(data->GetNumberOfPoints());
}

}  // namespace atcg3d
