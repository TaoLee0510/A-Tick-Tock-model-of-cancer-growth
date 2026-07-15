#include "io/vtkhdf_writer.hpp"

#include <filesystem>
#include <limits>
#include <stdexcept>
#include <vector>

#include <vtkCellArray.h>
#include <vtkFieldData.h>
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

template <class Builder>
void write_polydata_atomic(const std::filesystem::path& path, Builder builder) {
    std::filesystem::create_directories(path.parent_path());
    const std::filesystem::path temporary = path.string() + ".tmp";
    std::filesystem::remove(temporary);

    try {
        vtkSmartPointer<vtkPolyData> data = builder();
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

}  // namespace

bool vtkhdf_output_available() noexcept {
    return true;
}

void write_vtkhdf_points_atomic(const std::filesystem::path& path,
                                const SimulationSnapshotView3D& snapshot) {
    if (snapshot.size() > static_cast<std::size_t>(std::numeric_limits<vtkIdType>::max())) {
        throw std::overflow_error("snapshot has more points than this VTK build can index");
    }
    write_polydata_atomic(path, [&snapshot]() {
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
        auto total_cell_count =
            make_array<vtkUnsignedLongLongArray>("total_cell_count", 1);
        total_cell_count->SetValue(
            0, static_cast<unsigned long long>(snapshot.cells.alive_count()));

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

        vtkSmartPointer<vtkPolyData> data = vtkSmartPointer<vtkPolyData>::New();
        data->SetPoints(points);
        data->GetPointData()->AddArray(cell_id);
        data->GetPointData()->AddArray(clone_id);
        data->GetPointData()->AddArray(cell_type);
        data->GetPointData()->AddArray(stage);
        data->GetPointData()->AddArray(viability);
        data->GetPointData()->AddArray(display_radius);
        data->GetFieldData()->AddArray(total_cell_count);

        return data;
    });
}

void write_vtkhdf_vessels_atomic(const std::filesystem::path& path,
                                 const VesselNodeStore3D& nodes) {
    if (nodes.alive_count() >
        static_cast<std::size_t>(std::numeric_limits<vtkIdType>::max())) {
        throw std::overflow_error("vascular snapshot has more nodes than this VTK build can index");
    }

    write_polydata_atomic(path, [&nodes]() {
        const std::vector<VesselNodeSlot> slots = nodes.alive_slots();
        const vtkIdType count = static_cast<vtkIdType>(slots.size());
        vtkNew<vtkPoints> points;
        points->SetDataTypeToFloat();
        points->SetNumberOfPoints(count);
        auto node_id = make_array<vtkUnsignedLongLongArray>("node_id", count);
        auto vessel_id = make_array<vtkUnsignedLongLongArray>("vessel_id", count);
        auto branch_role = make_array<vtkUnsignedCharArray>("branch_role", count);
        auto perfused = make_array<vtkUnsignedCharArray>("perfused", count);
        auto diameter_voxels = make_array<vtkFloatArray>("diameter_voxels", count);
        auto radius_voxels = make_array<vtkFloatArray>("radius_voxels", count);

        std::vector<vtkIdType> point_for_slot(nodes.slot_count(), static_cast<vtkIdType>(-1));
        for (vtkIdType index = 0; index < count; ++index) {
            const VesselNodeSlot slot = slots[static_cast<std::size_t>(index)];
            point_for_slot[slot] = index;
            const Vec3i position = nodes.position(slot);
            // Cell snapshots use voxel-corner coordinates: a one-voxel cell
            // at anchor p is displayed at p+0.5. Vessel nodes occupy the same
            // lattice voxels, so use the identical centre convention.
            points->SetPoint(index, static_cast<float>(position.x) + 0.5F,
                             static_cast<float>(position.y) + 0.5F,
                             static_cast<float>(position.z) + 0.5F);
            node_id->SetValue(index, nodes.uid(slot));
            vessel_id->SetValue(index, nodes.vessel_id(slot));
            branch_role->SetValue(index, static_cast<unsigned char>(nodes.role(slot)));
            perfused->SetValue(index, nodes.perfused(slot) ? 1U : 0U);
            diameter_voxels->SetValue(index, nodes.diameter_voxels(slot));
            radius_voxels->SetValue(index, nodes.diameter_voxels(slot) * 0.5F);
        }

        vtkNew<vtkCellArray> lines;
        lines->AllocateEstimate(count, 2);
        for (const VesselNodeSlot slot : slots) {
            const VesselNodeSlot parent = nodes.parent_node_slot(slot);
            if (parent == kEmptyVesselNodeSlot) continue;
            if (parent >= point_for_slot.size() || point_for_slot[parent] < 0) {
                throw std::runtime_error("vascular snapshot contains a missing parent node");
            }
            const vtkIdType endpoints[2] = {point_for_slot[parent], point_for_slot[slot]};
            lines->InsertNextCell(2, endpoints);
        }

        vtkSmartPointer<vtkPolyData> data = vtkSmartPointer<vtkPolyData>::New();
        data->SetPoints(points);
        data->SetLines(lines);
        data->GetPointData()->AddArray(node_id);
        data->GetPointData()->AddArray(vessel_id);
        data->GetPointData()->AddArray(branch_role);
        data->GetPointData()->AddArray(perfused);
        data->GetPointData()->AddArray(diameter_voxels);
        data->GetPointData()->AddArray(radius_voxels);
        return data;
    });
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
