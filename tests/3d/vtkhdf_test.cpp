#include <cassert>
#include <filesystem>
#include <vector>

#include <vtkDataArray.h>
#include <vtkFieldData.h>
#include <vtkHDFReader.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>

#include "core/cell_store.hpp"
#include "config/model_config.hpp"
#include "engine/simulation.hpp"
#include "io/output_manager.hpp"
#include "io/snapshot.hpp"
#include "io/vtkhdf_writer.hpp"
#include "vasculature/vessel_store.hpp"

int main() {
    using namespace atcg3d;
    CellStore3D cells;
    CellInit large;
    large.uid = 42;
    large.clone_id = 7;
    large.anchor = {-2, 3, 4};
    large.stage = CellStage::large;
    const Slot first = cells.create(large);
    CellInit small;
    small.uid = 43;
    small.anchor = {8, -1, 2};
    small.type = CellType::K;
    const Slot second = cells.create(small);
    const std::vector<Slot> slots{first, second};
    DisplayRadiusConfig radii;
    radii.large = 1.75;
    radii.small = 0.625;
    radii.ultrasmall = 0.2;
    const SimulationSnapshotView3D snapshot{cells, slots, {3, 1.5}, radii};

    const auto path = std::filesystem::temp_directory_path() / "atcg3d_vtkhdf_test.vtkhdf";
    std::filesystem::remove(path);
    std::filesystem::remove(path.string() + ".tmp");
    write_vtkhdf_points_atomic(path, snapshot);
    assert(std::filesystem::is_regular_file(path));
    assert(!std::filesystem::exists(path.string() + ".tmp"));

    vtkNew<vtkHDFReader> reader;
    reader->SetFileName(path.string().c_str());
    reader->Update();
    vtkPolyData* data = vtkPolyData::SafeDownCast(reader->GetOutputDataObject(0));
    assert(data != nullptr);
    assert(data->GetNumberOfPoints() == 2);
    assert(data->GetNumberOfCells() == 0);
    assert(data->GetPointData()->GetArray("cell_id")->GetDataTypeSize() == 8);
    assert(data->GetPointData()->GetArray("clone_id")->GetDataTypeSize() == 4);
    assert(data->GetPointData()->GetArray("cell_type")->GetDataTypeSize() == 1);
    assert(data->GetPointData()->GetArray("stage")->GetDataTypeSize() == 1);
    assert(data->GetPointData()->GetArray("viability")->GetDataTypeSize() == 1);
    assert(data->GetPointData()->GetArray("display_radius")->GetDataTypeSize() == 4);
    assert(data->GetPointData()->GetArray("display_radius")->GetTuple1(0) == 1.75);
    assert(data->GetPointData()->GetArray("display_radius")->GetTuple1(1) == 0.625);
    assert(data->GetFieldData()->GetArray("total_cell_count") != nullptr);
    assert(data->GetFieldData()->GetArray("total_cell_count")->GetDataTypeSize() == 8);
    assert(data->GetFieldData()->GetArray("total_cell_count")->GetTuple1(0) == 2.0);
    double point[3]{};
    data->GetPoint(0, point);
    assert(point[0] == -1.0 && point[1] == 4.0 && point[2] == 5.0);
    data->GetPoint(1, point);
    assert(point[0] == 8.5 && point[1] == -0.5 && point[2] == 2.5);
    std::filesystem::remove(path);

    VesselNodeStore3D vessels;
    VesselNodeInit3D root;
    root.position = {0, 0, 0};
    root.uid = 101;
    root.vessel_id = 77;
    root.role = VesselBranchRole::root;
    root.perfused = true;
    root.diameter_voxels = 2.0F;
    const VesselNodeSlot root_slot = vessels.create(root);
    VesselNodeInit3D child;
    child.position = {1, 1, 0};
    child.uid = 102;
    child.parent_uid = root.uid;
    child.parent_node_slot = root_slot;
    child.vessel_id = root.vessel_id;
    child.role = VesselBranchRole::outward;
    child.perfused = false;
    child.diameter_voxels = 3.0F;
    vessels.create(child);

    const auto vessel_path =
        std::filesystem::temp_directory_path() / "atcg3d_vessels_vtkhdf_test.vtkhdf";
    std::filesystem::remove(vessel_path);
    std::filesystem::remove(vessel_path.string() + ".tmp");
    write_vtkhdf_vessels_atomic(vessel_path, vessels);
    assert(std::filesystem::is_regular_file(vessel_path));
    assert(!std::filesystem::exists(vessel_path.string() + ".tmp"));
    vtkNew<vtkHDFReader> vessel_reader;
    vessel_reader->SetFileName(vessel_path.string().c_str());
    vessel_reader->Update();
    vtkPolyData* vessel_data =
        vtkPolyData::SafeDownCast(vessel_reader->GetOutputDataObject(0));
    assert(vessel_data != nullptr);
    assert(vessel_data->GetNumberOfPoints() == 2);
    assert(vessel_data->GetNumberOfLines() == 1);
    assert(vessel_data->GetNumberOfCells() == 1);
    assert(vessel_data->GetNumberOfVerts() == 0);
    assert(vessel_data->GetNumberOfPolys() == 0);
    assert(vessel_data->GetPointData()->GetArray("node_id")->GetDataTypeSize() == 8);
    assert(vessel_data->GetPointData()->GetArray("vessel_id")->GetDataTypeSize() == 8);
    assert(vessel_data->GetPointData()->GetArray("branch_role")->GetDataTypeSize() == 1);
    assert(vessel_data->GetPointData()->GetArray("perfused")->GetDataTypeSize() == 1);
    assert(vessel_data->GetPointData()->GetArray("diameter_voxels")->GetDataTypeSize() == 4);
    assert(vessel_data->GetPointData()->GetArray("radius_voxels")->GetDataTypeSize() == 4);
    assert(vessel_data->GetPointData()->GetArray("diameter_voxels")->GetTuple1(0) == 2.0);
    assert(vessel_data->GetPointData()->GetArray("radius_voxels")->GetTuple1(0) == 1.0);
    assert(vessel_data->GetPointData()->GetArray("diameter_voxels")->GetTuple1(1) == 3.0);
    assert(vessel_data->GetPointData()->GetArray("radius_voxels")->GetTuple1(1) == 1.5);
    vessel_data->GetPoint(0, point);
    assert(point[0] == 0.5 && point[1] == 0.5 && point[2] == 0.5);
    vessel_data->GetPoint(1, point);
    assert(point[0] == 1.5 && point[1] == 1.5 && point[2] == 0.5);
    std::filesystem::remove(vessel_path);

    // Snapshot output is observational: enabling it must not consume RNG or
    // change the final model state, and a 3D run directory contains no PNG.
    Model3DConfig baseline_config;
    baseline_config.output_enabled = false;
    baseline_config.initial_r_cells = 4;
    baseline_config.initial_K_cells = 4;
    baseline_config.initial_radius = 7;
    baseline_config.end_time_hours = 8.0;
    baseline_config.max_events = 10000;
    Simulation3D baseline(baseline_config);
    baseline.run();

    const auto run_directory = std::filesystem::temp_directory_path() / "atcg3d_vtk_output_test";
    std::filesystem::remove_all(run_directory);
    Model3DConfig output_config = baseline_config;
    output_config.output_enabled = true;
    output_config.output_directory = run_directory;
    output_config.preview_every_hours = 1.0;
    output_config.full_every_hours = 4.0;
    output_config.checkpoint_every_hours = 0.0;
    Simulation3D with_output(output_config);
    OutputManager3D output(output_config);
    with_output.run([&output](const Simulation3D& current) { output.observe(current); });
    output.finalize(with_output);
    assert(with_output.state_checksum() == baseline.state_checksum());
    assert(std::filesystem::exists(run_directory / "preview.vtkhdf.series"));
    assert(std::filesystem::exists(run_directory / "full.vtkhdf.series"));
    assert(std::filesystem::exists(run_directory / "vessels.vtkhdf.series"));
    assert(std::filesystem::is_directory(run_directory / "viz" / "vessels"));
    assert(output.vessel_entries().size() == output.preview_entries().size());
    for (std::size_t index = 0; index < output.vessel_entries().size(); ++index) {
        assert(output.vessel_entries()[index].time_hours ==
               output.preview_entries()[index].time_hours);
    }
    for (const auto& entry : std::filesystem::recursive_directory_iterator(run_directory)) {
        assert(entry.path().extension() != ".png");
        assert(entry.path().extension() != ".tmp");
    }
    std::filesystem::remove_all(run_directory);
}
