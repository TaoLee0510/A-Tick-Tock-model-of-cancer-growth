#include <algorithm>
#include <cassert>
#include <filesystem>
#include <fstream>
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
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
#include <H5Cpp.h>
#include "io/checkpoint_hdf5.hpp"
#endif
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
    CellInit supporting;
    supporting.uid = 44;
    supporting.anchor = {-1, 3, 4};
    const Slot third = cells.create(supporting);
    const std::vector<Slot> slots{first, second};
    LesionIndexConfig3D lesion_config;
    lesion_config.block_edge = 8;
    lesion_config.core_activation_occupied_fraction = 0.0;
    lesion_config.core_deactivation_occupied_fraction = 0.0;
    lesion_config.minimum_cells_per_core_block = 2;
    lesion_config.halo_blocks = 0;
    LesionIndex3D lesions(lesion_config);
    for (const Slot slot : {first, second, third}) {
        lesions.add_cell_anchor(cells.anchor(slot));
        lesions.add_occupied_site(cells.anchor(slot));
    }
    lesions.refresh_topology();
    const LesionId first_lesion = lesions.lesion_for_anchor(large.anchor)
        .value_or(kNoLesionId);
    assert(first_lesion != kNoLesionId);
    assert(lesions.lesion_for_anchor(supporting.anchor) == first_lesion);
    assert(!lesions.lesion_for_anchor(small.anchor).has_value());
    DisplayRadiusConfig radii;
    radii.large = 1.75;
    radii.small = 0.625;
    radii.ultrasmall = 0.2;
    const SimulationSnapshotView3D snapshot{
        cells, slots, lesions, {3, 1.5}, radii};

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
    assert(data->GetPointData()->GetArray("cell_slot")->GetDataTypeSize() == 4);
    assert(data->GetPointData()->GetArray("cell_slot")->GetTuple1(0) == first);
    assert(data->GetPointData()->GetArray("cell_slot")->GetTuple1(1) == second);
    assert(data->GetPointData()->GetArray("lesion_id")->GetDataTypeSize() == 8);
    assert(data->GetPointData()->GetArray("lesion_id")->GetTuple1(0) ==
           static_cast<double>(first_lesion));
    assert(data->GetPointData()->GetArray("lesion_id")->GetTuple1(1) == 0.0);
    assert(data->GetPointData()->GetArray("clone_id")->GetDataTypeSize() == 4);
    assert(data->GetPointData()->GetArray("cell_type")->GetDataTypeSize() == 1);
    assert(data->GetPointData()->GetArray("stage")->GetDataTypeSize() == 1);
    assert(data->GetPointData()->GetArray("viability")->GetDataTypeSize() == 1);
    assert(data->GetPointData()->GetArray("display_radius")->GetDataTypeSize() == 4);
    assert(data->GetPointData()->GetArray("display_radius")->GetTuple1(0) == 1.75);
    assert(data->GetPointData()->GetArray("display_radius")->GetTuple1(1) == 0.625);
    assert(data->GetFieldData()->GetArray("total_cell_count") != nullptr);
    assert(data->GetFieldData()->GetArray("total_cell_count")->GetDataTypeSize() == 8);
    assert(data->GetFieldData()->GetArray("total_cell_count")->GetTuple1(0) == 3.0);
    assert(data->GetFieldData()->GetArray("total_slot_count") != nullptr);
    assert(data->GetFieldData()->GetArray("total_slot_count")->GetTuple1(0) == 3.0);
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
    root.source_lesion_id = 19;
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
    child.source_lesion_id = root.source_lesion_id;
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
    assert(vessel_data->GetPointData()->GetArray("source_lesion_id")->GetDataTypeSize() == 8);
    assert(vessel_data->GetPointData()->GetArray("source_lesion_id")->GetTuple1(0) == 19.0);
    assert(vessel_data->GetPointData()->GetArray("source_lesion_id")->GetTuple1(1) == 19.0);
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
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
    output_config.checkpoint_every_hours = 4.0;
#else
    output_config.checkpoint_every_hours = 0.0;
#endif
    output_config.output_async_enabled = true;
    output_config.output_async_queue_depth = 1;
    output_config.live_preview_when_attached = true;
    output_config.live_preview_wall_interval_seconds = 0.001;
    Simulation3D with_output(output_config);
    OutputManager3D output(output_config);
    std::filesystem::create_directories(run_directory / "control");
    {
        std::ofstream marker(
            run_directory / "control" / "viewer.attached",
            std::ios::binary | std::ios::trunc);
        marker << "1\n";
    }
    with_output.run([&output](const Simulation3D& current) { output.observe(current); });
    output.finalize(with_output);
    assert(with_output.state_checksum() == baseline.state_checksum());
    assert(std::filesystem::exists(run_directory / "preview.vtkhdf.series"));
    assert(std::filesystem::exists(run_directory / "full.vtkhdf.series"));
    assert(std::filesystem::exists(run_directory / "vessels.vtkhdf.series"));
    assert(std::filesystem::is_directory(run_directory / "viz" / "vessels"));
    assert(std::filesystem::is_regular_file(
        run_directory / "viz" / "live" / "current.vtkhdf"));
    assert(std::filesystem::is_regular_file(
        run_directory / "viz" / "live" / "vessels.vtkhdf"));
    assert(std::filesystem::is_regular_file(
        run_directory / "live.vtkhdf.series"));
    assert(std::filesystem::is_regular_file(
        run_directory / "live-vessels.vtkhdf.series"));
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
    assert(std::filesystem::is_directory(run_directory / "checkpoints"));
    assert(!std::filesystem::is_empty(run_directory / "checkpoints"));
    std::vector<std::filesystem::path> evolving_checkpoints;
    for (const auto& entry : std::filesystem::directory_iterator(
             run_directory / "checkpoints")) {
        if (entry.path().extension() == ".h5") {
            evolving_checkpoints.push_back(entry.path());
        }
    }
    std::sort(evolving_checkpoints.begin(), evolving_checkpoints.end());
    assert(evolving_checkpoints.size() >= 3);
    {
        H5::H5File second(
            evolving_checkpoints[1].string(), H5F_ACC_RDONLY);
        std::uint32_t schema{};
        second.openGroup("/meta")
            .openAttribute("schema_version")
            .read(H5::PredType::NATIVE_UINT32, &schema);
        assert(schema == kCheckpointJournalDeltaSchemaVersion3D);
    }
    const CheckpointData3D evolving_reconstructed =
        read_hdf5_checkpoint(evolving_checkpoints.back(), output_config);
    assert(evolving_reconstructed.state_checksum ==
           with_output.state_checksum());
#endif
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

#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
    // Incremental mode keeps hourly logical state while bounding self-contained
    // full frames. The second checkpoint changes only global clock state, so it
    // must be a schema-v7 field delta with zero updated cell rows and
    // reconstruct exactly.
    const auto incremental_directory =
        std::filesystem::temp_directory_path() /
        "atcg3d_incremental_output_test";
    std::filesystem::remove_all(incremental_directory);
    Model3DConfig incremental_config;
    incremental_config.output_enabled = true;
    incremental_config.output_directory = incremental_directory;
    incremental_config.storage_mode = "journal_delta_hdf5_v2";
    incremental_config.checkpoint_format =
        "hdf5_base_v6_slot_journal_v8";
    incremental_config.preview_every_hours = 1.0;
    incremental_config.full_every_hours = 1.0;
    incremental_config.checkpoint_every_hours = 1.0;
    incremental_config.preview_keyframe_every_hours = 1.0;
    incremental_config.full_keyframe_every_hours = 24.0;
    incremental_config.checkpoint_base_every_hours = 168.0;
    incremental_config.checkpoint_max_delta_chain = 168;
    incremental_config.delta_full_ratio = 0.70;
    incremental_config.output_async_enabled = false;
    incremental_config.initial_r_cells = 0;
    incremental_config.initial_K_cells = 0;
    incremental_config.validate();

    CellInit quiet_cell;
    quiet_cell.uid = 1;
    quiet_cell.stage = CellStage::small;
    quiet_cell.density_growth_rate = 1.0F;
    quiet_cell.division_work_remaining = 100.0F;
    quiet_cell.next_division_time = 100.0;
    quiet_cell.migration_schedule_generation = 1;
    quiet_cell.division_schedule_generation = 1;
    quiet_cell.death_schedule_generation = 1;
    Simulation3D at_zero(incremental_config);
    at_zero.restore({quiet_cell}, 2, {}, {}, {});
    SimulationClock3D one_hour;
    one_hour.time_hours = 1.0;
    one_hour.completed_events = 1;
    Simulation3D at_one(incremental_config);
    at_one.restore({quiet_cell}, 2, one_hour, {}, {});
    OutputManager3D incremental_output(incremental_config);
    incremental_output.observe(at_zero);
    incremental_output.observe(at_one);
    incremental_output.finalize(at_one);
    assert(incremental_output.preview_entries().size() == 2);
    // full_every_hours is the actual full-frame sampling interval. Storage no
    // longer silently suppresses logical frames behind a second keyframe
    // interval.
    assert(incremental_output.full_entries().size() == 2);

    std::vector<std::filesystem::path> checkpoint_paths;
    for (const auto& entry : std::filesystem::directory_iterator(
             incremental_directory / "checkpoints")) {
        if (entry.path().extension() == ".h5") {
            checkpoint_paths.push_back(entry.path());
        }
    }
    std::sort(checkpoint_paths.begin(), checkpoint_paths.end());
    assert(checkpoint_paths.size() == 2);
    std::uint32_t schemas[2]{};
    for (std::size_t index = 0; index < checkpoint_paths.size(); ++index) {
        H5::H5File file(checkpoint_paths[index].string(), H5F_ACC_RDONLY);
        file.openGroup("/meta").openAttribute("schema_version").read(
            H5::PredType::NATIVE_UINT32, &schemas[index]);
    }
    assert(schemas[0] == kCheckpointBaseSchemaVersion3D &&
           schemas[1] == kCheckpointJournalDeltaSchemaVersion3D);
    const CheckpointData3D reconstructed =
        read_hdf5_checkpoint(checkpoint_paths.back(), incremental_config);
    assert(reconstructed.state_checksum == at_one.state_checksum());
    {
        H5::H5File delta(checkpoint_paths.back().string(), H5F_ACC_RDONLY);
        assert(delta.openDataSet("/cells/changed/uid")
                   .getSpace().getSimpleExtentNpoints() == 0);
    }
    std::filesystem::remove_all(incremental_directory);
#endif
}
