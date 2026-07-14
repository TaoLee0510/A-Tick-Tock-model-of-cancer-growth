#include <cassert>
#include <filesystem>
#include <vector>

#include <vtkDataArray.h>
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
    const SimulationSnapshotView3D snapshot{cells, slots, {3, 1.5}};

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
    std::filesystem::remove(path);

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
    for (const auto& entry : std::filesystem::recursive_directory_iterator(run_directory)) {
        assert(entry.path().extension() != ".png");
        assert(entry.path().extension() != ".tmp");
    }
    std::filesystem::remove_all(run_directory);
}
