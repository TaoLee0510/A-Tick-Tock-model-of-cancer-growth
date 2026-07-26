#include <cassert>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include "config/model_config.hpp"
#include "core/cell_store.hpp"
#include "engine/simulation.hpp"
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
#include "io/checkpoint_hdf5.hpp"
#endif
#include "io/output_manager.hpp"
#include "io/preview_sampler.hpp"
#include "io/run_manifest.hpp"
#include "io/vtkhdf_writer.hpp"
#include "vasculature/vessel_store.hpp"

namespace {

void touch_frame(const std::filesystem::path& path) {
    std::filesystem::create_directories(path.parent_path());
    std::ofstream stream(path, std::ios::binary | std::ios::trunc);
    stream << "test";
    assert(stream.good());
}

template <class Action>
bool rejects(Action action) {
    try {
        action();
    } catch (const std::runtime_error&) {
        return true;
    }
    return false;
}

std::string read_all(const std::filesystem::path& path) {
    std::ifstream stream(path, std::ios::binary);
    return std::string((std::istreambuf_iterator<char>(stream)),
                       std::istreambuf_iterator<char>());
}

}  // namespace

int main() {
    using namespace atcg3d;
    CellStore3D cells;
    for (std::uint64_t uid = 1; uid <= 100; ++uid) {
        CellInit cell;
        cell.uid = uid;
        cell.anchor = {static_cast<int>(uid), 0, 0};
        cells.create(cell);
    }
    const auto first = stable_preview_sample(cells, 10, 42);
    const auto second = stable_preview_sample(cells, 10, 42);
    const auto different = stable_preview_sample(cells, 10, 43);
    assert(first == second);
    assert(first.size() == 10);
    assert(first != different);
    assert(stable_preview_sample(cells, 200, 42).size() == 100);

    if (!vtkhdf_output_available()) {
        const auto unavailable_path =
            std::filesystem::temp_directory_path() / "atcg3d_unavailable_vessels.vtkhdf";
        std::filesystem::remove(unavailable_path);
        bool rejected = false;
        try {
            VesselNodeStore3D empty_vessels;
            write_vtkhdf_vessels_atomic(unavailable_path, empty_vessels);
        } catch (const std::runtime_error& error) {
            rejected = std::string(error.what()).find("VTK-HDF vessel output is unavailable") !=
                       std::string::npos;
        }
        assert(rejected);
        assert(!std::filesystem::exists(unavailable_path));
    }

    const auto directory = std::filesystem::temp_directory_path() / "atcg3d_manifest_test";
    std::filesystem::remove_all(directory);
    std::filesystem::create_directories(directory);
    const std::vector<SeriesEntry3D> preview{{"viz/preview/frame_00000000.vtkhdf", 0.0},
                                             {"viz/preview/frame_00000001.vtkhdf", 1.0}};
    const std::vector<SeriesEntry3D> full{{"viz/full/frame_00000000.vtkhdf", 0.0}};
    const std::vector<SeriesEntry3D> vessels{{"viz/vessels/frame_00000000.vtkhdf", 0.0},
                                             {"viz/vessels/frame_00000001.vtkhdf", 1.0}};
    write_series_atomic(directory / "preview.vtkhdf.series", preview);
    write_series_atomic(directory / "full.vtkhdf.series", full);
    write_series_atomic(directory / "vessels.vtkhdf.series", vessels);
    for (const SeriesEntry3D& entry : preview) touch_frame(directory / entry.name);
    for (const SeriesEntry3D& entry : full) touch_frame(directory / entry.name);
    for (const SeriesEntry3D& entry : vessels) touch_frame(directory / entry.name);
    Model3DConfig config;
    write_run_manifest_atomic(directory, config, preview, full, vessels, true);
    assert(std::filesystem::exists(directory / "preview.vtkhdf.series"));
    assert(std::filesystem::exists(directory / "vessels.vtkhdf.series"));
    assert(std::filesystem::exists(directory / "run.json"));
    assert(!std::filesystem::exists(directory / "run.json.tmp"));
    std::ifstream stream(directory / "preview.vtkhdf.series");
    const std::string contents((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());
    assert(contents.find("file-series-version") != std::string::npos);
    assert(contents.find(".tmp") == std::string::npos);
    std::ifstream manifest_stream(directory / "run.json");
    const std::string manifest((std::istreambuf_iterator<char>(manifest_stream)),
                               std::istreambuf_iterator<char>());
    assert(manifest.find("\"schema_version\": 4") != std::string::npos);
    assert(manifest.find("\"vessel_series\": \"vessels.vtkhdf.series\"") !=
           std::string::npos);
    assert(manifest.find("\"vessel_frames\": 2") != std::string::npos);
    assert(manifest.find("\"cell_field_arrays\": [\"total_cell_count\",\"total_slot_count\"]") !=
           std::string::npos);
    assert(manifest.find(
        "\"cell_point_arrays\": [\"cell_id\",\"cell_slot\",\"lesion_id\",\"clone_id\",\"cell_type\",\"stage\",\"viability\",\"display_radius\"]") !=
           std::string::npos);
    assert(manifest.find(
        "\"vessel_point_arrays\": [\"node_id\",\"vessel_id\",\"source_lesion_id\",\"branch_role\",\"perfused\",\"diameter_voxels\",\"radius_voxels\"]") !=
           std::string::npos);
    assert(manifest.find("\"dynamics_config_json\"") != std::string::npos);

    const ExistingRunOutput3D loaded = load_existing_run_output(directory, config);
    assert(loaded.preview.size() == 2);
    assert(loaded.full.size() == 1);
    assert(loaded.vessels.size() == 2);

    Model3DConfig mismatched = config;
    mismatched.alpha += 0.25;
    assert(rejects([&] { (void)load_existing_run_output(directory, mismatched); }));

    {
        std::string incompatible_manifest = manifest;
        const std::string required_array = "\"lesion_id\",";
        const std::size_t position = incompatible_manifest.find(required_array);
        assert(position != std::string::npos);
        incompatible_manifest.erase(position, required_array.size());
        std::ofstream output(directory / "run.json",
                             std::ios::binary | std::ios::trunc);
        output << incompatible_manifest;
        assert(output.good());
    }
    assert(rejects([&] { (void)load_existing_run_output(directory, config); }));
    {
        std::ofstream output(directory / "run.json",
                             std::ios::binary | std::ios::trunc);
        output << manifest;
        assert(output.good());
    }

    {
        std::ofstream malformed_series(directory / "preview.vtkhdf.series",
                                       std::ios::binary | std::ios::trunc);
        malformed_series
            << "{\"file-series-version\":\"1.0\",\"files\":["
            << "{\"name\":\"viz/preview/frame_00000000.vtkhdf\",\"time\":1},"
            << "{\"name\":\"viz/preview/frame_00000001.vtkhdf.tmp\",\"time\":0}]}";
    }
    assert(rejects([&] { (void)load_existing_run_output(directory, config); }));
    std::filesystem::remove_all(directory);

    // OutputManager resume loads existing catalogs, rejects fresh overwrite,
    // and aligns lineage with the checkpoint-restored in-memory prefix before
    // it appends anything. This test does not require VTK because all periodic
    // output intervals are disabled.
    const auto resume_directory =
        std::filesystem::temp_directory_path() / "atcg3d_output_resume_test";
    std::filesystem::remove_all(resume_directory);
    Model3DConfig fresh;
    fresh.output_enabled = true;
    fresh.output_directory = resume_directory;
    fresh.preview_every_hours = 0.0;
    fresh.full_every_hours = 0.0;
    fresh.checkpoint_every_hours = 0.0;
    fresh.end_time_hours = 0.0;
    fresh.max_events = 1;
    fresh.initial_r_cells = 1;
    fresh.initial_K_cells = 1;
    fresh.initial_large_fraction = 0.0;
    Simulation3D simulation(fresh);
    CellInit restored_cell;
    restored_cell.uid = 1;
    restored_cell.clone_id = 1;
    restored_cell.anchor = {0, 0, 0};
    LineageEdge restored_edge;
    restored_edge.child_uid = 1;
    restored_edge.clone_id = 1;
    restored_edge.type = CellType::r;
    simulation.restore({restored_cell}, 2, {}, {}, {restored_edge});
    {
        OutputManager3D output(fresh);
        output.observe(simulation);
    }
    const std::filesystem::path lineage_path =
        resume_directory / "lineage" / "edges.csv";
    const std::string original_lineage = read_all(lineage_path);
    assert(!original_lineage.empty());
    assert(rejects([&] { OutputManager3D duplicate(fresh); }));

    const auto unrelated_directory =
        std::filesystem::temp_directory_path() / "atcg3d_output_nonempty_test";
    std::filesystem::remove_all(unrelated_directory);
    std::filesystem::create_directories(unrelated_directory);
    touch_frame(unrelated_directory / "unrelated.file");
    Model3DConfig nonempty = fresh;
    nonempty.output_directory = unrelated_directory;
    assert(rejects([&] { OutputManager3D overwrite(nonempty); }));
    std::filesystem::remove_all(unrelated_directory);

    Model3DConfig resume = fresh;
    resume.run_mode = "resume";
    resume.resume_checkpoint = "checkpoint.h5";
    {
        OutputManager3D output(resume);
        assert(output.preview_entries().empty());
        output.observe(simulation);
    }
    assert(read_all(lineage_path) == original_lineage);

    {
        std::ofstream lineage(lineage_path, std::ios::binary | std::ios::app);
        lineage << "1,2,1,1,1\n";
        assert(lineage.good());
    }
    {
        OutputManager3D lineage_recovery(resume);
        lineage_recovery.observe(simulation);
    }
    const std::filesystem::path lineage_recovery_directory =
        resume_directory / "recovery" / "checkpoint_0000000000000000";
    assert(read_all(lineage_path) == original_lineage);
    assert(read_all(lineage_recovery_directory / "lineage" / "edges_tail.csv").find(
               "1,2,1,1,1") != std::string::npos);
    assert(std::filesystem::is_regular_file(lineage_recovery_directory /
                                            "recovery.json"));

    {
        std::ofstream lineage(lineage_path, std::ios::binary | std::ios::trunc);
        lineage << "birth_time,child_uid,parent_uid,clone_id,type\n"
                << "0,999999,0,1,1\n";
    }
    OutputManager3D mismatched_lineage(resume);
    assert(rejects([&] { mismatched_lineage.observe(simulation); }));
    std::filesystem::remove_all(resume_directory);

    if (vtkhdf_output_available()) {
        const auto vtk_resume_directory =
            std::filesystem::temp_directory_path() / "atcg3d_vtk_output_resume_test";
        std::filesystem::remove_all(vtk_resume_directory);
        Model3DConfig vtk_fresh = fresh;
        vtk_fresh.output_directory = vtk_resume_directory;
        vtk_fresh.preview_every_hours = 1.0;
        vtk_fresh.full_every_hours = 2.0;

        Simulation3D at_zero(vtk_fresh);
        at_zero.restore({restored_cell}, 2, {0, 0.0}, {}, {restored_edge});
        {
            OutputManager3D output(vtk_fresh);
            output.observe(at_zero);
            assert(output.preview_entries().size() == 1);
            assert(output.full_entries().size() == 1);
            assert(output.vessel_entries().size() == 1);
        }
        const std::filesystem::path original_frame =
            vtk_resume_directory / "viz" / "preview" / "frame_00000000.vtkhdf";
        const std::uintmax_t original_size = std::filesystem::file_size(original_frame);

        Model3DConfig vtk_resume = vtk_fresh;
        vtk_resume.run_mode = "resume";
        vtk_resume.resume_checkpoint = "checkpoint.h5";
        Simulation3D at_checkpoint(vtk_resume);
        at_checkpoint.restore({restored_cell}, 2, {10, 1.0}, {}, {restored_edge});
        {
            OutputManager3D no_progress(vtk_resume);
            no_progress.observe(at_checkpoint);
            no_progress.finalize(at_checkpoint);
            assert(no_progress.preview_entries().size() == 1);
        }
        OutputManager3D output(vtk_resume);
        output.observe(at_checkpoint);
        assert(output.preview_entries().size() == 1);

        Simulation3D after_resume(vtk_resume);
        after_resume.restore({restored_cell}, 2, {20, 2.0}, {}, {restored_edge});
        output.observe(after_resume);
        assert(output.preview_entries().size() == 2);
        assert(output.full_entries().size() == 2);
        assert(output.vessel_entries().size() == 2);
        assert(output.preview_entries().back().name ==
               "viz/preview/frame_00000001.vtkhdf");
        assert(output.preview_entries().back().time_hours == 2.0);
        assert(std::filesystem::file_size(original_frame) == original_size);
        assert(!std::filesystem::exists(original_frame.string() + ".tmp"));
        std::filesystem::remove_all(vtk_resume_directory);

        // Simulate a process that wrote visualization frames, lineage, final
        // metrics, and a later checkpoint after the checkpoint selected for
        // resume. Resume must publish the selected checkpoint prefix, preserve
        // every superseded artifact under recovery/, and safely reuse the next
        // canonical frame indices.
        const auto crash_directory =
            std::filesystem::temp_directory_path() / "atcg3d_crash_ahead_resume_test";
        std::filesystem::remove_all(crash_directory);
        Model3DConfig crash_fresh = vtk_fresh;
        crash_fresh.output_directory = crash_directory;
        // Exercise resume recovery with the production asynchronous path. The
        // output catalogs are loaded before validation trims artifacts newer
        // than the selected checkpoint, so the scheduled-time de-duplication
        // markers must be rolled back with the catalogs.
        crash_fresh.output_async_enabled = true;
        crash_fresh.output_async_queue_depth = 1;

        LineageEdge future_edge;
        future_edge.birth_time = 2.0;
        future_edge.child_uid = 2;
        future_edge.parent_uid = 1;
        future_edge.clone_id = 1;
        future_edge.type = CellType::r;

        Simulation3D crash_at_zero(crash_fresh);
        crash_at_zero.restore({restored_cell}, 2, {0, 0.0}, {}, {restored_edge});
        Simulation3D crash_at_checkpoint(crash_fresh);
        crash_at_checkpoint.restore({restored_cell}, 2, {10, 1.0}, {},
                                    {restored_edge});
        Simulation3D crash_ahead(crash_fresh);
        crash_ahead.restore({restored_cell}, 3, {20, 2.0}, {},
                            {restored_edge, future_edge});
        const std::filesystem::path selected_checkpoint =
            crash_directory / "checkpoints" /
            "checkpoint_0000000000000010.h5";
        {
            OutputManager3D output(crash_fresh);
            output.observe(crash_at_zero);
            output.observe(crash_at_checkpoint);
            output.observe(crash_ahead);
            output.finalize(crash_ahead);
            assert(output.preview_entries().size() == 3);
            assert(output.full_entries().size() == 2);
            assert(output.vessel_entries().size() == 3);
        }
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
        // Homebrew HDF5 is not built thread-safe. Production serializes VTK-HDF
        // and checkpoint writes through the same output worker, so construct
        // this synthetic selected checkpoint only after that worker has joined.
        write_hdf5_checkpoint(selected_checkpoint, crash_at_checkpoint);
#endif

        touch_frame(crash_directory / "viz" / "preview" /
                    "frame_00000003.vtkhdf");
        touch_frame(crash_directory / "viz" / "full" /
                    "frame_00000001.vtkhdf.tmp");
        touch_frame(crash_directory / "checkpoints" /
                    "checkpoint_0000000000000030.h5");

        Model3DConfig crash_resume = crash_fresh;
        crash_resume.run_mode = "resume";
        crash_resume.resume_checkpoint = selected_checkpoint;
        Simulation3D restored_checkpoint(crash_resume);
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
        const CheckpointData3D selected =
            read_hdf5_checkpoint(selected_checkpoint, crash_resume);
        restored_checkpoint.restore(
            selected.cells, selected.next_uid, selected.clock, selected.stats,
            selected.lineage, selected.vasculature, selected.cell_slot_count,
            selected.cell_slots, selected.cell_free_slots);
#else
        restored_checkpoint.restore({restored_cell}, 2, {10, 1.0}, {},
                                    {restored_edge});
#endif
        OutputManager3D recovered(crash_resume);
        recovered.observe(restored_checkpoint);
        assert(recovered.preview_entries().size() == 2);
        assert(recovered.full_entries().size() == 1);
        assert(recovered.vessel_entries().size() == 2);

        const std::filesystem::path recovery =
            crash_directory / "recovery" / "checkpoint_0000000000000010";
        assert(std::filesystem::is_regular_file(recovery / "recovery.json"));
        assert(std::filesystem::is_regular_file(
            recovery / "viz" / "preview" / "frame_00000002.vtkhdf"));
        assert(std::filesystem::is_regular_file(
            recovery / "viz" / "preview" / "frame_00000003.vtkhdf"));
        assert(std::filesystem::is_regular_file(
            recovery / "viz" / "full" / "frame_00000001.vtkhdf"));
        assert(std::filesystem::is_regular_file(
            recovery / "viz" / "full" / "frame_00000001.vtkhdf.tmp"));
        assert(std::filesystem::is_regular_file(
            recovery / "viz" / "vessels" / "frame_00000002.vtkhdf"));
        assert(std::filesystem::is_regular_file(
            recovery / "checkpoints" / "checkpoint_0000000000000030.h5"));
        assert(std::filesystem::is_regular_file(
            recovery / "lineage" / "edges_tail.csv"));
        assert(std::filesystem::is_regular_file(
            recovery / "metrics" / "final.json"));
        assert(read_all(recovery / "lineage" / "edges_tail.csv").find(
                   "2,2,1,1,1") != std::string::npos);
        assert(read_all(crash_directory / "lineage" / "edges.csv").find(
                   "2,2,1,1,1") == std::string::npos);
        assert(!std::filesystem::exists(
            crash_directory / "viz" / "preview" / "frame_00000002.vtkhdf"));
        assert(!std::filesystem::exists(crash_directory / "metrics" / "final.json"));

        const ExistingRunOutput3D recovered_catalogs =
            load_existing_run_output(crash_directory, crash_resume);
        assert(recovered_catalogs.preview.size() == 2);
        assert(recovered_catalogs.full.size() == 1);
        assert(recovered_catalogs.vessels.size() == 2);

        Simulation3D continued(crash_resume);
        continued.restore({restored_cell}, 3, {20, 2.0}, {},
                          {restored_edge, future_edge});
        recovered.observe(continued);
        recovered.finalize(continued);
        assert(recovered.preview_entries().size() == 3);
        assert(recovered.full_entries().size() == 2);
        assert(recovered.vessel_entries().size() == 3);
        assert(recovered.preview_entries().back().name ==
               "viz/preview/frame_00000002.vtkhdf");
        assert(recovered.full_entries().back().name ==
               "viz/full/frame_00000001.vtkhdf");
        assert(std::filesystem::is_regular_file(
            crash_directory / "viz" / "preview" / "frame_00000002.vtkhdf"));
        assert(std::filesystem::is_regular_file(
            recovery / "viz" / "preview" / "frame_00000002.vtkhdf"));
        std::filesystem::remove_all(crash_directory);
    }
}
