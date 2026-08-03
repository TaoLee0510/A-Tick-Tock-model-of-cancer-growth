#include <cassert>
#include <cmath>
#include <filesystem>
#include <numeric>

#include "config/continuum_config.hpp"
#include "engine/simulation.hpp"
#include "model/continuum_model.hpp"
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
#include "io/checkpoint_hdf5.hpp"
#endif

namespace {

bool close(double lhs, double rhs, double tolerance = 1.0e-9) {
    return std::abs(lhs - rhs) <= tolerance *
        std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

}  // namespace

int main() {
    using namespace atcg3d;
    using namespace atcg3d::continuum;

    const auto wrapper = std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_Continuum/config/continuum_smoke_v1.yaml";
    ContinuumModelConfig3D loaded = ContinuumModelConfig3D::load(wrapper);
    assert(loaded.schema_version == 1);
    assert(loaded.profile == "continuum_smoke_v1");
    assert(loaded.base.profile == "smoke_test_v3");
    assert(loaded.to_json().find("dynamics_fingerprint") != std::string::npos);

    // The exact ABM initialization is coarse-grained conservatively into four
    // population fields. Large footprints are distributed over their occupied
    // sites, so integrated cell number remains 64 rather than occupied volume.
    Model3DConfig base = loaded.base;
    base.output_enabled = false;
    base.control_enabled = false;
    Simulation3D abm(base);
    abm.initialize();
    ContinuumModel3D mapped(loaded);
    mapped.initialize_from_abm(abm);
    const auto mapped_diagnostics = mapped.diagnostics();
    assert(close(mapped_diagnostics.type_mass[0], 32.0));
    assert(close(mapped_diagnostics.type_mass[1], 32.0));
    assert(mapped_diagnostics.maximum_occupied_fraction <= 1.0 + 1.0e-12);
    assert(mapped_diagnostics.maximum_nutrient > 0.0);
    assert(mapped_diagnostics.vessel_volume > 0.0);

#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
    // The same conservative mapping accepts a materialized ATCG3D HDF5 state,
    // which is the production path for importing evolved vessels and cells.
    const auto abm_checkpoint = std::filesystem::temp_directory_path() /
        "atcg3d_continuum_abm_import.h5";
    std::filesystem::remove(abm_checkpoint);
    write_hdf5_checkpoint(abm_checkpoint, abm);
    const CheckpointData3D imported = read_hdf5_checkpoint(abm_checkpoint, base);
    Simulation3D restored_abm(base);
    restored_abm.restore(
        imported.cells, imported.next_uid, imported.clock, imported.stats,
        imported.lineage, imported.vasculature, imported.cell_slot_count,
        imported.cell_slots, imported.cell_free_slots);
    ContinuumModelConfig3D import_config = loaded;
    import_config.initialization_mode = "abm_checkpoint";
    import_config.abm_checkpoint = abm_checkpoint;
    ContinuumModel3D imported_continuum(import_config);
    imported_continuum.initialize_from_abm(restored_abm);
    assert(close(imported_continuum.diagnostics().type_mass[0], 32.0));
    assert(close(imported_continuum.diagnostics().type_mass[1], 32.0));
    std::filesystem::remove(abm_checkpoint);
#endif

    // A single supplied voxel produces a bounded nutrient maximum at the
    // source and a lower value at a distant corner.
    ContinuumModelConfig3D small = loaded;
    small.output.enabled = false;
    small.grid.shape = {5, 5, 5};
    small.grid.origin = {-5.0, -5.0, -5.0};
    small.grid.spacing_voxels = 2.0;
    small.vascular.source_mode = "abm_perfusion";
    small.end_time_hours = 2.0;
    small.validate();
    const std::size_t voxel_count = 125;
    std::array<std::vector<double>, kPopulationFieldCount3D> empty;
    for (auto& field : empty) field.assign(voxel_count, 0.0);
    std::vector<double> vessel(voxel_count, 0.0);
    const std::size_t center = (2U * 5U + 2U) * 5U + 2U;
    vessel[center] = 1.0;
    ContinuumModel3D supplied(small);
    supplied.initialize_from_arrays(empty, vessel);
    assert(supplied.nutrient()[center] > supplied.nutrient().front());
    assert(supplied.nutrient()[center] <= small.nutrient.vessel_value);

    // Binary restart is exact across a nutrient refresh. Both the population
    // PDE fields and the fixed-iteration nutrient state remain bit identical.
    auto populations = empty;
    populations[static_cast<std::size_t>(PopulationField3D::r_small)][center] = 0.1;
    ContinuumModel3D uninterrupted(small);
    uninterrupted.initialize_from_arrays(populations, vessel);
    for (int step = 0; step < 10; ++step) assert(uninterrupted.step());
    const auto checkpoint = std::filesystem::temp_directory_path() /
        "atcg3d_continuum_roundtrip.continuum.bin";
    std::filesystem::remove(checkpoint);
    uninterrupted.save_checkpoint(checkpoint);
    ContinuumModel3D resumed(small);
    resumed.load_checkpoint(checkpoint);
    assert(resumed.state_checksum() == uninterrupted.state_checksum());
    for (int step = 0; step < 20; ++step) {
        assert(uninterrupted.step());
        assert(resumed.step());
    }
    assert(resumed.state_checksum() == uninterrupted.state_checksum());
    std::filesystem::remove(checkpoint);

    // No migration/reaction operation may violate the volume-filling bound.
    const auto final_diagnostics = uninterrupted.diagnostics();
    assert(final_diagnostics.maximum_occupied_fraction <=
           small.reaction.maximum_occupied_fraction + 1.0e-10);
    assert(final_diagnostics.type_mass[0] > 0.0);
    return 0;
}
