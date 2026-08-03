#include <cassert>
#include <cmath>
#include <filesystem>
#include <memory>

#include "config/nutrient_config.hpp"
#include "core/cell_store.hpp"
#include "engine/simulation.hpp"
#include "field/nutrient_field.hpp"
#include "space/domain.hpp"
#include "vasculature/vessel_grid.hpp"
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
#include "io/checkpoint_hdf5.hpp"
#endif

namespace {

atcg3d::SparseVesselGrid3D one_vessel_grid() {
    atcg3d::Model3DConfig base;
    atcg3d::SparseVesselGrid3D vessels(4, atcg3d::DomainPolicy(base));
    const auto placed = vessels.add(
        {0, 0, 0}, atcg3d::VesselBranchRole::root, true, 1);
    assert(placed.placed && placed.newly_occupied_voxels == 1);
    return vessels;
}

atcg3d::nutrient::NutrientFieldConfig3D test_field_config() {
    atcg3d::nutrient::NutrientFieldConfig3D config;
    config.block_edge = 2;
    config.halo_voxels = 6;
    config.diffusion_voxels2_per_hour = 1.0;
    config.decay_per_hour = 0.02;
    config.vessel_exchange_per_hour = 20.0;
    config.vessel_value = 1.0;
    config.r_consumption_per_voxel_hour = 0.2;
    config.K_consumption_per_voxel_hour = 0.2;
    config.solver_iterations = 256;
    config.relaxation = 0.8;
    config.refresh_every_hours = 0.25;
    config.validate();
    return config;
}

}  // namespace

int main() {
    using namespace atcg3d;
    using namespace atcg3d::nutrient;

    const auto wrapper_path = std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_Nutrient/config/nutrient_smoke_v1.yaml";
    const NutrientModelConfig3D loaded = NutrientModelConfig3D::load(wrapper_path);
    assert(loaded.schema_version == 1);
    assert(loaded.profile == "nutrient_smoke_v1");
    assert(loaded.base.profile == "smoke_test_v3");
    assert(loaded.nutrient.maximum_capacity_multiplier == 2.0);
    assert(loaded.to_json().find("effective_resource_surplus_v1") !=
           std::string::npos);

    // A mirrored bundle can place the wrapper and referenced base config in
    // one directory even when the source-tree relative path no longer exists.
    const auto bundle_root = std::filesystem::temp_directory_path() /
        "atcg3d_nutrient_config_bundle_test";
    const auto bundle_directory = bundle_root / "model" / "config";
    const auto bundled_wrapper = bundle_directory / wrapper_path.filename();
    const auto bundled_base = bundle_directory /
        loaded.base_config_path.filename();
    std::filesystem::remove_all(bundle_root);
    std::filesystem::create_directories(bundle_directory);
    std::filesystem::copy_file(wrapper_path, bundled_wrapper);
    std::filesystem::copy_file(loaded.base_config_path, bundled_base);
    const NutrientModelConfig3D bundled =
        NutrientModelConfig3D::load(bundled_wrapper);
    assert(bundled.base_config_path == bundled_base);
    assert(bundled.base.profile == loaded.base.profile);
    std::filesystem::remove_all(bundle_root);

    const NutrientFieldConfig3D field_config = test_field_config();
    SparseVesselGrid3D vessels = one_vessel_grid();
    CellStore3D empty_cells;
    NutrientEnvironment3D supplied(field_config, false);
    const auto initialized = supplied.initialize(0.0, empty_cells, vessels);
    assert(initialized.refresh_cell_rates);
    assert(supplied.value({0, 0, 0}) > supplied.value({3, 0, 0}));
    assert(supplied.value({3, 0, 0}) > 0.0F);
    assert(supplied.value({30, 0, 0}) == 0.0F);
    assert(supplied.capacity_multiplier({0, 0, 0}) > 1.0);
    assert(supplied.capacity_multiplier({0, 0, 0}) <= 2.0);
    assert(supplied.retained_density({0, 0, 0}) < 1.0);
    assert(supplied.diagnostics().perfused_source_voxels == 1);

    CellStore3D consuming_cells;
    CellInit cell;
    cell.uid = 1;
    cell.anchor = {2, 0, 0};
    cell.type = CellType::r;
    cell.stage = CellStage::small;
    consuming_cells.create(cell);
    NutrientEnvironment3D consumed(field_config, false);
    consumed.initialize(0.0, consuming_cells, vessels);
    assert(consumed.value(cell.anchor) < supplied.value(cell.anchor));
    assert(consumed.diagnostics().consuming_voxels == 1);

    // No perfused source means the vascular-surplus field is exactly zero and
    // leaves the legacy carrying capacity unchanged.
    Model3DConfig domain_config;
    SparseVesselGrid3D no_vessels(
        4, DomainPolicy(domain_config));
    NutrientEnvironment3D no_supply(field_config, false);
    no_supply.initialize(0.0, consuming_cells, no_vessels);
    assert(no_supply.value(cell.anchor) == 0.0F);
    assert(no_supply.retained_density(cell.anchor) == 1.0);

    // The binary sidecar preserves the exact continuous field and its refresh
    // schedule independently of the base HDF5 checkpoint.
    const std::filesystem::path sidecar =
        std::filesystem::temp_directory_path() /
        "atcg3d_nutrient_roundtrip.nutrient.bin";
    std::filesystem::remove(sidecar);
    supplied.save_checkpoint(sidecar, 1234, 0.0, 0);
    NutrientEnvironment3D roundtrip(field_config, false);
    roundtrip.load_checkpoint(sidecar, 1234, 0.0, 0);
    const auto restored = roundtrip.initialize(0.0, empty_cells, vessels);
    assert(!restored.refresh_cell_rates);
    assert(roundtrip.field_checksum() == supplied.field_checksum());
    std::filesystem::remove(sidecar);

    // Integration: even with no biological actors, nutrient refreshes are
    // first-class deterministic events in the Simulation3D queue.
    Model3DConfig simulation_config;
    simulation_config.output_enabled = false;
    simulation_config.angiogenesis.enabled = false;
    simulation_config.migration_activation_enabled = false;
    simulation_config.end_time_hours = 0.6;
    simulation_config.max_events = 10;
    auto integration_environment =
        std::make_unique<NutrientEnvironment3D>(field_config, false);
    NutrientEnvironment3D* integration_field = integration_environment.get();
    Simulation3D simulation(
        simulation_config, std::move(integration_environment));
    simulation.restore({}, 1, {}, {}, {});
    assert(simulation.pending_event_count() == 1);
    assert(simulation.step());
    assert(std::abs(simulation.clock().time_hours - 0.25) < 1.0e-12);
    assert(simulation.clock().completed_events == 1);
    assert(integration_field->refresh_count() == 2);

#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
    const std::filesystem::path base_checkpoint =
        std::filesystem::temp_directory_path() /
        "atcg3d_nutrient_integration.h5";
    const std::filesystem::path nutrient_checkpoint =
        std::filesystem::temp_directory_path() /
        "atcg3d_nutrient_integration.nutrient.bin";
    std::filesystem::remove(base_checkpoint);
    std::filesystem::remove(nutrient_checkpoint);
    write_hdf5_checkpoint(base_checkpoint, simulation);
    integration_field->save_checkpoint(
        nutrient_checkpoint, simulation.state_checksum(),
        simulation.clock().time_hours,
        simulation.clock().completed_events);
    const CheckpointData3D checkpoint =
        read_hdf5_checkpoint(base_checkpoint, simulation_config);
    auto resumed_environment =
        std::make_unique<NutrientEnvironment3D>(field_config, false);
    NutrientEnvironment3D* resumed_field = resumed_environment.get();
    resumed_field->load_checkpoint(
        nutrient_checkpoint, checkpoint.state_checksum,
        checkpoint.clock.time_hours, checkpoint.clock.completed_events);
    Simulation3D resumed(simulation_config, std::move(resumed_environment));
    resumed.restore(
        checkpoint.cells, checkpoint.next_uid, checkpoint.clock,
        checkpoint.stats, checkpoint.lineage, checkpoint.vasculature,
        checkpoint.cell_slot_count, checkpoint.cell_slots,
        checkpoint.cell_free_slots);
    assert(resumed.state_checksum() == simulation.state_checksum());
    assert(resumed_field->field_checksum() ==
           integration_field->field_checksum());
    assert(simulation.step());
    assert(resumed.step());
    assert(resumed.state_checksum() == simulation.state_checksum());
    assert(resumed_field->field_checksum() ==
           integration_field->field_checksum());
    std::filesystem::remove(base_checkpoint);
    std::filesystem::remove(nutrient_checkpoint);
#else
    assert(simulation.step());
#endif
    assert(std::abs(simulation.clock().time_hours - 0.5) < 1.0e-12);
    assert(integration_field->refresh_count() == 3);

    return 0;
}
