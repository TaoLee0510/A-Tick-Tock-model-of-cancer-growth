#include <cassert>
#include <cmath>
#include <filesystem>
#include <memory>
#include <numeric>

#include "engine/simulation.hpp"
#include "model/shared_resource_environment.hpp"
#include "model/structured_pde_model.hpp"
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
#include "io/checkpoint_hdf5.hpp"
#endif

namespace {
atcg3d::structured_pde::StructuredPdeConfig3D config(bool thin) {
    auto result = atcg3d::structured_pde::StructuredPdeConfig3D::load(
        std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_SharedRules/config/angiogenesis_controlled_v15.yaml");
    auto& base = result.continuum.base;
    base.thin_layer = thin;
    base.initial_r_cells = base.initial_K_cells = 1;
    base.initial_large_fraction = 0.5;
    base.initial_radius = 6;
    base.initial_shell_inner_radius = base.initial_inner_small_radius = 4;
    base.initial_shell_thickness = 2;
    base.threads = 1;
    result.continuum.grid.shape = {17, 17, thin ? 1 : 17};
    result.continuum.grid.origin = {-8.0, -8.0, thin ? -0.5 : -8.0};
    result.continuum.reaction.small_daughter_vacancy_exponent = thin ? 8.0 : 26.0;
    result.continuum.end_time_hours = 4.0;
    result.continuum.angiogenesis.seed_tips_per_hour = 100.0;
    result.continuum.angiogenesis.exclusion_fraction = 1.0e-6;
    result.continuum.angiogenesis.tip_diffusion_voxels2_per_hour = 0.5;
    result.migration.vessel_exclusion = true;
    result.validate();
    return result;
}

void check_coupling(bool thin) {
    auto rules = config(thin);
    auto environment = std::make_unique<atcg3d::shared_rules::SharedResourceEnvironment3D>(rules);
    auto* resource = environment.get();
    atcg3d::Simulation3D original(atcg3d::shared_rules::abm_config(rules), std::move(environment));
    original.initialize();
    assert(original.cells().alive_count() == 2);
    bool has_large = false;
    for (const auto slot : original.cells().alive_slots()) {
        has_large = has_large || original.cells().stage(slot) == atcg3d::CellStage::large;
    }
    assert(has_large);
    if (thin) {
        assert(resource->destination_available({0, 0, 0}) == resource->destination_available({0, 0, 1}));
    }
    atcg3d::structured_pde::StructuredPdeModel3D pde(rules);
    pde.initialize_from_abm(original);
    const auto initial = pde.diagnostics();
    assert(initial.r_total + initial.K_total == 2.0);
    while (original.clock().time_hours < 1.0) assert(original.step());
    const auto sidecar = std::filesystem::temp_directory_path() /
        (thin ? "atcg_shared_vascular_2d.resource.bin" : "atcg_shared_vascular_3d.resource.bin");
    std::filesystem::remove(sidecar);
    resource->save_checkpoint(sidecar, original.state_checksum());
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
    const auto hdf5 = sidecar.string() + ".h5";
    std::filesystem::remove(hdf5);
    atcg3d::write_hdf5_checkpoint(hdf5, original);
    const auto saved = atcg3d::read_hdf5_checkpoint(hdf5, atcg3d::shared_rules::abm_config(rules));
    assert(saved.state_checksum == original.state_checksum());
    std::filesystem::remove(hdf5);
#endif
    rules.continuum.base.threads = 8;
    auto reloaded = std::make_unique<atcg3d::shared_rules::SharedResourceEnvironment3D>(rules);
    auto* restarted_resource = reloaded.get();
    restarted_resource->load_checkpoint(sidecar, original.state_checksum());
    atcg3d::Simulation3D restarted(atcg3d::shared_rules::abm_config(rules), std::move(reloaded));
    restarted.restore(original.snapshot_cells(), original.next_uid(), original.clock(),
        original.stats(), original.lineage(), original.snapshot_vasculature(),
        original.cells().slot_count(), original.snapshot_cell_slots(), original.cells().free_slots());
    assert(original.state_checksum() == restarted.state_checksum());
    while (original.step()) {
        assert(restarted.step());
        assert(original.state_checksum() == restarted.state_checksum());
        assert(resource->field_checksum() == restarted_resource->field_checksum());
    }
    assert(!restarted.step());
    const auto removed = std::accumulate(resource->vascular_removed_mass().begin(),
        resource->vascular_removed_mass().end(), 0.0);
    assert(removed > 0.0);
    assert(original.cells().alive_count() + removed == 2.0);
    assert(resource->vascular_removed_mass() == restarted_resource->vascular_removed_mass());
    std::filesystem::remove(sidecar);

    rules.continuum.base.threads = 1;
    while (pde.time_hours() < 1.0) assert(pde.step());
    const auto checkpoint = std::filesystem::temp_directory_path() /
        (thin ? "atcg_shared_vascular_2d.pde.bin" : "atcg_shared_vascular_3d.pde.bin");
    std::filesystem::remove(checkpoint);
    pde.save_checkpoint(checkpoint);
    rules.continuum.base.threads = 8;
    atcg3d::structured_pde::StructuredPdeModel3D resumed(rules);
    resumed.load_checkpoint(checkpoint);
    assert(pde.state_checksum() == resumed.state_checksum());
    while (pde.step()) {
        assert(resumed.step());
        assert(pde.state_checksum() == resumed.state_checksum());
    }
    assert(!resumed.step());
    const auto final = pde.diagnostics();
    assert(final.vascular_removed_mass.total() > 0.0);
    assert(std::abs(final.r_total + final.K_total + final.vascular_removed_mass.total() - 2.0) < 1.0e-12);
    std::filesystem::remove(checkpoint);
}
}  // namespace

int main() {
    check_coupling(true);
    check_coupling(false);
}
