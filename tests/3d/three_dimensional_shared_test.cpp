#include <cassert>
#include <filesystem>
#include <memory>

#include "engine/simulation.hpp"
#include "model/shared_angiogenesis.hpp"
#include "model/shared_resource_environment.hpp"
#include "model/structured_pde_model.hpp"

namespace {
using namespace atcg3d;

structured_pde::StructuredPdeConfig3D config() {
    auto c = structured_pde::StructuredPdeConfig3D::load(
        std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_SharedRules/config/three_dimensional_r20_v16.yaml");
    auto& base = c.continuum.base;
    base.initial_r_cells = base.initial_K_cells = 32;
    base.initial_radius = 6;
    base.initial_shell_inner_radius = base.initial_inner_small_radius = 4;
    base.initial_shell_thickness = 2;
    base.migration_activation_threshold = 1e-5;
    c.migration.reactivation_density_threshold = 5e-6;
    c.continuum.grid.shape = {17, 17, 17};
    c.continuum.grid.origin = {-8.0, -8.0, -8.0};
    c.continuum.vascular.source_mode = "static_voxels";
    c.continuum.vascular.static_sources.clear();
    c.continuum.base.threads = 1;
    c.validate();
    return c;
}

void equal_fields(const structured_pde::StructuredPdeModel3D& a,
                  const structured_pde::StructuredPdeModel3D& b) {
    assert(a.time_hours() == b.time_hours());
    assert(a.nutrient() == b.nutrient());
    assert(a.vessel_fraction() == b.vessel_fraction());
    for (std::size_t i = 0; i < a.voxel_count(); ++i) {
        for (const auto stage : {structured_pde::StructuredStage3D::small,
                                structured_pde::StructuredStage3D::large}) {
            assert(a.r_normal(stage, i) == b.r_normal(stage, i));
            assert(a.r_active(stage, i) == b.r_active(stage, i));
            assert(a.K(stage, i) == b.K(stage, i));
            assert(a.activation_density(stage, i) == b.activation_density(stage, i));
        }
    }
    assert(a.angiogenesis()->checksum() == b.angiogenesis()->checksum());
}
}  // namespace

int main() {
    auto rules = config();
    auto environment = std::make_unique<shared_rules::SharedResourceEnvironment3D>(rules);
    auto* resource = environment.get();
    Simulation3D original(shared_rules::abm_config(rules), std::move(environment));
    original.initialize();
    assert(original.cells().alive_count() == 64);
    assert(resource->prepared_sector_tiles() > 0);
    std::size_t active = 0;
    for (const auto slot : original.cells().alive_slots())
        active += bool(original.cells().flags(slot) & kMigrationActive);
    assert(active > 0);

    structured_pde::StructuredPdeModel3D sparse(rules);
    sparse.initialize_from_abm(original);
    auto dense_rules = rules;
    dense_rules.storage_model = "dense_v1";
    structured_pde::StructuredPdeModel3D dense(dense_rules);
    dense.initialize_from_abm(original);
    equal_fields(sparse, dense);
    while (sparse.time_hours() < 0.5) {
        assert(sparse.step());
        assert(dense.step());
        equal_fields(sparse, dense);
    }
    const auto pde_checkpoint = std::filesystem::current_path() / "shared-3d-pde.bin";
    std::filesystem::remove(pde_checkpoint);
    sparse.save_checkpoint(pde_checkpoint);
    rules.continuum.base.threads = 8;
    structured_pde::StructuredPdeModel3D resumed(rules);
    resumed.load_checkpoint(pde_checkpoint);
    assert(sparse.state_checksum() == resumed.state_checksum());
    while (sparse.step()) {
        assert(dense.step());
        assert(resumed.step());
        equal_fields(sparse, dense);
        assert(sparse.state_checksum() == resumed.state_checksum());
    }
    assert(sparse.angiogenesis()->shared_diagnostics()->seeded_tips > 0.0);
    assert(sparse.angiogenesis()->shared_diagnostics()->centerline_growth > 0.0);
    std::filesystem::remove(pde_checkpoint);

    while (original.clock().time_hours < 0.5) assert(original.step());
    const auto sidecar = std::filesystem::current_path() / "shared-3d-resource.bin";
    std::filesystem::remove(sidecar);
    resource->save_checkpoint(sidecar, original.state_checksum());
    auto reload = std::make_unique<shared_rules::SharedResourceEnvironment3D>(rules);
    auto* restored_resource = reload.get();
    restored_resource->load_checkpoint(sidecar, original.state_checksum());
    Simulation3D restarted(shared_rules::abm_config(rules), std::move(reload));
    restarted.restore(original.snapshot_cells(), original.next_uid(), original.clock(),
        original.stats(), original.lineage(), original.snapshot_vasculature(),
        original.cells().slot_count(), original.snapshot_cell_slots(), original.cells().free_slots());
    assert(original.state_checksum() == restarted.state_checksum());
    assert(restored_resource->prepared_sector_tiles() > 0);
    while (original.step()) {
        assert(restarted.step());
        assert(original.state_checksum() == restarted.state_checksum());
        assert(resource->field_checksum() == restored_resource->field_checksum());
    }
    assert(original.stats().migration_attempts > 0);
    assert(resource->angiogenesis()->shared_diagnostics()->seeded_tips > 0);
    assert(resource->angiogenesis()->shared_diagnostics()->centerline_growth > 0);
    std::filesystem::remove(sidecar);

    rules.schema_version = 15;
    bool rejected = false;
    try { rules.validate(); } catch (const std::invalid_argument&) { rejected = true; }
    assert(rejected);
}
