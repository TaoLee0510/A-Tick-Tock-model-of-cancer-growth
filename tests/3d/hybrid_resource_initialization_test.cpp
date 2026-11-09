#include <algorithm>
#include <cassert>
#include <filesystem>
#include <memory>

#include "engine/simulation.hpp"
#include "model/hybrid_model.hpp"
#include "model/shared_resource_environment.hpp"

using namespace atcg3d;

namespace {
void compare_initial_resources(bool thin) {
    auto config = hybrid::HybridConfig3D::load(std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_Hybrid/config/production_2d_2000_r200_v5.yaml");
    auto& rules = config.rules;
    auto& base = rules.continuum.base;
    base.thin_layer = thin;
    base.initial_r_cells = 32;
    base.initial_K_cells = 32;
    base.initial_radius = 8;
    base.initial_shell_thickness = 2;
    base.initial_shell_inner_radius = 6;
    base.initial_inner_small_radius = 6;
    base.migration_activation_threshold = thin ? 0.001 : 1.0e-5;
    base.activated_r_normal_multiplier = 20.0;
    base.threads = 1;
    rules.migration.reactivation_density_threshold = thin ? 0.0005 : 5.0e-6;
    rules.continuum.reaction.small_daughter_vacancy_exponent = thin ? 8.0 : 26.0;
    rules.continuum.migration.activated_r_mobility_multiplier = 20.0;
    rules.continuum.grid.shape = {32, 32, thin ? 1 : 32};
    rules.continuum.grid.origin = {-16, -16, thin ? -0.5 : -16};
    rules.continuum.nutrient.boundary_mode = thin
        ? "moving_tumor_front_and_vessels_dirichlet_v2" : "planar_edges_dirichlet_v1";
    rules.continuum.vascular.source_mode = "static_voxels";
    rules.continuum.vascular.static_sources.clear();
    config.validate();
    auto environment = std::make_unique<shared_rules::SharedResourceEnvironment3D>(rules);
    auto* resource = environment.get();
    Simulation3D reference(shared_rules::abm_config(rules), std::move(environment));
    reference.initialize();
    hybrid::HybridModel3D mixed(config);
    mixed.initialize();
    assert(reference.cells().alive_count() == 64);
    assert(mixed.diagnostics().total_mass == 64.0);
    assert(resource->nutrient() == mixed.pde().nutrient());
    assert(resource->vessel_fraction() == mixed.pde().vessel_fraction());
    assert(*std::min_element(resource->nutrient().begin(), resource->nutrient().end()) < 0.2);
    std::size_t active = 0;
    for (const auto slot : reference.cells().alive_slots()) {
        if (reference.cells().type(slot) != CellType::r) {
            continue;
        }
        const auto uid = reference.cells().uid(slot);
        const auto& agents = mixed.abm().cells();
        Slot other = kEmptySlot;
        for (const auto candidate : agents.alive_slots()) {
            if (agents.uid(candidate) == uid) {
                other = candidate;
                break;
            }
        }
        assert(other != kEmptySlot);
        assert(reference.cells().anchor(slot) == agents.anchor(other));
        assert(reference.cells().division_work_remaining(slot) == agents.division_work_remaining(other));
        assert(reference.cells().migration_activation_end_time(slot) == agents.migration_activation_end_time(other));
        active += bool(reference.cells().flags(slot) & kMigrationActive);
    }
    assert(active > 0);
}
}  // namespace

int main() {
    compare_initial_resources(true);
    compare_initial_resources(false);
}
