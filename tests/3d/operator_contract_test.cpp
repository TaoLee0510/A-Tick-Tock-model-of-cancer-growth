#include <algorithm>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <memory>
#include <numeric>

#include "engine/simulation.hpp"
#include "model/feasible_jump.hpp"
#include "model/shared_resource_environment.hpp"
#include "model/structured_pde_model.hpp"
#include "rules/initial_rates.hpp"

namespace {

void feasible_choice() {
    using atcg3d::structured_pde::uniform_feasible_jump_probabilities;
    const std::vector<double> availability{0.0, 0.2, 0.8, 1.0, 0.1, 0.7, 0.4, 0.9};
    const auto probabilities = uniform_feasible_jump_probabilities(availability);
    std::vector<double> enumerated(8, 0.0);
    for (unsigned subset = 1; subset < 256; ++subset) {
        double chance = 1.0;
        unsigned count = 0;
        for (unsigned direction = 0; direction < 8; ++direction) {
            const bool selected = subset & (1U << direction);
            chance *= selected ? availability[direction] : 1.0 - availability[direction];
            count += selected;
        }
        for (unsigned direction = 0; direction < 8; ++direction) {
            if (subset & (1U << direction)) enumerated[direction] += chance / count;
        }
    }
    for (std::size_t direction = 0; direction < 8; ++direction) {
        assert(std::abs(probabilities[direction] - enumerated[direction]) < 1.0e-14);
    }
    for (const std::size_t count : {8U, 26U}) {
        const std::vector<double> uniform(count, 0.5);
        const auto result = uniform_feasible_jump_probabilities(uniform);
        const double expected = (1.0 - std::pow(0.5, count)) / count;
        for (double value : result) assert(std::abs(value - expected) < 1.0e-14);
    }
}

void lattice_transport(bool planar) {
    using namespace atcg3d;
    using namespace atcg3d::structured_pde;
    auto config = StructuredPdeConfig3D::load(std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_SharedRules/config/regular_cycle_v14.yaml");
    config.continuum.base.thin_layer = planar;
    config.continuum.grid.shape = {12, 12, planar ? 1 : 12};
    config.continuum.grid.origin = {-6.0, -6.0, planar ? -0.5 : -6.0};
    config.continuum.reaction.small_daughter_vacancy_exponent = planar ? 8.0 : 26.0;
    config.continuum.vascular.source_mode = "static_voxels";
    config.continuum.vascular.static_sources.clear();
    config.continuum.nutrient.boundary_mode = "vessels_dirichlet_v1";
    config.continuum.nutrient.initial_value = 1.0;
    config.continuum.base.migration_activation_threshold = 1.0;
    config.continuum.end_time_hours = 0.5;
    config.operator_model = "shared_operator_switches_v1";
    config.division_operator_enabled = false;
    const std::size_t count = 144 * (planar ? 1 : 12);
    StructuredInitialFields3D fields;
    fields.vessel_fraction.assign(count, 0.0);
    for (std::size_t stage = 0; stage < 2; ++stage) {
        fields.r_normal[stage].assign(count, 0.0);
        fields.r_active[stage].assign(count, 0.0);
        fields.K[stage].assign(count, 0.0);
        fields.active_remaining_hours[stage].assign(count, 0.0);
    }
    const std::size_t center = ((planar ? 0 : 6) * 12 + 6) * 12 + 6;
    fields.r_normal[0][center] = 0.25;
    StructuredPdeModel3D model(config);
    model.initialize_from_arrays(fields);
    assert(model.step());
    const double rate = config.continuum.base.normal_r_migration_beta.scale *
        config.continuum.base.normal_r_migration_beta.alpha /
        (config.continuum.base.normal_r_migration_beta.alpha + config.continuum.base.normal_r_migration_beta.beta);
    const double expected = 0.25 * config.continuum.time_step_hours * rate / (planar ? 8.0 : 26.0);
    const std::size_t diagonal = center + 13 + (planar ? 0 : 144);
    assert(std::abs(model.r_normal(StructuredStage3D::small, diagonal) - expected) < 1.0e-14);
    assert(std::abs(model.diagnostics().r_total - 0.25) < 1.0e-13);
    const auto path = std::filesystem::current_path() / "operator-checkpoint.bin";
    model.save_checkpoint(path);
    auto incompatible = config;
    incompatible.migration_operator_enabled = false;
    StructuredPdeModel3D wrong_operator(incompatible);
    bool rejected = false;
    try {
        wrong_operator.load_checkpoint(path);
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    assert(rejected);
    config.continuum.base.threads = 4;
    StructuredPdeModel3D restored(config);
    restored.load_checkpoint(path);
    assert(model.step() && restored.step());
    assert(model.state_checksum() == restored.state_checksum());
    std::filesystem::remove(path);

    config.migration_operator_enabled = false;
    config.activation_operator_enabled = false;
    fields.r_normal[0][center] = 0.0;
    fields.r_active[0][center] = 0.25;
    fields.active_remaining_hours[0][center] = 0.125;
    StructuredPdeModel3D stationary(config);
    stationary.initialize_from_arrays(fields);
    assert(stationary.step());
    assert(stationary.r_active(StructuredStage3D::small, center) == 0.0);
    assert(stationary.r_normal(StructuredStage3D::small, center) == 0.25);
    assert(stationary.r_normal(StructuredStage3D::small, diagonal) == 0.0);
    assert(stationary.diagnostics().r_total == 0.25);

    config.division_operator_enabled = true;
    config.division_clock_model = "mean_rate_v1";
    config.small_daughter_placement = "feasible_neighbor_birth_v2";
    fields.r_active[0][center] = 0.0;
    fields.r_normal[0][center] = 0.25;
    fields.active_remaining_hours[0][center] = 0.0;
    StructuredPdeModel3D birth(config);
    birth.initialize_from_arrays(fields);
    assert(birth.step());
    assert(birth.r_normal(StructuredStage3D::small, center) == 0.25);
    assert(birth.r_normal(StructuredStage3D::small, diagonal) > 0.0);
    double daughters = 0.0;
    for (std::size_t here = 0; here < count; ++here) {
        if (here != center) daughters += birth.r_normal(StructuredStage3D::small, here);
    }
    assert(std::abs(birth.diagnostics().r_total - 0.25 - daughters) < 1.0e-13);
    birth.save_checkpoint(path);
    config.continuum.base.threads = 1;
    StructuredPdeModel3D continued(config);
    continued.load_checkpoint(path);
    assert(birth.step() && continued.step());
    assert(birth.state_checksum() == continued.state_checksum());
    std::filesystem::remove(path);

    config.division_clock_model = "transported_shifted_geometric_v1";
    config.division_work_bin_width = 0.05;
    config.continuum.base.division_timing.base_cycle_hours = 0.25;
    config.continuum.base.division_timing.stochastic_time_quantum_hours = 0.01;
    config.continuum.end_time_hours = 1.0;
    StructuredPdeModel3D renewal_birth(config);
    renewal_birth.initialize_from_arrays(fields);
    assert(renewal_birth.step());
    assert(renewal_birth.step());
    assert(renewal_birth.r_normal(StructuredStage3D::small, diagonal) > 0.0);
    for (std::size_t here = 0; here < count; ++here) {
        assert(std::abs(renewal_birth.division_renewal()->mass(here, 0) -
            renewal_birth.r_normal(StructuredStage3D::small, here)) < 1.0e-13);
    }
    renewal_birth.save_checkpoint(path);
    config.continuum.base.threads = 4;
    StructuredPdeModel3D renewal_continued(config);
    renewal_continued.load_checkpoint(path);
    assert(renewal_birth.step() && renewal_continued.step());
    assert(renewal_birth.state_checksum() == renewal_continued.state_checksum());
    std::filesystem::remove(path);
}

void abm_switches() {
    using namespace atcg3d;
    using namespace atcg3d::shared_rules;
    auto config = structured_pde::StructuredPdeConfig3D::load(std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_SharedRules/config/regular_cycle_v14.yaml");
    config.operator_model = "shared_operator_switches_v1";
    config.migration_operator_enabled = false;
    config.division_operator_enabled = false;
    config.activation_operator_enabled = false;
    auto environment = std::make_unique<SharedResourceEnvironment3D>(config);
    Simulation3D simulation(abm_config(config), std::move(environment));
    simulation.initialize();
    const auto initial = simulation.snapshot_cells();
    simulation.run();
    assert(simulation.stats().migration_attempts == 0);
    assert(simulation.stats().divisions == 0);
    assert(simulation.cells().alive_slots().size() == initial.size());
    for (const auto slot : simulation.cells().alive_slots()) {
        assert(simulation.cells().anchor(slot) == initial[slot].anchor);
        assert(!(simulation.cells().flags(slot) & kMigrationActive));
    }
    config.schema_version = 13;
    bool rejected = false;
    try {
        config.validate();
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    assert(rejected);
}

}  // namespace

int main() {
    feasible_choice();
    lattice_transport(true);
    lattice_transport(false);
    abm_switches();
}
