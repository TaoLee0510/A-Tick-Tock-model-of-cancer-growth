#include <cassert>
#include <cmath>
#include <filesystem>

#include "model/hybrid_model.hpp"

namespace {
atcg3d::structured_pde::StructuredInitialFields3D empty_fields(std::size_t size) {
    atcg3d::structured_pde::StructuredInitialFields3D f;
    for (int s = 0; s < 2; ++s) {
        f.r_normal[s].resize(size);
        f.r_active[s].resize(size);
        f.K[s].resize(size);
        f.active_remaining_hours[s].resize(size);
    }
    f.vessel_fraction.resize(size);
    return f;
}
void remove_checkpoint(const std::filesystem::path& path) {
    for (auto suffix : {"", ".pde.bin", ".resource.bin"})
        std::filesystem::remove(path.string() + suffix);
}
double work_moment(const atcg3d::hybrid::HybridModel3D& model) {
    double moment = 0.0;
    const auto* bank = model.pde().division_renewal();
    for (std::size_t i = 0; i < model.pde().voxel_count(); ++i)
        for (int c = 0; c < 4; ++c) moment += bank->mass(i, c) * bank->mean_work(i, c);
    for (auto s : model.abm().cells().alive_slots())
        moment += model.abm().cells().division_work_remaining(s);
    return moment;
}
}

int main() {
    using namespace atcg3d;
    using namespace atcg3d::hybrid;
    const auto root = std::filesystem::path(ATCG_SOURCE_DIR);
    auto c = HybridConfig3D::load(root / "ATCG3D_Hybrid/config/hybrid_regular_cycle_v2.yaml");
    c.rules.continuum.grid.shape = {12, 12, 1};
    c.rules.continuum.grid.origin = {-6, -6, -0.5};
    c.rules.continuum.end_time_hours = 8;
    c.exchange_every_hours = 12;
    c.smoothing_radius = 0;
    c.core_on = 0.9;
    c.core_off = 0.8;
    c.rules.continuum.base.initial_growth_rate_model = "fixed";
    c.rules.continuum.base.initial_r_growth_rate = 1;
    c.rules.continuum.base.initial_K_growth_rate = 1;
    c.rules.continuum.base.growth_density_window_edge = 1;
    c.rules.continuum.base.migration_activation_threshold = 1;
    c.rules.continuum.base.r_death_delay_hours = 1e30;
    c.rules.continuum.migration.diffusion_scale = 1e-12;
    c.rules.continuum.nutrient.K_consumption_rate_per_hour = 0;
    c.rules.continuum.nutrient.r_consumption_rate_per_hour = 0;
    c.rules.continuum.nutrient.decay_per_hour = 0;
    c.rules.continuum.vascular.source_mode = "static_voxels";
    c.rules.continuum.vascular.static_sources.clear();
    auto dilute = empty_fields(144);
    for (int i = 0; i < 4; ++i) dilute.r_normal[0][78 + i] = 0.2;
    HybridModel3D tail(c);
    tail.initialize_from_arrays(std::move(dilute));
    const auto untouched = tail.pde().division_renewal()->checksum();
    for (int n = 0; n < 10; ++n) tail.exchange();
    assert(tail.diagnostics().to_abm == 0 && tail.diagnostics().to_pde == 0);
    assert(tail.pde().division_renewal()->checksum() == untouched);
    for (int i = 0; i < 4; ++i)
        assert(tail.pde().r_normal(structured_pde::StructuredStage3D::small, 78 + i) == 0.2);
    auto fields = empty_fields(144);
    for (int i = 0; i < 4; ++i) fields.r_normal[0][78 + i] = 0.3;
    HybridModel3D aged(c);
    aged.initialize_from_arrays(std::move(fields));
    for (int i = 0; i < 24; ++i) assert(aged.step());
    assert(std::abs(aged.diagnostics().total_mass - 1.2) < 1e-9);
    aged.exchange();
    assert(aged.diagnostics().to_abm == 1);
    const auto slot = aged.abm().cells().alive_slots().front();
    // A freshly drawn cycle cannot be below its minimum work of 21.6.
    assert(aged.abm().cells().division_work_remaining(slot) < 21.6);
    const double carried_work = work_moment(aged);
    aged.exchange();
    assert(aged.diagnostics().to_pde == 1);
    assert(std::abs(work_moment(aged) - carried_work) < 1e-5);
    assert(std::abs(aged.diagnostics().total_mass - 1.2) < 1e-9);

    // Full activation/time/rate distributions survive mixed restart.
    auto a = c;
    a.rules = structured_pde::StructuredPdeConfig3D::load(root /
        "ATCG3D_SharedRules/config/active_r200_ci_v12.yaml");
    a.rules.division_clock_model = "transported_shifted_geometric_v1";
    a.rules.continuum.grid.shape = {12, 12, 1};
    a.rules.continuum.grid.origin = {-6, -6, -0.5};
    a.rules.continuum.end_time_hours = 2;
    a.rules.continuum.vascular.source_mode = "static_voxels";
    a.rules.continuum.vascular.static_sources.clear();
    a.exchange_every_hours = 0.5;
    auto active = empty_fields(144);
    for (int i = 0; i < 4; ++i) {
        active.r_active[0][78 + i] = 0.3;
        active.active_remaining_hours[0][78 + i] = 5;
    }
    HybridModel3D mixed(a);
    mixed.initialize_from_arrays(std::move(active));
    const double active_mass = mixed.diagnostics().total_mass;
    mixed.exchange();
    assert(std::abs(mixed.diagnostics().total_mass - active_mass) < 1e-12);
    assert(mixed.diagnostics().to_abm == 1);
    const auto active_slot = mixed.abm().cells().alive_slots().front();
    assert(mixed.abm().cells().migration_activation_end_time(active_slot) == 5);
    assert(mixed.abm().cells().normal_migration_rate(active_slot) > 0);
    assert(mixed.step());
    const auto path = std::filesystem::current_path() / "hybrid-distribution.bin";
    remove_checkpoint(path);
    mixed.save_checkpoint(path);
    a.rules.continuum.base.threads = 4;
    HybridModel3D resumed(a);
    resumed.load_checkpoint(path);
    assert(mixed.state_checksum() == resumed.state_checksum());
    while (mixed.step()) assert(resumed.step());
    assert(!resumed.step());
    assert(mixed.state_checksum() == resumed.state_checksum());
    // Native float directional fields accumulate substep roundoff; exchange
    // itself is checked above with a double-precision mass budget.
    assert(std::abs(mixed.diagnostics().total_mass - active_mass) < 1e-5);
    remove_checkpoint(path);

    a.rules.continuum.base.thin_layer = false;
    a.rules.continuum.nutrient.boundary_mode = "vessels_dirichlet_v1";
    a.rules.continuum.reaction.small_daughter_vacancy_exponent = 26;
    a.rules.continuum.grid.shape = {12, 12, 12};
    a.rules.continuum.grid.origin = {-6, -6, -6};
    auto volume = empty_fields(1728);
    for (int z = 6; z < 8; ++z) for (int y = 6; y < 8; ++y) for (int x = 6; x < 10; ++x) {
        const auto i = std::size_t((z * 12 + y) * 12 + x);
        volume.r_active[1][i] = 0.08;
        volume.active_remaining_hours[1][i] = 5;
    }
    HybridModel3D three(a);
    three.initialize_from_arrays(std::move(volume));
    three.exchange();
    assert(three.diagnostics().to_abm == 1);
    assert(std::abs(three.diagnostics().total_mass - 1.28) < 1e-7);
    const auto large = three.abm().cells().alive_slots().front();
    assert(three.abm().cells().stage(large) == CellStage::large);
    assert(three.abm().cells().migration_activation_end_time(large) == 5);
    three.save_checkpoint(path);
    HybridModel3D restored(a);
    restored.load_checkpoint(path);
    assert(three.state_checksum() == restored.state_checksum());
    remove_checkpoint(path);

    auto native = c;
    native.rules.continuum.end_time_hours = 0.5;
    native.rules.continuum.base.initial_r_cells = 8;
    native.rules.continuum.base.initial_K_cells = 8;
    native.rules.continuum.base.initial_radius = 4;
    native.rules.continuum.base.initial_shell_thickness = 1;
    native.rules.continuum.base.initial_shell_inner_radius = 3;
    native.rules.continuum.base.initial_inner_small_radius = 3;
    native.mode = "all_abm";
    HybridModel3D all_abm(native);
    while (all_abm.step()) {}
    Simulation3D reference(shared_rules::abm_config(native.rules),
        std::make_unique<shared_rules::SharedResourceEnvironment3D>(native.rules));
    reference.run();
    assert(all_abm.state_checksum() == reference.state_checksum());
    native.mode = "all_pde";
    HybridModel3D all_pde(native);
    while (all_pde.step()) {}
    Simulation3D initial(shared_rules::abm_config(native.rules),
        std::make_unique<shared_rules::SharedResourceEnvironment3D>(native.rules));
    initial.initialize();
    structured_pde::StructuredPdeModel3D density(native.rules);
    density.initialize_from_abm(initial);
    while (density.step()) {}
    assert(all_pde.state_checksum() == density.state_checksum());
}
