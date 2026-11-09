#include <cassert>
#include <cmath>
#include <filesystem>
#include <memory>

#include "model/hybrid_model.hpp"

int main() {
    using namespace atcg3d;
    using namespace atcg3d::hybrid;
    auto c = HybridConfig3D::load(std::filesystem::path(ATCG_SOURCE_DIR) /
                                  "ATCG3D_Hybrid/config/hybrid_smoke_v1.yaml");
    c.rules.continuum.end_time_hours = 2;
    c.exchange_every_hours = 0.5;
    c.rules.continuum.angiogenesis.model = "disabled";
    c.rules.continuum.base.angiogenesis.enabled = false;
    c.rules.continuum.vascular.source_mode = "static_voxels";
    c.rules.continuum.vascular.static_sources.clear();
    c.rules.continuum.base.initial_r_cells = 32;
    c.rules.continuum.base.initial_K_cells = 32;
    c.rules.continuum.base.initial_large_fraction = 0;
    c.rules.continuum.base.initial_radius = 8;
    c.rules.continuum.base.initial_shell_inner_radius = 4;
    c.rules.continuum.base.initial_inner_small_radius = 5;
    c.rules.continuum.base.division_timing.base_cycle_hours = 1e30;
    c.rules.continuum.base.growth_density_window_edge = 1;
    c.rules.continuum.base.r_death_delay_hours = 1e30;
    c.rules.continuum.base.K_death_delay_hours = 1e30;
    c.rules.continuum.base.migration_activation_threshold = 1.0;
    c.rules.continuum.migration.diffusion_scale = 1e-12;
    c.smoothing_radius = 2;
    c.core_on = 0.30;
    c.core_off = 0.20;
    HybridModel3D mixed(c);
    mixed.initialize();
    const auto mass = mixed.diagnostics().total_mass;
    assert(mass == 64);
    assert(mixed.diagnostics().to_pde > 0);
    assert(mixed.diagnostics().abm_mass > 0);
    for (int i = 0; i < 2; ++i) {
        assert(mixed.step());
        assert(std::abs(mixed.diagnostics().total_mass - mass) < 1e-7);
    }
    const auto path = std::filesystem::current_path() / "atcg-hybrid-test.bin";
    for (auto suffix : {"", ".pde.bin", ".resource.bin"})
        std::filesystem::remove(path.string() + suffix);
    mixed.save_checkpoint(path);
    auto parallel = c;
    parallel.rules.continuum.base.threads = 4;
    HybridModel3D resumed(parallel);
    resumed.load_checkpoint(path);
    assert(resumed.state_checksum() == mixed.state_checksum());
    while (mixed.step()) {
    }
    while (resumed.step()) {
    }
    assert(mixed.state_checksum() == resumed.state_checksum());
    assert(std::abs(mixed.diagnostics().total_mass - mass) < 1e-7);
    // Exact reductions use the existing engines and preserve all native state.
    auto full = c;
    full.mode = "all_abm";
    HybridModel3D all_abm(full);
    while (all_abm.step()) {
    }
    auto env =
        std::make_unique<shared_rules::SharedResourceEnvironment3D>(c.rules);
    Simulation3D reference(shared_rules::abm_config(c.rules), std::move(env));
    reference.run();
    assert(reference.state_checksum() == all_abm.state_checksum());
    full.mode = "all_pde";
    HybridModel3D all_pde(full);
    while (all_pde.step()) {
    }
    auto seed_env =
        std::make_unique<shared_rules::SharedResourceEnvironment3D>(c.rules);
    Simulation3D seed(shared_rules::abm_config(c.rules), std::move(seed_env));
    seed.initialize();
    structured_pde::StructuredPdeModel3D pde(c.rules);
    pde.initialize_from_abm(seed);
    while (pde.step()) {
    }
    assert(pde.state_checksum() == all_pde.state_checksum());
    auto dilute = c;
    dilute.smoothing_radius = 0;
    dilute.core_on = 0.9;
    dilute.core_off = 0.8;
    HybridModel3D rounding(dilute);
    structured_pde::StructuredInitialFields3D fields;
    const auto voxels = rounding.pde().voxel_count();
    for (int stage = 0; stage < 2; ++stage) {
        fields.r_normal[stage].resize(voxels);
        fields.r_active[stage].resize(voxels);
        fields.K[stage].resize(voxels);
        fields.active_remaining_hours[stage].resize(voxels);
        fields.r_refractory[stage].resize(voxels);
        fields.refractory_remaining_hours[stage].resize(voxels);
    }
    fields.vessel_fraction.resize(voxels);
    for (int i = 0; i < 4; ++i) {
        auto here = std::size_t(20 * 48 + 20 + i);
        fields.r_normal[0][here] = 0.6;
        fields.r_refractory[0][here] = 0.3;
        fields.refractory_remaining_hours[0][here] = 10;
    }
    rounding.initialize_from_arrays(std::move(fields));
    const double before = rounding.diagnostics().total_mass;
    rounding.exchange();
    assert(rounding.diagnostics().to_abm == 2);
    assert(std::abs(rounding.diagnostics().total_mass - before) < 1e-12);
    for (auto suffix : {"", ".pde.bin", ".resource.bin"})
        std::filesystem::remove(path.string() + suffix);
    rounding.save_checkpoint(path);
    HybridModel3D rounding_resume(dilute);
    rounding_resume.load_checkpoint(path);
    assert(rounding_resume.state_checksum() == rounding.state_checksum());
    for (auto suffix : {"", ".pde.bin", ".resource.bin"})
        std::filesystem::remove(path.string() + suffix);
    auto volume = dilute;
    volume.rules.continuum.base.thin_layer = false;
    volume.rules.continuum.reaction.small_daughter_vacancy_exponent = 26;
    volume.rules.continuum.grid.shape = {12, 12, 12};
    volume.rules.continuum.grid.origin = {-6, -6, -6};
    HybridModel3D three(volume);
    structured_pde::StructuredInitialFields3D active;
    for (int stage = 0; stage < 2; ++stage) {
        active.r_normal[stage].resize(1728);
        active.r_active[stage].resize(1728);
        active.K[stage].resize(1728);
        active.active_remaining_hours[stage].resize(1728);
    }
    active.vessel_fraction.resize(1728);
    for (int z = 6; z < 8; ++z)
        for (int y = 6; y < 8; ++y)
            for (int x = 6; x < 10; ++x) {
                auto i = std::size_t((z * 12 + y) * 12 + x);
                active.r_active[1][i] = 0.08;
                active.active_remaining_hours[1][i] = 5;
            }
    three.initialize_from_arrays(std::move(active));
    const double active_before = three.diagnostics().total_mass;
    three.exchange();
    assert(three.diagnostics().to_abm == 1);
    assert(std::abs(three.diagnostics().total_mass - active_before) < 1e-7);
    const auto agent = three.abm().cells().alive_slots().front();
    assert(three.abm().cells().stage(agent) == CellStage::large);
    assert(three.abm().cells().migration_activation_end_time(agent) > 4.99);
    three.save_checkpoint(path);
    HybridModel3D three_resume(volume);
    three_resume.load_checkpoint(path);
    assert(three_resume.state_checksum() == three.state_checksum());
    for (auto suffix : {"", ".pde.bin", ".resource.bin"})
        std::filesystem::remove(path.string() + suffix);
    auto vascular =
        HybridConfig3D::load(std::filesystem::path(ATCG_SOURCE_DIR) /
                             "ATCG3D_Hybrid/config/hybrid_smoke_v1.yaml");
    HybridModel3D coupled(vascular);
    while (coupled.step()) {
    }
    assert(coupled.diagnostics().to_pde > 0 &&
           coupled.diagnostics().to_abm > 0);
    assert(coupled.pde().angiogenesis() != nullptr);
    coupled.save_checkpoint(path);
    HybridModel3D coupled_resume(vascular);
    coupled_resume.load_checkpoint(path);
    assert(coupled.state_checksum() == coupled_resume.state_checksum());
    for (auto suffix : {"", ".pde.bin", ".resource.bin"})
        std::filesystem::remove(path.string() + suffix);
    vascular.mode = "all_abm";
    vascular.rules.continuum.end_time_hours = 2;
    HybridModel3D vascular_abm(vascular);
    while (vascular_abm.step()) {
    }
    vascular_abm.save_checkpoint(path);
    vascular.rules.continuum.end_time_hours = 8;
    HybridModel3D vascular_resume(vascular);
    vascular_resume.load_checkpoint(path);
    while (vascular_resume.step()) {
    }
    HybridModel3D vascular_whole(vascular);
    while (vascular_whole.step()) {
    }
    assert(vascular_whole.state_checksum() == vascular_resume.state_checksum());
    for (auto suffix : {"", ".pde.bin", ".resource.bin"})
        std::filesystem::remove(path.string() + suffix);
}
