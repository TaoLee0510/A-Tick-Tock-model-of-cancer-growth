#include <algorithm>
#include <cassert>
#include <cmath>
#include <filesystem>

#include "ode/ode_model.hpp"
#include "model/structured_pde_model.hpp"

int main() {
    using namespace atcg3d;
    using namespace atcg3d::ode;
    auto config=OdeConfig3D::load(std::filesystem::path(ATCG_SOURCE_DIR)/"ATCG3D_ODE/config/ode_smoke_v1.yaml");
    auto& c=config.rules.continuum;
    c.grid.shape={6,6,1}; c.grid.origin={-3,-3,-0.5};
    c.base.growth_density_window_edge=1;
    c.base.r_to_K_conversion.enabled=false;
    c.end_time_hours=2; c.time_step_hours=0.01;
    c.nutrient.boundary_mode="vessels_dirichlet_v1";
    c.nutrient.initial_value=1;
    c.nutrient.vessel_exchange_per_hour=0;
    c.nutrient.decay_per_hour=0;
    c.nutrient.r_consumption_rate_per_hour=c.nutrient.K_consumption_rate_per_hour=0;
    c.vascular.source_mode="abm_perfusion";
    config.initial={0.03,0.002,0,0,0.02,0.003,0,0,1,0,0,0,0,0};
    OdeModel3D ode(config);
    PeriodicReactionPde3D periodic(config);
    while(ode.step()) assert(periodic.step());
    const auto mean=periodic.mean();
    for(std::size_t i=0;i<14;++i) assert(std::abs(mean[i]-ode.state()[i])<2.0e-6);
    for(const auto& field:periodic.fields()) for(std::size_t i=0;i<14;++i) assert(std::abs(field[i]-mean[i])<1.0e-14);

    OdeModel3D first(config);
    for(int i=0;i<7;++i) assert(first.step());
    const auto checkpoint=std::filesystem::temp_directory_path()/"atcg_ode_test.bin";
    std::filesystem::remove(checkpoint);
    first.save_checkpoint(checkpoint);
    auto parallel=config; parallel.rules.continuum.base.threads=4;
    OdeModel3D resumed(parallel); resumed.load_checkpoint(checkpoint);
    assert(resumed.state_checksum()==first.state_checksum());
    first.run(); resumed.run();
    assert(resumed.state()==first.state());
    assert(resumed.state_checksum()==ode.state_checksum());
    std::filesystem::remove(checkpoint);

    // Independently compare the published structured-v7 Euler reaction, at
    // a small time step, against the continuous RK45 mean-rate limit.
    c.time_step_hours=0.001;
    OdeModel3D fine(config);
    structured_pde::StructuredPdeModel3D legacy(config.rules);
    structured_pde::StructuredInitialFields3D fields;
    for(int stage=0;stage<2;++stage) {
        fields.r_normal[stage].assign(36,config.initial[stage]);
        fields.r_active[stage].assign(36,0);
        fields.K[stage].assign(36,config.initial[4+stage]);
        fields.active_remaining_hours[stage].assign(36,0);
    }
    fields.vessel_fraction.assign(36,0);
    legacy.initialize_from_arrays(fields);
    while(fine.step()) assert(legacy.step());
    const auto diagnostics=legacy.diagnostics();
    for(int stage=0;stage<2;++stage) {
        assert(std::abs(diagnostics.r_normal_mass[stage]/36-fine.state()[stage])<5.0e-5);
        assert(std::abs(diagnostics.K_mass[stage]/36-fine.state()[4+stage])<5.0e-5);
    }

    // RK45 accuracy against the analytic nutrient decay solution.
    auto frozen=config;
    frozen.rules.continuum.base.division_timing.base_cycle_hours=1.0e30;
    frozen.rules.continuum.base.death_growth_rate_threshold=-1.0e30;
    frozen.rules.continuum.nutrient.decay_per_hour=2.0;
    State state{}; state[8]=1;
    Reaction3D decay(frozen.rules); decay.advance(state,3,config.solver);
    assert(std::abs(state[8]-std::exp(-6.0))<1.0e-10);

    auto stiff=frozen;
    stiff.rules.continuum.base.death_growth_rate_threshold=10;
    stiff.rules.continuum.base.r_death_delay_hours=0.01;
    state={}; state[0]=0.1;
    Reaction3D mortality(stiff.rules); mortality.advance(state,0.2,config.solver);
    assert(state[0]>=0 && std::abs(state[0]-0.1*std::exp(-20.0))<1.0e-11);

    // Expiry/refractory transitions conserve cells and keep nonnegative time
    // mass, including multiple transitions inside one requested interval.
    frozen.rules.continuum.nutrient.decay_per_hour=0;
    frozen.rules.migration.reactivation_cooldown_hours=0.5;
    Reaction3D reaction(frozen.rules);
    state={}; state[2]=0.1; state[6]=0.04; state[8]=1;
    reaction.advance(state,1.0,config.solver);
    assert(std::abs(reaction.total_mass(state)-0.1)<1.0e-12);
    assert(state[2]==0 && state[6]==0 && state[10]==0 && state[12]==0);
    assert(std::abs(state[0]-0.1)<1.0e-12);

    // Conservative periodic transport carries refractory clock with its mass.
    frozen.rules.continuum.end_time_hours=0.5;
    frozen.rules.continuum.time_step_hours=0.1;
    frozen.initial={}; frozen.initial[8]=1;
    PeriodicReactionPde3D transport(frozen);
    std::vector<State> heterogeneous(36,frozen.initial);
    heterogeneous[0][0]=heterogeneous[0][10]=0.1;
    heterogeneous[0][12]=0.1;
    transport.initialize(heterogeneous);
    while(transport.step()) {}
    const auto conserved=transport.mean();
    assert(std::abs(reaction.total_mass(conserved)-0.1/36)<1.0e-12);
    assert(std::abs(conserved[10]-0.1/36)<1.0e-12);
    assert(std::abs(conserved[12]-0.05/36)<1.0e-12);
    for(const auto& field:transport.fields()) for(double v:field) assert(std::isfinite(v) && v>=0);
}
