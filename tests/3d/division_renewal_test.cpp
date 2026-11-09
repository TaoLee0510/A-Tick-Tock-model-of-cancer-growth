#include <cassert>
#include <cmath>
#include <filesystem>
#include <sstream>

#include "model/division_renewal.hpp"
#include "model/structured_pde_model.hpp"

int main() {
    using namespace atcg3d;
    using namespace atcg3d::structured_pde;
    DivisionTimingConfig timing;
    timing.minimum_fraction = 0.9;
    timing.stochastic_tail_fraction = 0.1;
    DivisionRenewal3D renewal(timing, 0.5, 128.0, {1.0, 0.6});
    renewal.add_fresh(0, 0, 1.0);
    assert(std::abs(renewal.mass(0, 0) - 1.0) < 1.0e-12);
    // The literal shifted-geometric work law has mean 24 unit-rate hours.
    assert(std::abs(renewal.mean_work(0, 0) - 24.0) < 1.0e-9);
    assert(renewal.advance(0, 0, 10.0) == 0.0);
    assert(std::abs(renewal.mean_work(0, 0) - 14.0) < 1.0e-9);
    renewal.begin_transport();
    renewal.transfer(0, 1, 0, 0.4);
    renewal.finish_transport([](auto location, auto channel) {
        return channel == 0 ? (location == 0 ? 0.6 : 0.4) : 0.0;
    });
    assert(std::abs(renewal.mass(1, 0) - 0.4) < 1.0e-12);
    assert(std::abs(renewal.mean_work(1, 0) - 14.0) < 1.0e-9);
    std::stringstream checkpoint(std::ios::in | std::ios::out | std::ios::binary);
    renewal.save(checkpoint);
    DivisionRenewal3D restored(timing, 0.5, 128.0, {1.0, 0.6});
    restored.load(checkpoint, 2);
    assert(restored.checksum() == renewal.checksum());
    const double completed = renewal.advance(1, 0, 128.0);
    assert(std::abs(completed - 0.4) < 1.0e-12);
    assert(renewal.mass(1, 0) == 0.0);

    // Independent renewal equation: E Z(t) = 1 + sum_g 2^(g-1) F_g(t).
    // With constant unit growth, shifted geometric cycles have a closed-form
    // negative-binomial convolution. Aligning work/time bins removes numerical
    // advection dispersion from this check.
    DivisionRenewal3D well_mixed(timing, 0.05, 128.0, {1.0, 0.6});
    well_mixed.add_fresh(0, 0, 1.0);
    for (int step = 0; step < 960; ++step) {
        const double events = well_mixed.advance(0, 0, 0.05);
        well_mixed.add_fresh(0, 0, 2.0 * events);
    }
    const double p = 1.0 / 2.4;
    double exact = 1.0;
    for (int generation = 1; generation <= 2; ++generation) {
        const int maximum = int(std::floor(48.0 - generation * 21.6));
        double cdf = 0.0;
        for (int k = generation; k <= maximum; ++k) {
            const double choose = generation == 1 ? 1.0 : k - 1.0;
            cdf += choose * std::pow(p, generation) * std::pow(1.0 - p, k - generation);
        }
        exact += std::pow(2.0, generation - 1) * cdf;
    }
    assert(std::abs(well_mixed.mass(0, 0) - exact) < 1.0e-9);

    auto config = StructuredPdeConfig3D::load(std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_SharedRules/config/regular_cycle_v10.yaml");
    config.continuum.grid.shape = {12, 12, 1};
    config.continuum.grid.origin = {-6.0, -6.0, -0.5};
    config.continuum.vascular.source_mode = "static_voxels";
    config.continuum.vascular.static_sources.clear();
    config.continuum.end_time_hours = 16;
    config.continuum.base.growth_density_window_edge = 1;
    config.continuum.base.r_to_K_conversion.enabled = false;
    StructuredInitialFields3D fields;
    for (std::size_t s = 0; s < 2; ++s) {
        fields.r_normal[s].resize(144);
        fields.r_active[s].resize(144);
        fields.K[s].resize(144);
        fields.active_remaining_hours[s].resize(144);
    }
    fields.vessel_fraction.resize(144);
    fields.r_normal[0][78] = 0.01;
    fields.K[0][78] = 0.01;
    StructuredPdeModel3D original(config);
    original.initialize_from_arrays(fields);
    for (int tick = 0; tick < 32; ++tick) assert(original.step());
    const auto d = original.diagnostics();
    assert(std::abs(d.r_total + d.K_total - 0.02) < 1.0e-9);
    const auto path = std::filesystem::current_path() / "division-renewal.bin";
    std::filesystem::remove(path);
    original.save_checkpoint(path);
    auto parallel = config;
    parallel.continuum.base.threads = 4;
    StructuredPdeModel3D resumed(parallel);
    resumed.load_checkpoint(path);
    assert(original.state_checksum() == resumed.state_checksum());
    while (original.step()) assert(resumed.step());
    assert(!resumed.step());
    assert(original.state_checksum() == resumed.state_checksum());
    std::filesystem::remove(path);
}
