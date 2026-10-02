#include <cassert>
#include <cmath>
#include <filesystem>
#include <numeric>

#include "model/beta_duration.hpp"
#include "model/division_renewal.hpp"
#include "model/structured_pde_model.hpp"

int main() {
    using namespace atcg3d::structured_pde;
    for (int i = 1; i < 100; ++i) {
        const double x = i / 100.0;
        assert(std::abs(beta_cdf(x, 1, 1) - x) < 1e-13);
        assert(std::abs(beta_cdf(x, 1, 3) - (1 - std::pow(1 - x, 3))) < 1e-13);
        assert(std::abs(beta_cdf(x, 3, 1) - std::pow(x, 3)) < 1e-13);
        assert(std::abs(beta_cdf(x, 0.005, 0.011666666666666667) +
                        beta_cdf(1 - x, 0.011666666666666667, 0.005) - 1) < 1e-12);
    }
    for (const auto [a, b] : {std::pair{1.0, 1.0}, std::pair{2.0, 4.666666666666667},
                              std::pair{0.005, 0.011666666666666667}}) {
        const auto law = beta_duration_kernel(a, b, 24, 0.5, 32);
        assert(std::abs(std::accumulate(law.begin(), law.end(), 0.0) - 1) < 1e-12);
        double mean = 0;
        for (std::size_t i = 0; i < law.size(); ++i) mean += i * 0.5 * law[i];
        assert(std::abs(mean - 24 * a / (a + b)) < 1e-10);
    }
    const auto law = beta_duration_kernel(1, 1, 24, 0.5, 32);
    DivisionRenewal3D mixture(0.5, 32, {law, law});
    mixture.add(0, 0, 0.3, 1);
    mixture.add(0, 0, 0.7, 9);
    assert(std::abs(mixture.advance(0, 0, 4) - 0.3) < 1e-12);
    assert(std::abs(mixture.mass(0, 0) - 0.7) < 1e-12);
    assert(std::abs(mixture.mean_work(0, 0) - 5) < 1e-12);
    mixture.begin_transport();
    mixture.transfer(0, 1, 0, 0.2);
    mixture.finish_transport([](auto i, auto c) { return c ? 0.0 : i ? 0.2 : 0.5; });
    assert(std::abs(mixture.mean_work(1, 0) - 5) < 1e-12);

    DivisionRenewal3D selective(0.5, 32, {law, law});
    selective.add(0, 0, 0.5, 1);
    selective.add(0, 0, 0.5, 2);
    std::array<std::vector<double>, 4> weights;
    weights[0].resize(law.size());
    weights[0][2] = 0.2;
    weights[0][4] = 0.4;
    selective.begin_weighted_transport(weights);
    assert(std::abs(selective.weighted_transport_mass(0, 0) - 0.3) < 1e-12);
    selective.transfer(0, 1, 0, 0.15);
    selective.finish_transport([](auto i, auto c) { return c ? 0.0 : i ? 0.15 : 0.85; });
    assert(std::abs(selective.mean_work(1, 0) - 5.0 / 3.0) < 1e-12);
    assert(std::abs(selective.mass(0, 0) + selective.mass(1, 0) - 1) < 1e-12);

    auto config = StructuredPdeConfig3D::load(std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_SharedRules/config/active_r200_ci_v12.yaml");
    config.division_clock_model = "transported_shifted_geometric_v1";
    config.continuum.grid.shape = {12, 12, 1};
    config.continuum.grid.origin = {-6, -6, -0.5};
    config.continuum.end_time_hours = 2;
    config.continuum.base.growth_density_window_edge = 1;
    config.continuum.vascular.source_mode = "static_voxels";
    config.continuum.vascular.static_sources.clear();
    StructuredInitialFields3D fields;
    for (int stage = 0; stage < 2; ++stage) {
        fields.r_normal[stage].resize(144);
        fields.r_active[stage].resize(144);
        fields.K[stage].resize(144);
        fields.active_remaining_hours[stage].resize(144);
    }
    fields.vessel_fraction.resize(144);
    fields.r_active[0][78] = fields.r_active[0][79] = 0.01;
    fields.active_remaining_hours[0][78] = 0.5;
    fields.active_remaining_hours[0][79] = 5;
    StructuredPdeModel3D original(config);
    original.initialize_from_arrays(std::move(fields));
    assert(original.step());
    const auto path = std::filesystem::current_path() / "activation-distribution.bin";
    std::filesystem::remove(path);
    original.save_checkpoint(path);
    config.continuum.base.threads = 4;
    StructuredPdeModel3D resumed(config);
    resumed.load_checkpoint(path);
    assert(original.state_checksum() == resumed.state_checksum());
    while (original.step()) assert(resumed.step());
    assert(!resumed.step());
    assert(original.state_checksum() == resumed.state_checksum());
    const auto d = original.diagnostics();
    assert(std::abs(d.r_total - 0.02) < 1e-7);
    assert(d.r_active_total > 0.001 && d.r_active_total < 0.015);
    assert(d.r_normal_mass[0] > 0.005);
    std::filesystem::remove(path);
}
