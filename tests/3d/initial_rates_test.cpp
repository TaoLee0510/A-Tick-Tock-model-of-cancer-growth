#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>

#include "config/model_config.hpp"
#include "rules/initial_rates.hpp"

namespace {

double normal_cdf(double value) {
    constexpr double inverse_sqrt_two = 0.70710678118654752440;
    return 0.5 * std::erfc(-value * inverse_sqrt_two);
}

double normal_pdf(double value) {
    constexpr double inverse_sqrt_two_pi = 0.39894228040143267794;
    return inverse_sqrt_two_pi * std::exp(-0.5 * value * value);
}

double truncated_normal_expected_mean(const atcg3d::TruncatedNormalRateConfig& parameters) {
    const double lower = (parameters.minimum - parameters.mean) /
                         parameters.standard_deviation;
    const double upper = (parameters.maximum - parameters.mean) /
                         parameters.standard_deviation;
    return parameters.mean + parameters.standard_deviation *
        (normal_pdf(lower) - normal_pdf(upper)) /
        (normal_cdf(upper) - normal_cdf(lower));
}

}  // namespace

int main() {
    using namespace atcg3d;

    const Model3DConfig smoke = Model3DConfig::load(
        std::filesystem::path("configs") / "atcg3d_smoke_test_v2.yaml");
    for (CellUid uid = 1; uid <= 100; ++uid) {
        const InitialCellRates3D r = sample_initial_cell_rates(
            smoke, CellType::r, smoke.seed, uid);
        const InitialCellRates3D K = sample_initial_cell_rates(
            smoke, CellType::K, smoke.seed, uid);
        assert(r.inherent_growth_rate == 1.0 && r.migration_rate == 0.25);
        assert(K.inherent_growth_rate == 1.0 && K.migration_rate == 0.25);
    }

    const Model3DConfig production = Model3DConfig::load(
        std::filesystem::path("configs") / "atcg3d_legacy_2d_mapped_v2.yaml");
    constexpr std::uint64_t sample_count = 100000;
    double r_growth_sum = 0.0;
    double K_growth_sum = 0.0;
    double r_migration_sum = 0.0;
    double K_migration_sum = 0.0;
    std::uint64_t r_clamped_count = 0;
    std::uint64_t r_unclamped_count = 0;
    bool seed_changes_growth = false;
    bool seed_changes_migration = false;

    for (CellUid uid = 1; uid <= sample_count; ++uid) {
        const InitialCellRates3D first = sample_initial_cell_rates(
            production, CellType::r, production.seed, uid);
        const InitialCellRates3D repeated = sample_initial_cell_rates(
            production, CellType::r, production.seed, uid);
        assert(first.inherent_growth_rate == repeated.inherent_growth_rate);
        assert(first.migration_rate == repeated.migration_rate);

        // Calling the two independent streams in the opposite order cannot
        // perturb either value because no mutable RNG state is consumed.
        const double migration_first = sample_initial_migration_rate(
            production, CellType::r, production.seed, uid);
        const double growth_second = sample_initial_growth_rate(
            production, CellType::r, production.seed, uid);
        assert(migration_first == first.migration_rate);
        assert(growth_second == first.inherent_growth_rate);

        const InitialCellRates3D K = sample_initial_cell_rates(
            production, CellType::K, production.seed, uid);
        assert(first.inherent_growth_rate >=
               production.initial_r_growth_truncated_normal.minimum);
        assert(first.inherent_growth_rate <=
               production.initial_r_growth_truncated_normal.maximum);
        assert(K.inherent_growth_rate >=
               production.initial_K_growth_truncated_normal.minimum);
        assert(K.inherent_growth_rate <=
               production.initial_K_growth_truncated_normal.maximum);
        assert(first.migration_rate == 0.25 || first.migration_rate > 0.5);
        assert(first.migration_rate <= 200.0);
        assert(K.migration_rate >= 0.0 && K.migration_rate <= 0.25);

        r_clamped_count += first.migration_rate == 0.25 ? 1 : 0;
        r_unclamped_count += first.migration_rate > 0.5 ? 1 : 0;
        r_growth_sum += first.inherent_growth_rate;
        K_growth_sum += K.inherent_growth_rate;
        r_migration_sum += first.migration_rate;
        K_migration_sum += K.migration_rate;

        if (uid <= 128) {
            seed_changes_growth = seed_changes_growth ||
                sample_initial_growth_rate(production, CellType::r,
                                           production.seed + 1, uid) !=
                    first.inherent_growth_rate;
            seed_changes_migration = seed_changes_migration ||
                sample_initial_migration_rate(production, CellType::r,
                                              production.seed + 1, uid) !=
                    first.migration_rate;
        }
    }

    assert(r_clamped_count > 0 && r_unclamped_count > 0);
    assert(seed_changes_growth && seed_changes_migration);
    const double r_growth_mean = r_growth_sum / sample_count;
    const double K_growth_mean = K_growth_sum / sample_count;
    assert(std::abs(r_growth_mean - truncated_normal_expected_mean(
        production.initial_r_growth_truncated_normal)) < 0.001);
    assert(std::abs(K_growth_mean - truncated_normal_expected_mean(
        production.initial_K_growth_truncated_normal)) < 0.002);

    // K has E[Beta(5,5)]*0.25 = 0.125. The very low-shape r distribution
    // converges more slowly and includes the configured <=0.5 replacement.
    assert(std::abs(K_migration_sum / sample_count - 0.125) < 0.001);
    assert(r_migration_sum / sample_count > 28.0);
    assert(r_migration_sum / sample_count < 32.0);
}
