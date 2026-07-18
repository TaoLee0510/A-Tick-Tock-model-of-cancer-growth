#include "rules/initial_rates.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include "core/stateless_rng.hpp"

namespace atcg3d {
namespace {

constexpr std::uint64_t kInitialGrowthDomain = 0x4154434733475257ULL;
constexpr std::uint64_t kInitialMigrationDomain = 0x41544347334d4947ULL;
constexpr double kInverseSqrtTwo = 0.70710678118654752440;
constexpr double kTwoPi = 6.28318530717958647693;

double standard_normal_cdf(double value) {
    return 0.5 * std::erfc(-value * kInverseSqrtTwo);
}

// Peter J. Acklam's inverse-normal rational approximation. The configured
// legacy truncation intervals are central, but the tail branches keep this
// reusable for later calibrated profiles without an external statistics library.
double inverse_standard_normal_cdf(double probability) {
    if (!(probability > 0.0 && probability < 1.0)) {
        throw std::invalid_argument("normal quantile probability must be in (0,1)");
    }
    constexpr double a1 = -3.969683028665376e+01;
    constexpr double a2 = 2.209460984245205e+02;
    constexpr double a3 = -2.759285104469687e+02;
    constexpr double a4 = 1.383577518672690e+02;
    constexpr double a5 = -3.066479806614716e+01;
    constexpr double a6 = 2.506628277459239e+00;
    constexpr double b1 = -5.447609879822406e+01;
    constexpr double b2 = 1.615858368580409e+02;
    constexpr double b3 = -1.556989798598866e+02;
    constexpr double b4 = 6.680131188771972e+01;
    constexpr double b5 = -1.328068155288572e+01;
    constexpr double c1 = -7.784894002430293e-03;
    constexpr double c2 = -3.223964580411365e-01;
    constexpr double c3 = -2.400758277161838e+00;
    constexpr double c4 = -2.549732539343734e+00;
    constexpr double c5 = 4.374664141464968e+00;
    constexpr double c6 = 2.938163982698783e+00;
    constexpr double d1 = 7.784695709041462e-03;
    constexpr double d2 = 3.224671290700398e-01;
    constexpr double d3 = 2.445134137142996e+00;
    constexpr double d4 = 3.754408661907416e+00;
    constexpr double lower_tail = 0.02425;
    constexpr double upper_tail = 1.0 - lower_tail;

    if (probability < lower_tail) {
        const double q = std::sqrt(-2.0 * std::log(probability));
        return (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
               ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
    }
    if (probability > upper_tail) {
        const double q = std::sqrt(-2.0 * std::log1p(-probability));
        return -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
               ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
    }
    const double q = probability - 0.5;
    const double r = q * q;
    return (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
           (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0);
}

double sample_truncated_normal(const TruncatedNormalRateConfig& parameters,
                               std::uint64_t seed,
                               CellUid uid,
                               std::uint64_t type_domain) {
    const double lower_z = (parameters.minimum - parameters.mean) /
                           parameters.standard_deviation;
    const double upper_z = (parameters.maximum - parameters.mean) /
                           parameters.standard_deviation;
    const double lower_probability = standard_normal_cdf(lower_z);
    const double upper_probability = standard_normal_cdf(upper_z);
    if (!(upper_probability > lower_probability)) {
        throw std::invalid_argument("truncated-normal interval has no representable probability mass");
    }
    const double uniform = rng_unit(seed, uid, type_domain, 0, 0);
    double probability = lower_probability +
                         (upper_probability - lower_probability) * uniform;
    probability = std::clamp(
        probability, std::numeric_limits<double>::min(),
        std::nextafter(1.0, 0.0));
    const double sampled = parameters.mean + parameters.standard_deviation *
        inverse_standard_normal_cdf(probability);
    return std::clamp(sampled, parameters.minimum, parameters.maximum);
}

class StatelessDrawStream {
public:
    StatelessDrawStream(std::uint64_t seed, CellUid uid, std::uint64_t domain,
                        std::uint64_t event_sequence)
        : seed_(seed), uid_(uid), domain_(domain), event_sequence_(event_sequence) {}

    double unit_open() {
        const double value = rng_unit(seed_, uid_, domain_, event_sequence_, draw_++);
        return std::clamp(value, std::numeric_limits<double>::min(),
                          std::nextafter(1.0, 0.0));
    }

private:
    std::uint64_t seed_{};
    CellUid uid_{};
    std::uint64_t domain_{};
    std::uint64_t event_sequence_{};
    std::uint64_t draw_{};
};

double sample_standard_normal(StatelessDrawStream& stream) {
    const double radius = std::sqrt(-2.0 * std::log(stream.unit_open()));
    return radius * std::cos(kTwoPi * stream.unit_open());
}

double sample_log_gamma(double shape, StatelessDrawStream& stream) {
    if (!(shape > 0.0)) {
        throw std::invalid_argument("gamma shape must be positive");
    }
    if (shape < 1.0) {
        const double augmentation = stream.unit_open();
        return sample_log_gamma(shape + 1.0, stream) +
               std::log(augmentation) / shape;
    }

    const double d = shape - 1.0 / 3.0;
    const double c = 1.0 / std::sqrt(9.0 * d);
    for (int attempt = 0; attempt < 1024; ++attempt) {
        const double normal = sample_standard_normal(stream);
        const double base = 1.0 + c * normal;
        if (base <= 0.0) {
            continue;
        }
        const double value = base * base * base;
        const double uniform = stream.unit_open();
        const double normal_squared = normal * normal;
        if (uniform < 1.0 - 0.0331 * normal_squared * normal_squared ||
            std::log(uniform) < 0.5 * normal_squared +
                d * (1.0 - value + std::log(value))) {
            return std::log(d) + std::log(value);
        }
    }
    throw std::runtime_error("deterministic gamma sampler exceeded rejection limit");
}

double sample_beta_fraction(double alpha,
                            double beta,
                            std::uint64_t seed,
                            CellUid uid,
                            std::uint64_t domain,
                            std::uint64_t event_sequence = 0) {
    StatelessDrawStream alpha_stream(seed, uid, domain, event_sequence);
    StatelessDrawStream beta_stream(seed, uid, domain + 1U, event_sequence);
    const double log_alpha_gamma = sample_log_gamma(alpha, alpha_stream);
    const double log_beta_gamma = sample_log_gamma(beta, beta_stream);
    const double difference = log_alpha_gamma - log_beta_gamma;
    if (difference >= 0.0) {
        return 1.0 / (1.0 + std::exp(-difference));
    }
    const double exponential = std::exp(difference);
    return exponential / (1.0 + exponential);
}

double sample_beta_rate(const BetaRateConfig& parameters,
                        std::uint64_t seed,
                        CellUid uid,
                        std::uint64_t domain) {
    double rate = parameters.scale * sample_beta_fraction(
        parameters.alpha, parameters.beta, seed, uid, domain);
    if (parameters.lower_clamp_enabled && rate <= parameters.lower_clamp_threshold) {
        rate = parameters.lower_clamp_value;
    }
    return rate;
}

std::uint64_t type_offset(CellType type) {
    return type == CellType::r ? 0U : 16U;
}

}  // namespace

double sample_initial_growth_rate(const Model3DConfig& config,
                                  CellType type,
                                  std::uint64_t seed,
                                  CellUid uid) {
    if (config.initial_growth_rate_model == "fixed") {
        return type == CellType::r ? config.initial_r_growth_rate
                                   : config.initial_K_growth_rate;
    }
    if (config.initial_growth_rate_model == "legacy_truncated_normal_v1") {
        const TruncatedNormalRateConfig& parameters = type == CellType::r
            ? config.initial_r_growth_truncated_normal
            : config.initial_K_growth_truncated_normal;
        return sample_truncated_normal(
            parameters, seed, uid, kInitialGrowthDomain + type_offset(type));
    }
    throw std::invalid_argument("unsupported initial growth-rate model");
}

double sample_initial_migration_rate(const Model3DConfig& config,
                                     CellType type,
                                     std::uint64_t seed,
                                     CellUid uid) {
    if (type == CellType::r) {
        if (config.activated_r_migration_rate_model != "beta") {
            throw std::invalid_argument(
                "unsupported activated r migration-rate model");
        }
        return sample_beta_rate(config.activated_r_migration_beta, seed, uid,
                                kInitialMigrationDomain + type_offset(type));
    }
    if (config.initial_K_migration_rate_model == "fixed") {
        return config.initial_K_migration_rate;
    }
    if (config.initial_K_migration_rate_model == "legacy_beta_v1") {
        return sample_beta_rate(config.initial_K_migration_beta, seed, uid,
                                kInitialMigrationDomain + type_offset(type));
    }
    throw std::invalid_argument("unsupported initial K migration-rate model");
}

InitialCellRates3D sample_initial_cell_rates(const Model3DConfig& config,
                                             CellType type,
                                             std::uint64_t seed,
                                             CellUid uid) {
    return {sample_initial_growth_rate(config, type, seed, uid),
            sample_initial_migration_rate(config, type, seed, uid)};
}

double sample_stateless_beta_fraction(double alpha,
                                      double beta,
                                      std::uint64_t seed,
                                      CellUid uid,
                                      std::uint64_t event_domain,
                                      std::uint64_t event_sequence) {
    return sample_beta_fraction(alpha, beta, seed, uid, event_domain,
                                event_sequence);
}

}  // namespace atcg3d
