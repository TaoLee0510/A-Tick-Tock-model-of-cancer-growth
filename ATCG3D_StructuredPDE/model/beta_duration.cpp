#include "model/beta_duration.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace atcg3d::structured_pde {
namespace {
double fraction(double x, double a, double b) {
    // Modified Lentz evaluation of DLMF 8.17.22-23, using symmetry in
    // beta_cdf to select the rapidly convergent end of the interval.
    const auto nonzero = [](double value) {
        return std::abs(value) < 1e-280 ? std::copysign(1e-280, value) : value;
    };
    double c = 1, d = 1 / nonzero(1 - (a + b) * x / (a + 1)), h = d;
    const auto update = [&](double coefficient) {
        d = 1 / nonzero(1 + coefficient * d);
        c = nonzero(1 + coefficient / c);
        const double factor = c * d;
        h *= factor;
        return factor;
    };
    for (int m = 1; m <= 10000; ++m) {
        update(m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m)));
        const double factor = update(-(a + m) * (a + b + m) * x /
            ((a + 2 * m) * (a + 2 * m + 1)));
        if (std::abs(factor - 1) < 3e-14) return h;
    }
    throw std::runtime_error("beta duration continued fraction did not converge");
}
}
double beta_cdf(double x, double a, double b) {
    if (!std::isfinite(x) || !(a > 0) || !(b > 0) || !std::isfinite(a + b))
        throw std::invalid_argument("invalid beta duration distribution");
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    const double factor = std::exp(std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) +
        a * std::log(x) + b * std::log1p(-x));
    const double result = x < (a + 1) / (a + b + 2) ? factor * fraction(x, a, b) / a :
        1 - factor * fraction(1 - x, b, a) / b;
    if (!std::isfinite(result)) throw std::runtime_error("nonfinite beta duration CDF");
    return std::clamp(result, 0.0, 1.0);
}
std::vector<double> beta_duration_kernel(double a, double b, double cycle,
                                         double width, double maximum) {
    if (!(width > 0) || !(cycle > 0) || !std::isfinite(cycle) ||
        !std::isfinite(width) || !std::isfinite(maximum) || !(maximum > width) || maximum < cycle ||
        maximum / width > 4096)
        throw std::invalid_argument("activation duration grid does not cover the full cycle");
    std::vector<double> result(std::size_t(std::ceil(maximum / width)) + 1);
    for (std::size_t i = 0; i + 1 < result.size(); ++i) {
        const double lo = std::min(1.0, i * width / cycle);
        const double hi = std::min(1.0, (i + 1) * width / cycle);
        if (hi <= lo) break;
        const double mass = std::max(0.0, beta_cdf(hi, a, b) - beta_cdf(lo, a, b));
        if (mass == 0) continue;
        const double moment = a / (a + b) *
            (beta_cdf(hi, a + 1, b) - beta_cdf(lo, a + 1, b));
        const double coordinate = std::clamp(cycle * moment / (mass * width), double(i), double(i + 1));
        result[i] += mass * (i + 1 - coordinate);
        result[i + 1] += mass * (coordinate - i);
    }
    const double total = std::accumulate(result.begin(), result.end(), 0.0);
    for (auto& value : result) value /= total;
    return result;
}
} // namespace atcg3d::structured_pde
