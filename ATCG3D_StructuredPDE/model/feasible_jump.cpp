#include "model/feasible_jump.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace atcg3d::structured_pde {

std::vector<double> uniform_feasible_jump_probabilities(std::span<const double> availability) {
    if (availability.empty() || availability.size() > 26 ||
        std::any_of(availability.begin(), availability.end(), [](double value) {
            return !std::isfinite(value) || value < 0.0 || value > 1.0;
        })) {
        throw std::invalid_argument("feasible jump requires one to 26 probabilities in [0,1]");
    }
    if (std::all_of(availability.begin(), availability.end(), [](double value) { return value == 1.0; })) {
        return std::vector<double>(availability.size(), 1.0 / static_cast<double>(availability.size()));
    }
    std::vector<double> result(availability.size(), 0.0);
    std::vector<double> polynomial(availability.size(), 0.0);
    for (std::size_t selected = 0; selected < availability.size(); ++selected) {
        if (availability[selected] == 0.0) continue;
        std::fill(polynomial.begin(), polynomial.end(), 0.0);
        polynomial[0] = 1.0;
        std::size_t degree = 0;
        for (std::size_t other = 0; other < availability.size(); ++other) {
            if (other == selected) continue;
            const double chance = availability[other];
            for (std::size_t power = degree + 1; power > 0; --power) {
                polynomial[power] = polynomial[power] * (1.0 - chance) +
                    polynomial[power - 1] * chance;
            }
            polynomial[0] *= 1.0 - chance;
            ++degree;
        }
        double integral = 0.0;
        for (std::size_t power = 0; power <= degree; ++power) {
            integral += polynomial[power] / static_cast<double>(power + 1);
        }
        result[selected] = availability[selected] * integral;
    }
    return result;
}

}  // namespace atcg3d::structured_pde
