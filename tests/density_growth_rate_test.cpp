#include "density_growth_rate_calculation_1.hpp"

#include <cmath>
#include <iostream>

namespace
{
bool nearly_equal(double lhs, double rhs)
{
    return std::abs(lhs - rhs) < 1e-12;
}

bool expect_near(const char *name, double actual, double expected)
{
    if (nearly_equal(actual, expected))
    {
        return true;
    }
    std::cerr << name << ": expected " << expected << ", got " << actual << '\n';
    return false;
}
}

int main()
{
    const double r_limit = 5.0;
    const double K_limit = 5.0;
    const double alpha = 0.5;
    const double beta = 0.5;
    const double carrying_capacity_r = 20.0;
    const double carrying_capacity_K = 10.0;

    bool passed = true;
    DensityGrowthCounts below_limit{2, 2, 4};
    passed &= expect_near("r below limit", calculate_density_growth_rate(1, 1.2, below_limit, r_limit, K_limit, alpha, beta, carrying_capacity_r, carrying_capacity_K), 1.2);
    passed &= expect_near("K below limit", calculate_density_growth_rate(2, 0.8, below_limit, r_limit, K_limit, alpha, beta, carrying_capacity_r, carrying_capacity_K), 0.8);

    DensityGrowthCounts crowded{6, 2, 8};
    passed &= expect_near("r crowded", calculate_density_growth_rate(1, 1.2, crowded, r_limit, K_limit, alpha, beta, carrying_capacity_r, carrying_capacity_K), 0.72);
    passed &= expect_near("K crowded", calculate_density_growth_rate(2, 0.8, crowded, r_limit, K_limit, alpha, beta, carrying_capacity_r, carrying_capacity_K), -0.16);

    passed &= expect_near("unknown type", calculate_density_growth_rate(99, 1.0, crowded, r_limit, K_limit, alpha, beta, carrying_capacity_r, carrying_capacity_K), 1.0);
    return passed ? 0 : 1;
}
