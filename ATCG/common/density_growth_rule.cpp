#include "common/density_growth_rule.hpp"

double calculate_density_growth_rate(int cell_type,
                                     double inherent_growth_rate,
                                     const DensityGrowthCounts &counts,
                                     double r_limit,
                                     double K_limit,
                                     double alpha,
                                     double beta,
                                     double carrying_capacity_r,
                                     double carrying_capacity_K)
{
    return calculate_density_growth_rate_continuous(
        cell_type, inherent_growth_rate,
        static_cast<double>(counts.rc), static_cast<double>(counts.kc),
        static_cast<double>(counts.cells_number), r_limit, K_limit, alpha, beta,
        carrying_capacity_r, carrying_capacity_K);
}

double calculate_density_growth_rate_continuous(int cell_type,
                                                double inherent_growth_rate,
                                                double r_count,
                                                double K_count,
                                                double total_count,
                                                double r_limit,
                                                double K_limit,
                                                double alpha,
                                                double beta,
                                                double carrying_capacity_r,
                                                double carrying_capacity_K)
{
    switch (cell_type)
    {
        case 1:
        {
            if (total_count < r_limit)
            {
                return inherent_growth_rate;
            }
            double crowding = r_count + K_count + alpha * K_count - r_limit;
            return inherent_growth_rate - (inherent_growth_rate * 2 * crowding) / carrying_capacity_r;
        }
        case 2:
        {
            if (total_count < K_limit)
            {
                return inherent_growth_rate;
            }
            double crowding = beta * r_count + r_count + K_count - K_limit;
            return inherent_growth_rate - (inherent_growth_rate * 2 * crowding) / carrying_capacity_K;
        }
        default:
        {
            return inherent_growth_rate;
        }
    }
}
