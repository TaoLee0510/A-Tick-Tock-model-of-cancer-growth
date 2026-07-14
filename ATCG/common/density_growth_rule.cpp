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
    switch (cell_type)
    {
        case 1:
        {
            if (counts.cells_number < r_limit)
            {
                return inherent_growth_rate;
            }
            double crowding = counts.rc + counts.kc + alpha * counts.kc - r_limit;
            return inherent_growth_rate - (inherent_growth_rate * 2 * crowding) / carrying_capacity_r;
        }
        case 2:
        {
            if (counts.cells_number < K_limit)
            {
                return inherent_growth_rate;
            }
            double crowding = beta * counts.rc + counts.rc + counts.kc - K_limit;
            return inherent_growth_rate - (inherent_growth_rate * 2 * crowding) / carrying_capacity_K;
        }
        default:
        {
            return inherent_growth_rate;
        }
    }
}
