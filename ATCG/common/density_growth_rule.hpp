#ifndef ATCG_COMMON_DENSITY_GROWTH_RULE_HPP
#define ATCG_COMMON_DENSITY_GROWTH_RULE_HPP

struct DensityGrowthCounts
{
    long rc;
    long kc;
    long cells_number;
};

double calculate_density_growth_rate(int cell_type,
                                     double inherent_growth_rate,
                                     const DensityGrowthCounts &counts,
                                     double r_limit,
                                     double K_limit,
                                     double alpha,
                                     double beta,
                                     double carrying_capacity_r,
                                     double carrying_capacity_K);

// Geometry-independent continuous-count form used when a 3D environmental
// field (for example a vascular relief gradient) changes effective crowding.
// The legacy integer-count API above remains unchanged and delegates here.
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
                                                double carrying_capacity_K);

#endif
