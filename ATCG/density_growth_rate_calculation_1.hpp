#ifndef density_growth_rate_calculation_1_hpp
#define density_growth_rate_calculation_1_hpp

#include "cell_store.hpp"
#include "visual_range.hpp"
#include "common/density_growth_rule.hpp"

long unique_nonzero_count(long *values, int count);
bool is_r_density_label(long label, int N00, int N01);
DensityGrowthCounts density_growth_neighborhood_counts(int x1, int y1, int N00, int N01, const VisualRange &Visual_range);
void density_growth_rate_calculation_1(int Visual_range_x, int Visual_range_y, int N00, int N01, double r_limit, double K_limit, double lambda_r, double lambda_K, double alpha, double beta, double carrying_capacity_r, double carrying_capacity_K, double Cr, double CK, double death_time_range_r, double death_time_range_K, CellStore &cells, const VisualRange &Visual_range, long rng_time_step);

#endif
