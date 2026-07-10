#ifndef free_living_division_hpp
#define free_living_division_hpp

#include <cstdio>

#include "cell_store.hpp"
#include "cell_trace.hpp"
#include "visual_range.hpp"

void free_living_division(int i, double max_growth_rate_r, double max_growth_rate_K, CellStore &cell_array, VisualRange &Visual_range, CellRowBuffer cell_temp, int &cell_label, double &deltah, int utralsmall, double beta_distribution_alpha_for_normal_migration, double beta_distribution_beta_for_normal_migration, double migration_rate_K_mean, double uniup_K, double unilow_K, double sigmahatK, double muhatK, long &K_label, double beta_distribution_alpha, double beta_distribution_beta, double migration_rate_r_mean, double migration_rate_r_mean_quia, double beta_distribution_expected_for_normal_migration, CellTraceStore &cell_trace, CellTraceStore cell_trace_temp, long &cell_index, long &r_label, int Col, double K_formation_rate, FILE *fid2, int threads, CellTraceStore &cell_trace_ndcells, int &ndcells, CellRowBuffer &cell_array_ndcells, long rng_time_step);

#endif
