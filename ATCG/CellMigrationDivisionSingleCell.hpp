#ifndef CellMigrationDivisionSingleCell_hpp
#define CellMigrationDivisionSingleCell_hpp

#include <cstdio>

#include "cell_store.hpp"
#include "cell_trace.hpp"
#include "visual_range.hpp"

void CellMigrationDivisionSingleCell(int &DDM, int i, double &deltah, CellStore &cell_array, VisualRange &Visual_range, double &migration_judgement, double max_growth_rate_r, double max_growth_rate_K, CellRowBuffer cell_temp, int &cell_label, int utralsmall, double beta_distribution_alpha_for_normal_migration, double beta_distribution_beta_for_normal_migration, double migration_rate_K_mean, double uniup_K, double unilow_K, double sigmahatK, double muhatK, long &K_label, double beta_distribution_alpha, double beta_distribution_beta, double migration_rate_r_mean, double migration_rate_r_mean_quia, double beta_distribution_expected_for_normal_migration, CellTraceStore &cell_trace, CellTraceStore cell_trace_temp, long &cell_index, long &r_label, int Col, double K_formation_rate, double deathjudge, double beta_distribution_alpha_mig_time, double beta_distribution_beta_mig_time, int chemotaxis, double bunderD, int borderx, int bordery, FILE *fid2, int threads, long rng_time_step);

#endif
