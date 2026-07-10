#ifndef CellMigration_hpp
#define CellMigration_hpp

#include "cell_store.hpp"
#include "visual_range.hpp"

void CellMigration(int &DDM, int i, double &deltah, CellStore &cell_array, VisualRange &Visual_range, double &migration_judgement, double deathjudge, double beta_distribution_alpha_mig_time, double beta_distribution_beta_mig_time, int chemotaxis, double bunderD, int borderx, int bordery, double beta_distribution_alpha_for_normal_migration, double migration_rate_r_mean_quia, double beta_distribution_beta_for_normal_migration, long rng_time_step);

#endif
