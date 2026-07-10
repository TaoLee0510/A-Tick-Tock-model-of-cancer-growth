#ifndef migrate_activation_hpp
#define migrate_activation_hpp

#include "cell_store.hpp"
#include "visual_range.hpp"

void migrate_activation(CellStore &cells, double bunderD, const VisualRange &Visual_range, double migration_time_range, double migration_rate_r_mean_quia, double beta_distribution_alpha_for_normal_migration, double beta_distribution_beta_for_normal_migration, double beta_distribution_alpha_mig_time, double beta_distribution_beta_mig_time, int DDM, long rng_time_step);

#endif
