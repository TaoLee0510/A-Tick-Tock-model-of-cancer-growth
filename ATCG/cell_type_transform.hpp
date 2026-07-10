#ifndef cell_type_transform_hpp
#define cell_type_transform_hpp

#include "cell_store.hpp"
#include "visual_range.hpp"

void cell_type_transform(CellRowBuffer &cell_temp, double beta_distribution_alpha_for_normal_migration, double beta_distribution_beta_for_normal_migration, double migration_rate_K_mean, double uniup_K, double unilow_K, double sigmahatK, double muhatK, long &K_label, int i, const VisualRange &Visual_range, CellStore &cell_array, double beta_distribution_alpha, double beta_distribution_beta, double migration_rate_r_mean, double migration_rate_r_mean_quia, double beta_distribution_expected_for_normal_migration, long &r_label, double K_formation_rate, long cell_rng_id, long rng_time_step, long &rng_event);

#endif
