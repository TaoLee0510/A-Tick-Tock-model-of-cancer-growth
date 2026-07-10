#ifndef outer_initiation_array_low_density_hpp
#define outer_initiation_array_low_density_hpp

#include "cell_store.hpp"
#include "int_grid.hpp"

CellStore outer_initiation_low_density_cell_store(int N0, int Visual_range_x, int Visual_range_y, const IntGrid &A, double uniup_r, double unilow_r, double sigmahatr, double muhatr, double uniup_K, double unilow_K, double sigmahatK, double muhatK, int N0r, int N0K, double *migration_rate_r, double *migration_rate_K);

#endif
