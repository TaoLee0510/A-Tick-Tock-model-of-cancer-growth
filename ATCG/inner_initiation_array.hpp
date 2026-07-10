#ifndef inner_initiation_array_hpp
#define inner_initiation_array_hpp

#include "cell_store.hpp"
#include "visual_range.hpp"

void fill_inner_initiation_cells(int N0, int N01, int R0, int Visual_range_x, int Visual_range_y, CellStore &cell_array_inner, const VisualRange &Visual_range, double uniup_r1, double unilow_r1, double sigmahatr, double muhatr, double uniup_K1, double unilow_K1, double sigmahatK, double muhatK, int N0r1, int N0K1, double *migration_rate_r1, double *migration_rate_K1, int Col);
CellStore inner_initiation_cell_store(int N0, int N01, int R0, int Visual_range_x, int Visual_range_y, const VisualRange &Visual_range, double uniup_r1, double unilow_r1, double sigmahatr, double muhatr, double uniup_K1, double unilow_K1, double sigmahatK, double muhatK, int N0r1, int N0K1, double *migration_rate_r1, double *migration_rate_K1, int Col);

#endif
