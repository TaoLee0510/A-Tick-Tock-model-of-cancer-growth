#ifndef CellDivision_hpp
#define CellDivision_hpp

#include "cell_store.hpp"
#include "visual_range.hpp"

void CellDivision(int i, double max_growth_rate_r, double max_growth_rate_K, CellStore &cell_array, VisualRange &Visual_range, CellRowBuffer &cell_temp, int &cell_label, double &deltah, int utralsmall, int Col, double deathjudge, int borderx, int bordery, long rng_time_step);

#endif
