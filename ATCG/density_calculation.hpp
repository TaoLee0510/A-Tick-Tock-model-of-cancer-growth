#ifndef density_calculation_hpp
#define density_calculation_hpp

#include "cell_store.hpp"
#include "visual_range.hpp"

double density_calculation_from_position(int x1, int y1, int cell_stage, const VisualRange &Visual_range);
double density_calculation(int i, const VisualRange &Visual_range, const CellStore &cells);

#endif
