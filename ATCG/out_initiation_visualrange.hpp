#ifndef out_initiation_visualrange_hpp
#define out_initiation_visualrange_hpp

#include "cell_store.hpp"
#include "visual_range.hpp"

void set_outer_visual_range_row(VisualRange &Visual_range, int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4, int cell_array_index, int cell_array_stage, int cell_label);
VisualRange outer_initiation_visualrange(const CellStore &cells, int N0, int Vx, int Vy, int &cell_label);

#endif
