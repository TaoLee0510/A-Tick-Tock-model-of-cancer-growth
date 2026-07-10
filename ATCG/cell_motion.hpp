#ifndef cell_motion_hpp
#define cell_motion_hpp

#include "cell_store.hpp"
#include "visual_range.hpp"

int migration_direction_dx(int direction);
int migration_direction_dy(int direction);
void write_big_cell_to_visual_range(VisualRange &Visual_range, int x1, int y1, long cell_id, long cell_stage, long cell_label);
void write_small_cell_to_visual_range(VisualRange &Visual_range, int x1, int y1, long cell_id, long cell_stage, long cell_label);
void move_cell_store(int i, CellStore &cells, VisualRange &Visual_range, int direction);

#endif
