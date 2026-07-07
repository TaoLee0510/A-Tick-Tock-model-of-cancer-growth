//
//  cell_motion.hpp
//  ATCG
//

#ifndef cell_motion_hpp
#define cell_motion_hpp

#include <blitz/blitz.h>
#include <blitz/array.h>
#include "cell_columns.hpp"
#include "cell_store.hpp"
#include "visual_range.hpp"

using namespace blitz;

inline int migration_direction_dx(int direction)
{
    switch (direction)
    {
        case 1:
        case 2:
        case 3:
            return -1;
        case 5:
        case 6:
        case 7:
            return 1;
        default:
            return 0;
    }
}

inline int migration_direction_dy(int direction)
{
    switch (direction)
    {
        case 1:
        case 7:
        case 8:
            return -1;
        case 3:
        case 4:
        case 5:
            return 1;
        default:
            return 0;
    }
}

inline void write_big_cell_to_visual_range(VisualRange &Visual_range, int x1, int y1, long cell_id, long cell_stage, long cell_label)
{
    Visual_range.write_square(x1, y1, cell_id, cell_stage, cell_label);
}

inline void write_small_cell_to_visual_range(VisualRange &Visual_range, int x1, int y1, long cell_id, long cell_stage, long cell_label)
{
    Visual_range.write_site(x1, y1, cell_id, cell_stage, cell_label);
}

inline void move_cell_store(int i, CellStore &cells, VisualRange &Visual_range, int direction)
{
    int dx = migration_direction_dx(direction);
    int dy = migration_direction_dy(direction);
    int x1 = (int)cells(i,cell_col::kX1);
    int y1 = (int)cells(i,cell_col::kY1);
    int cell_stage = (int)cells(i,cell_col::kStage);
    long visual_stage = Visual_range.stage(x1,y1);
    long cell_label = Visual_range.cell_label(x1,y1);
    long cell_id = (long)cells(i,cell_col::kId);

    if (cell_stage == 0)
    {
        Visual_range.clear_square(x1, y1);
        for (int x_col = cell_col::kX1; x_col <= cell_col::kX4; ++x_col)
        {
            cells(i,x_col)=cells(i,x_col)+dx;
        }
        for (int y_col = cell_col::kY1; y_col <= cell_col::kY4; ++y_col)
        {
            cells(i,y_col)=cells(i,y_col)+dy;
        }

        write_big_cell_to_visual_range(Visual_range, (int)cells(i,cell_col::kX1), (int)cells(i,cell_col::kY1), cell_id, visual_stage, cell_label);
    }
    else
    {
        Visual_range.clear_site(x1, y1);
        cells(i,cell_col::kX1)=cells(i,cell_col::kX1)+dx;
        cells(i,cell_col::kY1)=cells(i,cell_col::kY1)+dy;
        write_small_cell_to_visual_range(Visual_range, (int)cells(i,cell_col::kX1), (int)cells(i,cell_col::kY1), cell_id, visual_stage, cell_label);
    }

    cells(i,cell_col::kMigrationDirection)=direction;
    cells(i,cell_col::kMigrationElapsed)=0;
}

#endif /* cell_motion_hpp */
