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

inline void write_big_cell_to_visual_range(Array<long, 3> &Visual_range, int x1, int y1, long cell_id, long cell_stage, long cell_label)
{
    Visual_range(Range(x1,x1+1),Range(y1,y1+1),1)=1;
    Visual_range(Range(x1,x1+1),Range(y1,y1+1),2)=cell_id;
    Visual_range(Range(x1,x1+1),Range(y1,y1+1),3)=cell_stage;
    Visual_range(Range(x1,x1+1),Range(y1,y1+1),4)=cell_label;
}

inline void write_small_cell_to_visual_range(Array<long, 3> &Visual_range, int x1, int y1, long cell_id, long cell_stage, long cell_label)
{
    Visual_range(x1,y1,1)=1;
    Visual_range(x1,y1,2)=cell_id;
    Visual_range(x1,y1,3)=cell_stage;
    Visual_range(x1,y1,4)=cell_label;
}

inline void move_cell_store(int i, CellStore &cells, Array<long, 3> &Visual_range, int direction)
{
    Range all = Range::all();
    int dx = migration_direction_dx(direction);
    int dy = migration_direction_dy(direction);
    int x1 = (int)cells(i,cell_col::kX1);
    int y1 = (int)cells(i,cell_col::kY1);
    int cell_stage = (int)cells(i,cell_col::kStage);
    long visual_stage = Visual_range(x1,y1,3);
    long cell_label = Visual_range(x1,y1,4);
    long cell_id = (long)cells(i,cell_col::kId);

    if (cell_stage == 0)
    {
        Visual_range(Range(x1,x1+1),Range(y1,y1+1),all)=0;
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
        Visual_range(x1,y1,all)=0;
        cells(i,cell_col::kX1)=cells(i,cell_col::kX1)+dx;
        cells(i,cell_col::kY1)=cells(i,cell_col::kY1)+dy;
        write_small_cell_to_visual_range(Visual_range, (int)cells(i,cell_col::kX1), (int)cells(i,cell_col::kY1), cell_id, visual_stage, cell_label);
    }

    cells(i,cell_col::kMigrationDirection)=direction;
    cells(i,cell_col::kMigrationElapsed)=0;
}

#endif /* cell_motion_hpp */
