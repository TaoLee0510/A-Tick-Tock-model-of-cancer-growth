//
//  cell_motion.hpp
//  ATCG
//

#include "cell_motion.hpp"

#include <blitz/blitz.h>
#include <blitz/array.h>
#include "cell_columns.hpp"
#include "cell_store.hpp"
#include "visual_range.hpp"

using namespace blitz;

int migration_direction_dx(int direction)
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

int migration_direction_dy(int direction)
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

void write_big_cell_to_visual_range(VisualRange &Visual_range, int x1, int y1, long cell_id, long cell_stage, long cell_label)
{
    Visual_range.write_square(x1, y1, cell_id, cell_stage, cell_label);
}

void write_small_cell_to_visual_range(VisualRange &Visual_range, int x1, int y1, long cell_id, long cell_stage, long cell_label)
{
    Visual_range.write_site(x1, y1, cell_id, cell_stage, cell_label);
}

void move_cell_store(int i, CellStore &cells, VisualRange &Visual_range, int direction)
{
    int row = i - 1;
    auto &x1_values = cells.x1();
    auto &x2_values = cells.x2();
    auto &x3_values = cells.x3();
    auto &x4_values = cells.x4();
    auto &y1_values = cells.y1();
    auto &y2_values = cells.y2();
    auto &y3_values = cells.y3();
    auto &y4_values = cells.y4();
    auto &stages = cells.stage();
    auto &ids = cells.id();
    auto &migration_directions = cells.migration_direction();
    auto &migration_elapsed = cells.migration_elapsed();
    int dx = migration_direction_dx(direction);
    int dy = migration_direction_dy(direction);
    int x1 = (int)x1_values[row];
    int y1 = (int)y1_values[row];
    int cell_stage = (int)stages[row];
    long visual_stage = Visual_range.stage(x1,y1);
    long cell_label = Visual_range.cell_label(x1,y1);
    long cell_id = (long)ids[row];

    if (cell_stage == 0)
    {
        Visual_range.clear_square(x1, y1);
        x1_values[row] += dx;
        x2_values[row] += dx;
        x3_values[row] += dx;
        x4_values[row] += dx;
        y1_values[row] += dy;
        y2_values[row] += dy;
        y3_values[row] += dy;
        y4_values[row] += dy;

        write_big_cell_to_visual_range(Visual_range, (int)x1_values[row], (int)y1_values[row], cell_id, visual_stage, cell_label);
    }
    else
    {
        Visual_range.clear_site(x1, y1);
        x1_values[row] += dx;
        y1_values[row] += dy;
        write_small_cell_to_visual_range(Visual_range, (int)x1_values[row], (int)y1_values[row], cell_id, visual_stage, cell_label);
    }

    migration_directions[row]=direction;
    migration_elapsed[row]=0;
}

