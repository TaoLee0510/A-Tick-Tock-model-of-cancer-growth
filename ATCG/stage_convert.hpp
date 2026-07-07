//
//  stage_convert.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef stage_convert_hpp
#define stage_convert_hpp

#include <stdio.h>
#include <random>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_matrix.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "stateless_rng.hpp"
#include "cell_store.hpp"
#include "visual_range.hpp"
using namespace blitz;
namespace stage_convert_detail
{
inline bool within_visual_range(int x, int y, int Visual_range_x, int Visual_range_y)
{
    return x >= 100 && y >= 100 && x <= Visual_range_x + 100 && y <= Visual_range_y + 100;
}

inline int choose_large_stage_square(int x1, int y1, const VisualRange &Visual_range, long cell_rng_id, long rng_time_step, long &rng_event)
{
    int candidates[4] = {0};
    int count = 0;
    long loci_cor[9] = {0};
    int loci_direction = 0;
    for (int xx = x1 - 1; xx <= x1 + 1; ++xx)
    {
        for (int yy = y1 - 1; yy <= y1 + 1; ++yy)
        {
            loci_cor[loci_direction] = Visual_range(xx, yy, 1);
            ++loci_direction;
        }
    }

    if (loci_cor[0] == 0 && loci_cor[1] == 0 && loci_cor[3] == 0)
    {
        candidates[count++] = 1;
    }
    if (loci_cor[1] == 0 && loci_cor[2] == 0 && loci_cor[5] == 0)
    {
        candidates[count++] = 2;
    }
    if (loci_cor[5] == 0 && loci_cor[7] == 0 && loci_cor[8] == 0)
    {
        candidates[count++] = 3;
    }
    if (loci_cor[3] == 0 && loci_cor[6] == 0 && loci_cor[7] == 0)
    {
        candidates[count++] = 4;
    }

    if (count == 0)
    {
        return 0;
    }
    stateless_shuffle(candidates, candidates + count, cell_rng_id, rng_time_step, rng_event++);
    return candidates[0];
}

inline void square_loci(int stage_square, int x1, int y1, int xs[4], int ys[4])
{
    if (stage_square == 1)
    {
        xs[0] = x1 - 1; xs[1] = x1 - 1; xs[2] = x1;     xs[3] = x1;
        ys[0] = y1 - 1; ys[1] = y1;     ys[2] = y1;     ys[3] = y1 - 1;
    }
    else if (stage_square == 2)
    {
        xs[0] = x1 - 1; xs[1] = x1 - 1; xs[2] = x1;     xs[3] = x1;
        ys[0] = y1;     ys[1] = y1 + 1; ys[2] = y1 + 1; ys[3] = y1;
    }
    else if (stage_square == 3)
    {
        xs[0] = x1;     xs[1] = x1;     xs[2] = x1 + 1; xs[3] = x1 + 1;
        ys[0] = y1;     ys[1] = y1 + 1; ys[2] = y1 + 1; ys[3] = y1;
    }
    else
    {
        xs[0] = x1;     xs[1] = x1;     xs[2] = x1 + 1; xs[3] = x1 + 1;
        ys[0] = y1 - 1; ys[1] = y1;     ys[2] = y1;     ys[3] = y1 - 1;
    }
}

inline void apply_large_stage_square(int row, CellStore &cells, VisualRange &Visual_range, int stage_square, int x1, int y1)
{
    int xs[4] = {0};
    int ys[4] = {0};
    square_loci(stage_square, x1, y1, xs, ys);

    cells.stage()[row - 1] = 0;
    cells.x1()[row - 1] = xs[0];
    cells.x2()[row - 1] = xs[1];
    cells.x3()[row - 1] = xs[2];
    cells.x4()[row - 1] = xs[3];
    cells.y1()[row - 1] = ys[0];
    cells.y2()[row - 1] = ys[1];
    cells.y3()[row - 1] = ys[2];
    cells.y4()[row - 1] = ys[3];

    const long cell_id = (long)cells.id()[row - 1];
    const long cell_label = Visual_range(x1, y1, 4);
    for (int idx = 0; idx < 4; ++idx)
    {
        Visual_range(xs[idx], ys[idx], 1) = 1;
        Visual_range(xs[idx], ys[idx], 2) = cell_id;
        Visual_range(xs[idx], ys[idx], 3) = 0;
        Visual_range(xs[idx], ys[idx], 4) = cell_label;
    }
}

inline int choose_small_stage_direction(int x1, int y1, const VisualRange &Visual_range, long cell_rng_id, long rng_time_step, long &rng_event)
{
    int candidates[8] = {0};
    int count = 0;
    if (Visual_range(x1 - 1, y1 - 1, 1) == 0) { candidates[count++] = 1; }
    if (Visual_range(x1 - 1, y1,     1) == 0) { candidates[count++] = 2; }
    if (Visual_range(x1 - 1, y1 + 1, 1) == 0) { candidates[count++] = 3; }
    if (Visual_range(x1,     y1 + 1, 1) == 0) { candidates[count++] = 4; }
    if (Visual_range(x1 + 1, y1 + 1, 1) == 0) { candidates[count++] = 5; }
    if (Visual_range(x1 + 1, y1,     1) == 0) { candidates[count++] = 6; }
    if (Visual_range(x1 + 1, y1 - 1, 1) == 0) { candidates[count++] = 7; }
    if (Visual_range(x1,     y1 - 1, 1) == 0) { candidates[count++] = 8; }

    if (count == 0)
    {
        return 0;
    }
    stateless_shuffle(candidates, candidates + count, cell_rng_id, rng_time_step, rng_event++);
    return candidates[0];
}

inline void direction_locus(int direction, int x1, int y1, int &x2, int &y2)
{
    static const int dx[9] = {0, -1, -1, -1, 0, 1, 1, 1, 0};
    static const int dy[9] = {0, -1, 0, 1, 1, 1, 0, -1, -1};
    x2 = x1 + dx[direction];
    y2 = y1 + dy[direction];
}

inline void apply_small_stage_direction(int row, CellStore &cells, VisualRange &Visual_range, int &cell_label, int direction, int x1, int y1)
{
    int x2 = x1;
    int y2 = y1;
    direction_locus(direction, x1, y1, x2, y2);

    const int cell_label_1 = cell_label + 1;
    Visual_range(x2, y2, 1) = 1;
    Visual_range(x2, y2, 2) = (long)cells.id()[row - 1];
    Visual_range(x2, y2, 3) = 1;
    Visual_range(x2, y2, 4) = cell_label_1;

    cells.x1()[row - 1] = x2;
    cells.y1()[row - 1] = y2;
    cells.stage()[row - 1] = 1;
    Visual_range(x1, y1, 3) = 1;

    const int row_count = cells.rows();
    CellStore::Column &xs = cells.x1();
    CellStore::Column &ys = cells.y1();
    CellStore::Column &stages = cells.stage();
    for (int other = 1; other <= row_count; ++other)
    {
        if ((int)xs[other - 1] == x1 && (int)ys[other - 1] == y1)
        {
            stages[other - 1] = 1;
        }
    }
}
}

inline void stage_convert(int Visual_range_x, int Visual_range_y, CellStore &cells, VisualRange &Visual_range, int &cell_label,int utralsmall,long rng_time_step)
{
    const int row_count = cells.rows();
    for (int row = 1; row <= row_count; ++row)
    {
        long cell_rng_id = (long)cells.id()[row - 1];
        if (cell_rng_id == 0)
        {
            cell_rng_id = row;
        }
        long rng_event = 100;
        const int x1 = (int)cells.x1()[row - 1];
        const int y1 = (int)cells.y1()[row - 1];
        if (!stage_convert_detail::within_visual_range(x1, y1, Visual_range_x, Visual_range_y))
        {
            continue;
        }

        if ((int)cells.stage()[row - 1] == 1)
        {
            int stage_square = stage_convert_detail::choose_large_stage_square(x1, y1, Visual_range, cell_rng_id, rng_time_step, rng_event);
            if (stage_square != 0)
            {
                stage_convert_detail::apply_large_stage_square(row, cells, Visual_range, stage_square, x1, y1);
            }
        }
        else if ((int)cells.stage()[row - 1] == 2)
        {
            int direction = stage_convert_detail::choose_small_stage_direction(x1, y1, Visual_range, cell_rng_id, rng_time_step, rng_event);
            if (direction != 0)
            {
                stage_convert_detail::apply_small_stage_direction(row, cells, Visual_range, cell_label, direction, x1, y1);
            }
        }
    }

    if (utralsmall == 1)
    {
        for (int row = 1; row <= row_count; ++row)
        {
            long cell_rng_id = (long)cells.id()[row - 1];
            if (cell_rng_id == 0)
            {
                cell_rng_id = row;
            }
            long rng_event = 1000;
            const int x1 = (int)cells.x1()[row - 1];
            const int y1 = (int)cells.y1()[row - 1];
            if (!stage_convert_detail::within_visual_range(x1, y1, Visual_range_x, Visual_range_y))
            {
                continue;
            }
            if ((int)cells.stage()[row - 1] == 1)
            {
                int stage_square = stage_convert_detail::choose_large_stage_square(x1, y1, Visual_range, cell_rng_id, rng_time_step, rng_event);
                if (stage_square != 0)
                {
                    stage_convert_detail::apply_large_stage_square(row, cells, Visual_range, stage_square, x1, y1);
                }
            }
        }
    }
}
#endif /* stage_convert_hpp */
