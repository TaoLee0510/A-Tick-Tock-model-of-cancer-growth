//
//  random_migration.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef random_migration_hpp
#define random_migration_hpp

#include <stdio.h>
#include <random>
#include <memory>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <ctime>
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
#define BZ_THREADSAFE
#define BZ_THREADSAFE_USE_OPENMP
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "cell_columns.hpp"
#include "cell_store.hpp"
#include "visual_range.hpp"
#include "deltah_calculation.hpp"
#include "cell_motion.hpp"
#include "stateless_rng.hpp"
using namespace blitz;

inline int select_random_migration_direction(int x1, int y1, int cell_stage, const VisualRange &Visual_range, long cell_rng_id, long rng_time_step, long rng_event)
{
    int direction[8]={0};
    auto is_empty = [&](int dx, int dy) {
        return Visual_range(x1 + dx, y1 + dy, 1) == 0;
    };
    if (cell_stage==0)
    {
        if (is_empty(0, -1) && is_empty(-1, -1) && is_empty(-1, 0))
        {
            direction[0]=1;
        }
        if (is_empty(-1, 0) && is_empty(-1, 1))
        {
            direction[1]=2;
        }
        if (is_empty(-1, 1) && is_empty(-1, 2) && is_empty(0, 2))
        {
            direction[2]=3;
        }
        if (is_empty(0, 2) && is_empty(1, 2))
        {
            direction[3]=4;
        }
        if (is_empty(1, 2) && is_empty(2, 2) && is_empty(2, 1))
        {
            direction[4]=5;
        }
        if (is_empty(2, 1) && is_empty(2, 0))
        {
            direction[5]=6;
        }
        if (is_empty(2, 0) && is_empty(2, -1) && is_empty(1, -1))
        {
            direction[6]=7;
        }
        if (is_empty(1, -1) && is_empty(0, -1))
        {
            direction[7]=8;
        }
    }
    else //small stage
    {
        if (is_empty(-1, -1))
        {
            direction[0]=1;
        }
        if (is_empty(-1, 0))
        {
            direction[1]=2;
        }
        if (is_empty(-1, 1))
        {
            direction[2]=3;
        }
        if (is_empty(0, 1))
        {
            direction[3]=4;
        }
        if (is_empty(1, 1))
        {
            direction[4]=5;
        }
        if (is_empty(1, 0))
        {
            direction[5]=6;
        }
        if (is_empty(1, -1))
        {
            direction[6]=7;
        }
        if (is_empty(0, -1))
        {
            direction[7]=8;
        }
    }

    int candidates[8]={0};
    int candidate_count=0;
    for (int loci=0; loci<8; loci++)
    {
        if (direction[loci]!=0)
        {
            candidates[candidate_count]=direction[loci];
            candidate_count++;
        }
    }
    if (candidate_count==0)
    {
        return 0;
    }

    stateless_shuffle(candidates, candidates + candidate_count, cell_rng_id, rng_time_step, rng_event);
    return candidates[0];
}

inline void random_migration(int i, double deltah, CellStore &cells, VisualRange &Visual_range, double &migration_judgement, long rng_time_step, long rng_event_base)
{
    (void)deltah;

    int row = i - 1;
    long cell_rng_id = (long)cells.id()[row];
    if (cell_rng_id == 0)
    {
        cell_rng_id = i;
    }
    int cell_stage=(int)cells.stage()[row];
    long rng_event = rng_event_base + ((long)cell_stage * 100);
    int x1=(int)cells.x1()[row];
    int y1=(int)cells.y1()[row];
    int order = select_random_migration_direction(x1, y1, cell_stage, Visual_range, cell_rng_id, rng_time_step, rng_event);
    if (order!=0)
    {
        move_cell_store(i, cells, Visual_range, order);
        cells.migration_elapsed()[row]=0;
        cells.migration_follow_flag()[row]=1;
    }
    migration_judgement=migration_judgement+0.0001;
}

#endif /* random_migration_hpp */
