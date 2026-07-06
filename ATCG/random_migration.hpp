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
#include "deltah_calculation.hpp"
#include "cell_motion.hpp"
#include "stateless_rng.hpp"
using namespace blitz;

inline int select_random_migration_direction(int x1, int y1, int cell_stage, Array<long, 3> &Visual_range, Array<int,2> &cor_big, Array<int, 2> &cor_small, long cell_rng_id, long rng_time_step, long rng_event)
{
    Range all = Range::all();
    int direction[8]={0};
    if (cell_stage==0)
    {
        cor_big.resize(4,4);
        cor_big=0;
        cor_big(all,all)=Visual_range(Range(x1-1,x1+2),Range(y1-1,y1+2),1);
        if (cor_big(2,1)==0 && cor_big(1,1)==0 && cor_big(1,2)==0)
        {
            direction[0]=1;
        }
        if (cor_big(1,2)==0 && cor_big(1,3)==0)
        {
            direction[1]=2;
        }
        if (cor_big(1,3)==0 && cor_big(1,4)==0 && cor_big(2,4)==0)
        {
            direction[2]=3;
        }
        if (cor_big(2,4)==0 && cor_big(3,4)==0)
        {
            direction[3]=4;
        }
        if (cor_big(3,4)==0 && cor_big(4,4)==0 && cor_big(4,3)==0)
        {
            direction[4]=5;
        }
        if (cor_big(4,3)==0 && cor_big(4,2)==0)
        {
            direction[5]=6;
        }
        if (cor_big(4,2)==0 && cor_big(4,1)==0 && cor_big(3,1)==0)
        {
            direction[6]=7;
        }
        if (cor_big(3,1)==0 && cor_big(2,1)==0)
        {
            direction[7]=8;
        }
    }
    else //small stage
    {
        cor_small.resize(3, 3);
        cor_small=0;
        cor_small(all,all)=Visual_range(Range(x1-1,x1+1),Range(y1-1,y1+1),1);
        if (cor_small(1,1)==0)
        {
            direction[0]=1;
        }
        if (cor_small(1,2)==0)
        {
            direction[1]=2;
        }
        if (cor_small(1,3)==0)
        {
            direction[2]=3;
        }
        if (cor_small(2,3)==0)
        {
            direction[3]=4;
        }
        if (cor_small(3,3)==0)
        {
            direction[4]=5;
        }
        if (cor_small(3,2)==0)
        {
            direction[5]=6;
        }
        if (cor_small(3,1)==0)
        {
            direction[6]=7;
        }
        if (cor_small(2,1)==0)
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

inline void random_migration(int i, double deltah,Array<double, 2> &cell_array, Array<long, 3> &Visual_range, Array<int,2> &cor_big, Array<int, 2> &area_square, Array<int, 2> &sub_area_square, Array<int, 2> &cor_small, Array<int, 2> &area_square_s, Array<int, 2>  &sub_area_square_s,double &migration_judgement, long rng_time_step, long rng_event_base)
{
    (void)deltah;
    (void)area_square;
    (void)sub_area_square;
    (void)area_square_s;
    (void)sub_area_square_s;

    long cell_rng_id = (long)cell_array(i,cell_col::kId);
    if (cell_rng_id == 0)
    {
        cell_rng_id = i;
    }
    int cell_stage=(int)cell_array(i,cell_col::kStage);
    long rng_event = rng_event_base + ((long)cell_stage * 100);
    int x1=(int)cell_array(i,cell_col::kX1);
    int y1=(int)cell_array(i,cell_col::kY1);
    int order = select_random_migration_direction(x1, y1, cell_stage, Visual_range, cor_big, cor_small, cell_rng_id, rng_time_step, rng_event);
    if (order!=0)
    {
        move_cell_array(i, cell_array, Visual_range, order);
        cell_array(i,cell_col::kMigrationElapsed)=0;
        cell_array(i,cell_col::kMigrationFollowFlag)=1;
    }
    migration_judgement=migration_judgement+0.0001;
}

inline void random_migration(int i, double deltah, CellStore &cells, Array<long, 3> &Visual_range, Array<int,2> &cor_big, Array<int, 2> &area_square, Array<int, 2> &sub_area_square, Array<int, 2> &cor_small, Array<int, 2> &area_square_s, Array<int, 2>  &sub_area_square_s,double &migration_judgement, long rng_time_step, long rng_event_base)
{
    (void)deltah;
    (void)area_square;
    (void)sub_area_square;
    (void)area_square_s;
    (void)sub_area_square_s;

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
    int order = select_random_migration_direction(x1, y1, cell_stage, Visual_range, cor_big, cor_small, cell_rng_id, rng_time_step, rng_event);
    if (order!=0)
    {
        move_cell_store(i, cells, Visual_range, order);
        cells.migration_elapsed()[row]=0;
        cells.migration_follow_flag()[row]=1;
    }
    migration_judgement=migration_judgement+0.0001;
}

#endif /* random_migration_hpp */
