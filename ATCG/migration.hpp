//
//  migration.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//
//
//
//
//90% probability to keep the migration direction as the last directioin.
//5% probability of each side-direction.
//50% probability of each side-side-direction.

//example:
//
// direction 1 is the last migration direction.
//
// 1 2 3
// 8 * 4
// 7 6 5
//
// if 1,2,8 are all aviable to migrate. 90% to 1, 5% to 2 and 5% to 8.
// if 2 and 8 are all aviable to migrate. 50% to 2 and 50% to 8.
// if 1 and 2 OR 1 and 8 are all aviable to migrate. 90% to 1, and 10% to 2 OR 8.
// if if 1,2,8 are all NOT aviable to migrate, migration stoped.



#ifndef migration_hpp
#define migration_hpp

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
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "visual_range.hpp"
#include "cell_columns.hpp"
#include "deltah_calculation.hpp"
#include "cell_motion.hpp"
#include "stateless_rng.hpp"
#include <chrono>

using std::chrono::high_resolution_clock;
using namespace blitz;

inline bool migration_visual_site_empty(const VisualRange &Visual_range, int x, int y)
{
    return Visual_range.occupied(x, y) == 0;
}

inline void fill_big_migration_directions(int x1, int y1, const VisualRange &Visual_range, int direction[8])
{
    if (migration_visual_site_empty(Visual_range, x1, y1 - 1) &&
        migration_visual_site_empty(Visual_range, x1 - 1, y1 - 1) &&
        migration_visual_site_empty(Visual_range, x1 - 1, y1))
    {
        direction[0] = 1;
    }
    if (migration_visual_site_empty(Visual_range, x1 - 1, y1) &&
        migration_visual_site_empty(Visual_range, x1 - 1, y1 + 1))
    {
        direction[1] = 2;
    }
    if (migration_visual_site_empty(Visual_range, x1 - 1, y1 + 1) &&
        migration_visual_site_empty(Visual_range, x1 - 1, y1 + 2) &&
        migration_visual_site_empty(Visual_range, x1, y1 + 2))
    {
        direction[2] = 3;
    }
    if (migration_visual_site_empty(Visual_range, x1, y1 + 2) &&
        migration_visual_site_empty(Visual_range, x1 + 1, y1 + 2))
    {
        direction[3] = 4;
    }
    if (migration_visual_site_empty(Visual_range, x1 + 1, y1 + 2) &&
        migration_visual_site_empty(Visual_range, x1 + 2, y1 + 2) &&
        migration_visual_site_empty(Visual_range, x1 + 2, y1 + 1))
    {
        direction[4] = 5;
    }
    if (migration_visual_site_empty(Visual_range, x1 + 2, y1 + 1) &&
        migration_visual_site_empty(Visual_range, x1 + 2, y1))
    {
        direction[5] = 6;
    }
    if (migration_visual_site_empty(Visual_range, x1 + 2, y1) &&
        migration_visual_site_empty(Visual_range, x1 + 2, y1 - 1) &&
        migration_visual_site_empty(Visual_range, x1 + 1, y1 - 1))
    {
        direction[6] = 7;
    }
    if (migration_visual_site_empty(Visual_range, x1 + 1, y1 - 1) &&
        migration_visual_site_empty(Visual_range, x1, y1 - 1))
    {
        direction[7] = 8;
    }
}

inline void fill_small_migration_directions(int x1, int y1, const VisualRange &Visual_range, int direction[8])
{
    if (migration_visual_site_empty(Visual_range, x1 - 1, y1 - 1))
    {
        direction[0] = 1;
    }
    if (migration_visual_site_empty(Visual_range, x1 - 1, y1))
    {
        direction[1] = 2;
    }
    if (migration_visual_site_empty(Visual_range, x1 - 1, y1 + 1))
    {
        direction[2] = 3;
    }
    if (migration_visual_site_empty(Visual_range, x1, y1 + 1))
    {
        direction[3] = 4;
    }
    if (migration_visual_site_empty(Visual_range, x1 + 1, y1 + 1))
    {
        direction[4] = 5;
    }
    if (migration_visual_site_empty(Visual_range, x1 + 1, y1))
    {
        direction[5] = 6;
    }
    if (migration_visual_site_empty(Visual_range, x1 + 1, y1 - 1))
    {
        direction[6] = 7;
    }
    if (migration_visual_site_empty(Visual_range, x1, y1 - 1))
    {
        direction[7] = 8;
    }
}

inline void add_unique_migration_label(int labels[100], int &count, int label)
{
    if (label == 0)
    {
        return;
    }
    for (int idx = 0; idx < count; ++idx)
    {
        if (labels[idx] == label)
        {
            return;
        }
    }
    labels[count++] = label;
}

inline void add_migration_density_site(const VisualRange &Visual_range, int x1, int y1, int local_x, int local_y, int labels[100], int &count)
{
    add_unique_migration_label(labels, count, (int)Visual_range.cell_label(x1 - 5 + local_x, y1 - 5 + local_y));
}

inline double big_migration_density(int x1, int y1, const VisualRange &Visual_range, int direction_index)
{
    int labels[100] = {0};
    int count = 0;
    int denominator = 25;
    switch (direction_index)
    {
        case 0:
            for (int xss = 1; xss <= 5; ++xss)
            {
                for (int yss = 1; yss <= 5; ++yss)
                {
                    add_migration_density_site(Visual_range, x1, y1, xss, yss, labels, count);
                }
            }
            break;
        case 2:
            for (int xss = 1; xss <= 5; ++xss)
            {
                for (int yss = 6; yss <= 10; ++yss)
                {
                    add_migration_density_site(Visual_range, x1, y1, xss, yss, labels, count);
                }
            }
            break;
        case 4:
            for (int xss = 6; xss <= 10; ++xss)
            {
                for (int yss = 6; yss <= 10; ++yss)
                {
                    add_migration_density_site(Visual_range, x1, y1, xss, yss, labels, count);
                }
            }
            break;
        case 6:
            for (int xss = 6; xss <= 10; ++xss)
            {
                for (int yss = 1; yss <= 5; ++yss)
                {
                    add_migration_density_site(Visual_range, x1, y1, xss, yss, labels, count);
                }
            }
            break;
        case 1:
            denominator = 30;
            for (int xss = 1, deltay = 0; xss <= 5; ++xss, ++deltay)
            {
                for (int yss = 1 + deltay; yss <= 10 - deltay; ++yss)
                {
                    add_migration_density_site(Visual_range, x1, y1, xss, yss, labels, count);
                }
            }
            break;
        case 3:
            denominator = 30;
            for (int yss = 10, deltay = 0; yss >= 6; --yss, ++deltay)
            {
                for (int xss = 10 - deltay; xss >= 1 + deltay; --xss)
                {
                    add_migration_density_site(Visual_range, x1, y1, xss, yss, labels, count);
                }
            }
            break;
        case 5:
            denominator = 30;
            for (int xss = 10, deltay = 0; xss >= 6; --xss, ++deltay)
            {
                for (int yss = 1 + deltay; yss <= 10 - deltay; ++yss)
                {
                    add_migration_density_site(Visual_range, x1, y1, xss, yss, labels, count);
                }
            }
            break;
        case 7:
            denominator = 30;
            for (int yss = 1, deltay = 0; yss <= 5; ++yss, ++deltay)
            {
                for (int xss = 1 + deltay; xss <= 10 - deltay; ++xss)
                {
                    add_migration_density_site(Visual_range, x1, y1, xss, yss, labels, count);
                }
            }
            break;
    }
    return (double)count / (double)denominator;
}

inline double small_migration_density(int x1, int y1, const VisualRange &Visual_range, int direction_index)
{
    int labels[100] = {0};
    int count = 0;
    switch (direction_index)
    {
        case 0:
            for (int xss = 1; xss <= 5; ++xss)
            {
                for (int yss = 1; yss <= 5; ++yss)
                {
                    add_migration_density_site(Visual_range, x1, y1, xss, yss, labels, count);
                }
            }
            break;
        case 2:
            for (int xss = 1; xss <= 5; ++xss)
            {
                for (int yss = 5; yss <= 9; ++yss)
                {
                    add_migration_density_site(Visual_range, x1, y1, xss, yss, labels, count);
                }
            }
            break;
        case 4:
            for (int xss = 5; xss <= 9; ++xss)
            {
                for (int yss = 5; yss <= 9; ++yss)
                {
                    add_migration_density_site(Visual_range, x1, y1, xss, yss, labels, count);
                }
            }
            break;
        case 6:
            for (int xss = 5; xss <= 9; ++xss)
            {
                for (int yss = 1; yss <= 5; ++yss)
                {
                    add_migration_density_site(Visual_range, x1, y1, xss, yss, labels, count);
                }
            }
            break;
        case 1:
            for (int xss = 1, deltay = 0; xss <= 5; ++xss, ++deltay)
            {
                for (int yss = 1 + deltay; yss <= 9 - deltay; ++yss)
                {
                    add_migration_density_site(Visual_range, x1, y1, xss, yss, labels, count);
                }
            }
            break;
        case 3:
            for (int yss = 9, deltay = 0; yss >= 5; --yss, ++deltay)
            {
                for (int xss = 1 + deltay; xss <= 9 - deltay; ++xss)
                {
                    add_migration_density_site(Visual_range, x1, y1, xss, yss, labels, count);
                }
            }
            break;
        case 5:
            for (int xss = 9, deltay = 0; xss >= 5; --xss, ++deltay)
            {
                for (int yss = 1 + deltay; yss <= 9 - deltay; ++yss)
                {
                    add_migration_density_site(Visual_range, x1, y1, xss, yss, labels, count);
                }
            }
            break;
        case 7:
            for (int yss = 1, deltay = 0; yss <= 5; ++yss, ++deltay)
            {
                for (int xss = 1 + deltay; xss <= 9 - deltay; ++xss)
                {
                    add_migration_density_site(Visual_range, x1, y1, xss, yss, labels, count);
                }
            }
            break;
    }
    return (double)count / 25.0;
}

template <typename CellArray>
inline void apply_big_migration_move(int i, CellArray &cell_array, VisualRange &Visual_range, int x1, int y1, int order)
{
    long cell_label_1=Visual_range.cell_label(x1,y1);
    long cellstage=Visual_range.stage(x1,y1);
    Visual_range.clear_square(x1,y1);
    int dx = migration_direction_dx(order);
    int dy = migration_direction_dy(order);
    for (int cor_cell=1; cor_cell<=4; cor_cell++)
    {
        int cor_cell_y=cor_cell+4;
        cell_array(i,cor_cell)=cell_array(i,cor_cell)+dx;
        cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)+dy;
    }
    Visual_range.write_square((int)cell_array(i,1), (int)cell_array(i,5), (long)cell_array(i,15), cellstage, cell_label_1);
    cell_array(i,23)=order;
    cell_array(i,20)=0;
}

template <typename CellArray>
inline void apply_small_migration_move(int i, CellArray &cell_array, VisualRange &Visual_range, int x1, int y1, int order)
{
    long cell_label_1=Visual_range.cell_label(x1,y1);
    long cellstage=Visual_range.stage(x1,y1);
    Visual_range.clear_site(x1,y1);
    int dx = migration_direction_dx(order);
    int dy = migration_direction_dy(order);
    cell_array(i,1)=cell_array(i,1)+dx;
    cell_array(i,5)=cell_array(i,5)+dy;
    Visual_range.write_site((int)cell_array(i,1), (int)cell_array(i,5), (long)cell_array(i,15), cellstage, cell_label_1);
    cell_array(i,23)=order;
    cell_array(i,20)=0;
}

template <typename CellArray>
inline void migration(int i, double deltah, CellArray &cell_array, VisualRange &Visual_range, double &migration_judgement, long rng_time_step, long rng_event_base)
{
    (void)deltah;
    int x1=cell_array(i,cell_col::kX1);
    int y1=cell_array(i,cell_col::kY1);
    long cell_rng_id = (long)cell_array(i,cell_col::kId);
    if (cell_rng_id == 0)
    {
        cell_rng_id = i;
    }
    long rng_event = rng_event_base + ((long)cell_array(i,cell_col::kType) * 10000) + ((long)cell_array(i,cell_col::kStage) * 1000);
    int cell_type=(int)cell_array(i,cell_col::kType);
    switch (cell_type)
    {
        case 1://r cells
        {
            int cell_shape=cell_array(i,cell_col::kStage);
            switch (cell_shape)
            {
                case 0: //big
                {
                    int direction[8]={0};
                    fill_big_migration_directions(x1, y1, Visual_range, direction);
                    int mloci=0;
                    for (int mlo=0; mlo<8; mlo++)
                    {
                        if (direction[mlo]!=0)
                        {
                            mloci++;
                        }
                    }
                    if (mloci>0)
                    {
                        int *direction1=new int[mloci];
                        int new_loci=0;
                        for (int loci=0; loci<8; loci++)
                        {
                            if (direction[loci]!=0)
                            {
                                direction1[new_loci]=direction[loci];
                                new_loci++;
                            }
                        }
                        double density[8]={0};
                        int order=0;
                        double mean_density=0.6;
                        ////////////////////////////////////////////////initial migration direction dudgement////////////////////////////////////////////
                        int migration_direction=cell_array(i,cell_col::kMigrationDirection);
                        switch (migration_direction)
                        {
                            case 0:
                            {
                                //////////////////////////////////////////////8 direction density calculation////////////////////////////////////
                                for (int loci_for_mig=0; loci_for_mig<8;loci_for_mig++)
                                {
                                    density[loci_for_mig] = big_migration_density(x1, y1, Visual_range, loci_for_mig);
                                }
                                
                                int new_direction_number=0;
                                for (int dl=0;dl<8;dl++)
                                {
                                    if (density[dl]<=mean_density)
                                    {
                                        new_direction_number=new_direction_number+1;
                                    }
                                }
                                if (new_direction_number>0)
                                {
                                    int *new_direction_for_migration=new int[new_direction_number];
                                    int new_direction_number_for_migration=0;
                                    for (int dl=0;dl<8;dl++)
                                    {
                                        if (density[dl]<=mean_density)
                                        {
                                            new_direction_for_migration[new_direction_number_for_migration]=direction[dl];
                                            new_direction_number_for_migration=new_direction_number_for_migration+1;
                                        }
                                    }
                                    int new_loci_for_migratio_number=0;
                                    for (int density_nonzero=0; density_nonzero<new_direction_number_for_migration; density_nonzero++)
                                    {
                                        if (new_direction_for_migration[density_nonzero]>0)
                                        {
                                            new_loci_for_migratio_number=new_loci_for_migratio_number+1;
                                        }
                                    }
                                    int *new_loci_for_migration=new int[new_loci_for_migratio_number];
                                    int locinumber=0;
                                    for (int loci_mig=0;loci_mig<new_direction_number_for_migration;loci_mig++)
                                    {
                                        if (new_direction_for_migration[loci_mig]>0)
                                        {
                                            new_loci_for_migration[locinumber]=new_direction_for_migration[loci_mig];
                                            locinumber=locinumber+1;
                                        }
                                    }
                                    if (locinumber>0)
                                    {
                                        stateless_shuffle(new_loci_for_migration, new_loci_for_migration + locinumber, cell_rng_id, rng_time_step, rng_event++);
                                        order=new_loci_for_migration[0];
                                    }
                                    delete[] new_direction_for_migration;
                                    new_direction_for_migration = NULL;
                                    delete[] new_loci_for_migration;
                                    new_loci_for_migration = NULL;
                                    cell_array(i,24)=1;
                                }
                                break;
                            }
                            case 1:///////////////////////////////////////////////////////following migration direction dudgement////////////////////////////////////////////
                            {
                                int pro_loci_left=0;
                                int pro_loci_right=0;
                                int pro_loci_mid=0;
//                                int pro_loci_right_riht=0;
//                                int pro_loci_left_left=0;
                                for (int x=0; x<new_loci; x++)
                                {
                                    if (direction1[x]==2)
                                    {
                                        pro_loci_right++;
                                    }
                                    else if (direction1[x]==8)
                                    {
                                        pro_loci_left++;
                                    }
                                    else if (direction1[x]==1)
                                    {
                                        pro_loci_mid++;
                                    }
//                                    else if (direction1[x]==7)
//                                    {
//                                        pro_loci_left_left++;
//                                    }
//                                    else if (direction1[x]==3)
//                                    {
//                                        pro_loci_right_riht++;
//                                    }
                                }
                                if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                                {
                                    int order_loci[20]={8,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
                                    int n=20;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    int order_loci[2]={2,8};
                                    int n=2;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                                {
                                    int order_loci[10]={1,2,1,1,1,1,1,1,1,1};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                                {
                                    int order_loci[10]={1,8,1,1,1,1,1,1,1,1};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                                {
                                    order=8;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                                {
                                    order=1;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    order=2;
                                }
//                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
//                                {
//                                    if(pro_loci_right_riht==1 && pro_loci_left_left==1)
//                                    {
//                                        int order_loci[2]={3,7};
//                                        int n=2;
//                                        stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
//                                        order=order_loci[0];
//                                    }
//                                    else if(pro_loci_right_riht==1 && pro_loci_left_left==0)
//                                    {
//                                        order=3;
//                                    }
//                                    else if(pro_loci_right_riht==0 && pro_loci_left_left==1)
//                                    {
//                                        order=7;
//                                    }
//                                    //                                    cell_array(i,20)=0;
//                                    //                                    cell_array(i,23)=0;
//                                }
                                break;
                            }
                            case 8:
                            {
                                int pro_loci_left=0;
                                int pro_loci_right=0;
                                int pro_loci_mid=0;
//                                int pro_loci_right_riht=0;
//                                int pro_loci_left_left=0;
                                for (int x=0; x<new_loci; x++)
                                {
                                    if (direction1[x]==7)
                                    {
                                        pro_loci_left++;
                                    }
                                    else if (direction1[x]==1)
                                    {
                                        pro_loci_right++;
                                    }
                                    else if (direction1[x]==8)
                                    {
                                        pro_loci_mid++;
                                    }
//                                    else if (direction1[x]==6)
//                                    {
//                                        pro_loci_left_left++;
//                                    }
//                                    else if (direction1[x]==2)
//                                    {
//                                        pro_loci_right_riht++;
//                                    }
                                }
                                if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                                {
                                    int order_loci[20]={7,8,1,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8};
                                    int n=20;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    int order_loci[2]={7,1};
                                    int n=2;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                                {
                                    int order_loci[10]={7,8,8,8,8,8,8,8,8,8};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                                {
                                    int order_loci[10]={1,8,8,8,8,8,8,8,8,8};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                                {
                                    order=7;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                                {
                                    order=8;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    order=1;
                                }
//                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
//                                {
//                                    if(pro_loci_right_riht==1 && pro_loci_left_left==1)
//                                    {
//                                        int order_loci[2]={2,6};
//                                        int n=2;
//                                        stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
//                                        order=order_loci[0];
//                                    }
//                                    else if(pro_loci_right_riht==1 && pro_loci_left_left==0)
//                                    {
//                                        order=2;
//                                    }
//                                    else if(pro_loci_right_riht==0 && pro_loci_left_left==1)
//                                    {
//                                        order=6;
//                                    }
//                                    //                                    cell_array(i,20)=0;
//                                    //                                    cell_array(i,23)=0;
//                                }
                                break;
                            }
                            case 2:
                            {
                                int pro_loci_left=0;
                                int pro_loci_right=0;
                                int pro_loci_mid=0;
//                                int pro_loci_right_riht=0;
//                                int pro_loci_left_left=0;
                                for (int x=0; x<new_loci; x++)
                                {
                                    if (direction1[x]==1)
                                    {
                                        pro_loci_left++;
                                    }
                                    else if (direction1[x]==3)
                                    {
                                        pro_loci_right++;
                                    }
                                    else if (direction1[x]==2)
                                    {
                                        pro_loci_mid++;
                                    }
//                                    else if (direction1[x]==8)
//                                    {
//                                        pro_loci_left_left++;
//                                    }
//                                    else if (direction1[x]==4)
//                                    {
//                                        pro_loci_right_riht++;
//                                    }
                                }
                                if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                                {
                                    int order_loci[20]={1,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
                                    int n=20;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    int order_loci[2]={1,3};
                                    int n=2;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                                {
                                    int order_loci[10]={1,2,2,2,2,2,2,2,2,2};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                                {
                                    int order_loci[10]={3,2,2,2,2,2,2,2,2,2};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                                {
                                    order=1;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                                {
                                    order=2;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    order=3;
                                }
//                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
//                                {
//                                    if(pro_loci_right_riht==1 && pro_loci_left_left==1)
//                                    {
//                                        int order_loci[2]={8,4};
//                                        int n=2;
//                                        stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
//                                        order=order_loci[0];
//                                    }
//                                    else if(pro_loci_right_riht==1 && pro_loci_left_left==0)
//                                    {
//                                        order=4;
//                                    }
//                                    else if(pro_loci_right_riht==0 && pro_loci_left_left==1)
//                                    {
//                                        order=8;
//                                    }
//                                    //                                    cell_array(i,20)=0;
//                                    //                                    cell_array(i,23)=0;
//                                }
                                break;
                            }
                            case 3:
                            {
                                int pro_loci_left=0;
                                int pro_loci_right=0;
                                int pro_loci_mid=0;
//                                int pro_loci_right_riht=0;
//                                int pro_loci_left_left=0;
                                for (int x=0; x<new_loci; x++)
                                {
                                    if (direction1[x]==2)
                                    {
                                        pro_loci_left++;
                                    }
                                    else if (direction1[x]==4)
                                    {
                                        pro_loci_right++;
                                    }
                                    else if (direction1[x]==3)
                                    {
                                        pro_loci_mid++;
                                    }
//                                    else if (direction1[x]==1)
//                                    {
//                                        pro_loci_right_riht++;
//                                    }
//                                    else if (direction1[x]==5)
//                                    {
//                                        pro_loci_left_left++;
//                                    }
                                }
                                if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                                {
                                    int order_loci[20]={2,3,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
                                    int n=20;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    int order_loci[2]={2,4};
                                    int n=2;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                                {
                                    int order_loci[10]={2,3,3,3,3,3,3,3,3,3};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                                {
                                    int order_loci[10]={4,3,3,3,3,3,3,3,3,3};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                                {
                                    order=2;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                                {
                                    order=3;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    order=4;
                                }
//                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
//                                {
//
//                                    if(pro_loci_right_riht==1 && pro_loci_left_left==1)
//                                    {
//                                        int order_loci[2]={1,5};
//                                        int n=2;
//                                        stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
//                                        order=order_loci[0];
//                                    }
//                                    else if(pro_loci_right_riht==1 && pro_loci_left_left==0)
//                                    {
//                                        order=1;
//                                    }
//                                    else if(pro_loci_right_riht==0 && pro_loci_left_left==1)
//                                    {
//                                        order=5;
//                                    }
//                                    //                                    cell_array(i,20)=0;
//                                    //                                    cell_array(i,23)=0;
//                                }
                                break;
                            }
                            case 4:
                            {
                                int pro_loci_left=0;
                                int pro_loci_right=0;
                                int pro_loci_mid=0;
//                                int pro_loci_right_riht=0;
//                                int pro_loci_left_left=0;
                                for (int x=0; x<new_loci; x++)
                                {
                                    if (direction1[x]==3)
                                    {
                                        pro_loci_left++;
                                    }
                                    else if (direction1[x]==5)
                                    {
                                        pro_loci_right++;
                                    }
                                    else if (direction1[x]==4)
                                    {
                                        pro_loci_mid++;
                                    }
//                                    else if (direction1[x]==2)
//                                    {
//                                        pro_loci_left_left++;
//                                    }
//                                    else if (direction1[x]==6)
//                                    {
//                                        pro_loci_right_riht++;
//                                    }
                                }
                                if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                                {
                                    int order_loci[20]={3,4,5,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};
                                    int n=20;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    int order_loci[2]={3,5};
                                    int n=2;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                                {
                                    int order_loci[10]={3,4,4,4,4,4,4,4,4,4};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                                {
                                    int order_loci[10]={5,4,4,4,4,4,4,4,4,4};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                                {
                                    order=3;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                                {
                                    order=4;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    order=5;
                                }
//                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
//                                {
//                                    if(pro_loci_right_riht==1 && pro_loci_left_left==1)
//                                    {
//                                        int order_loci[2]={2,6};
//                                        int n=2;
//                                        stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
//                                        order=order_loci[0];
//                                    }
//                                    else if(pro_loci_right_riht==1 && pro_loci_left_left==0)
//                                    {
//                                        order=6;
//                                    }
//                                    else if(pro_loci_right_riht==0 && pro_loci_left_left==1)
//                                    {
//                                        order=2;
//                                    }
//                                    //                                    cell_array(i,20)=0;
//                                    //                                    cell_array(i,23)=0;
//                                }
                                break;
                            }
                            case 5:
                            {
                                int pro_loci_left=0;
                                int pro_loci_right=0;
                                int pro_loci_mid=0;
//                                int pro_loci_right_riht=0;
//                                int pro_loci_left_left=0;
                                for (int x=0; x<new_loci; x++)
                                {
                                    if (direction1[x]==4)
                                    {
                                        pro_loci_left++;
                                    }
                                    else if (direction1[x]==6)
                                    {
                                        pro_loci_right++;
                                    }
                                    else if (direction1[x]==5)
                                    {
                                        pro_loci_mid++;
                                    }
//                                    else if (direction1[x]==3)
//                                    {
//                                        pro_loci_left_left++;
//                                    }
//                                    else if (direction1[x]==7)
//                                    {
//                                        pro_loci_right_riht++;
//                                    }
                                }
                                if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                                {
                                    int order_loci[20]={4,5,6,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5};
                                    int n=20;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    int order_loci[2]={4,6};
                                    int n=2;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                                {
                                    int order_loci[10]={4,5,5,5,5,5,5,5,5,5};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                                {
                                    int order_loci[10]={6,5,5,5,5,5,5,5,5,5};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                                {
                                    order=4;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                                {
                                    order=5;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    order=6;
                                }
//                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
//                                {
//
//                                    if(pro_loci_right_riht==1 && pro_loci_left_left==1)
//                                    {
//                                        int order_loci[2]={3,7};
//                                        int n=2;
//                                        stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
//                                        order=order_loci[0];
//                                    }
//                                    else if(pro_loci_right_riht==1 && pro_loci_left_left==0)
//                                    {
//                                        order=7;
//                                    }
//                                    else if(pro_loci_right_riht==0 && pro_loci_left_left==1)
//                                    {
//                                        order=3;
//                                    }
//                                    //                                    cell_array(i,20)=0;
//                                    //                                    cell_array(i,23)=0;
//                                }
                                break;
                            }
                            case 6:
                            {
                                int pro_loci_left=0;
                                int pro_loci_right=0;
                                int pro_loci_mid=0;
//                                int pro_loci_right_riht=0;
//                                int pro_loci_left_left=0;
                                for (int x=0; x<new_loci; x++)
                                {
                                    if (direction1[x]==5)
                                    {
                                        pro_loci_left++;
                                    }
                                    else if (direction1[x]==7)
                                    {
                                        pro_loci_right++;
                                    }
                                    else if (direction1[x]==6)
                                    {
                                        pro_loci_mid++;
                                    }
//                                    else if (direction1[x]==8)
//                                    {
//                                        pro_loci_right_riht++;
//                                    }
//                                    else if (direction1[x]==4)
//                                    {
//                                        pro_loci_left_left++;
//                                    }
                                }
                                if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                                {
                                    int order_loci[20]={5,6,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6};
                                    int n=20;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    int order_loci[2]={5,7};
                                    int n=2;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                                {
                                    int order_loci[10]={5,6,6,6,6,6,6,6,6,6};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                                {
                                    int order_loci[10]={7,6,6,6,6,6,6,6,6,6};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                                {
                                    order=5;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                                {
                                    order=6;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    order=7;
                                }
//                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
//                                {
//                                    if(pro_loci_right_riht==1 && pro_loci_left_left==1)
//                                    {
//                                        int order_loci[2]={4,8};
//                                        int n=2;
//                                        stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
//                                        order=order_loci[0];
//                                    }
//                                    else if(pro_loci_right_riht==1 && pro_loci_left_left==0)
//                                    {
//                                        order=8;
//                                    }
//                                    else if(pro_loci_right_riht==0 && pro_loci_left_left==1)
//                                    {
//                                        order=4;
//                                    }
//                                    //                                    cell_array(i,20)=0;
//                                    //                                    cell_array(i,23)=0;
//                                }
                                break;
                            }
                            case 7:
                            {
                                int pro_loci_left=0;
                                int pro_loci_right=0;
                                int pro_loci_mid=0;
//                                int pro_loci_right_riht=0;
//                                int pro_loci_left_left=0;
                                for (int x=0; x<new_loci; x++)
                                {
                                    if (direction1[x]==6)
                                    {
                                        pro_loci_left++;
                                    }
                                    else if (direction1[x]==8)
                                    {
                                        pro_loci_right++;
                                    }
                                    else if (direction1[x]==7)
                                    {
                                        pro_loci_mid++;
                                    }
//                                    else if (direction1[x]==1)
//                                    {
//                                        pro_loci_right_riht++;
//                                    }
//                                    else if (direction1[x]==5)
//                                    {
//                                        pro_loci_left_left++;
//                                    }
                                }
                                if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                                {
                                    int order_loci[20]={6,7,8,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7};
                                    int n=20;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    int order_loci[2]={6,8};
                                    int n=2;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                                {
                                    int order_loci[10]={6,7,7,7,7,7,7,7,7,7};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                                {
                                    int order_loci[10]={8,7,7,7,7,7,7,7,7,7};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                                {
                                    order=6;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                                {
                                    order=7;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    order=8;
                                }
//                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
//                                {
//                                    if(pro_loci_right_riht==1 && pro_loci_left_left==1)
//                                    {
//                                        int order_loci[2]={1,5};
//                                        int n=2;
//                                        stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
//                                        order=order_loci[0];
//                                    }
//                                    else if(pro_loci_right_riht==1 && pro_loci_left_left==0)
//                                    {
//                                        order=1;
//                                    }
//                                    else if(pro_loci_right_riht==0 && pro_loci_left_left==1)
//                                    {
//                                        order=5;
//                                    }
//                                    //                                    cell_array(i,20)=0;
//                                    //                                    cell_array(i,23)=0;
//                                }
                                break;
                            }
                        }
                        ////////////////////////////////////////////// migration ////////////////////////////////////
                                                if (order>=1 && order<=8)
                        {
                            apply_big_migration_move(i, cell_array, Visual_range, x1, y1, order);
                        }
                        else
                        {
                            cell_array(i,23)=0;
                            cell_array(i,20)=0;
                        }
                        delete[] direction1;
                        direction1 = NULL;
                        //                        cell_array(i,20)=0;
                    }
                    break;
                }
                default ://small
                {
                    int direction[8]={0};
                    fill_small_migration_directions(x1, y1, Visual_range, direction);
                    int mloci=0;
                    for (int mlo=0; mlo<8; mlo++)
                    {
                        if (direction[mlo]!=0)
                        {
                            mloci++;
                        }
                    }
                    if (mloci>0)
                    {
                        int *direction1=new int[mloci];
                        int new_loci=0;
                        for (int loci=0; loci<8; loci++)
                        {
                            if (direction[loci]!=0)
                            {
                                direction1[new_loci]=direction[loci];
                                new_loci++;
                            }
                        }
                        double density[8]={0};
                        int order=0;
                        double mean_density=0.6;
                        ////////////////////////////////////////////////////initial migration direction dudgement/////////////////////////////////////////////////////////////
                        int migration_direction=(int)cell_array(i,23);
                        switch (migration_direction)
                        {
                            case 0:
                            {
                                /////////////////////////////////////////////8 directions density dudgement/////////////////////////////////////
                                for (int loci_for_mig=0; loci_for_mig<8;loci_for_mig++)
                                {
                                    density[loci_for_mig] = small_migration_density(x1, y1, Visual_range, loci_for_mig);
                                }
                                int new_direction_number=0;
                                for (int dl=0;dl<8;dl++)
                                {
                                    if (density[dl]<=mean_density)
                                    {
                                        new_direction_number=new_direction_number+1;
                                    }
                                }
                                if (new_direction_number>0)
                                {
                                    int *new_direction_for_migration=new int[new_direction_number];
                                    int new_direction_number_for_migration=0;
                                    for (int dl=0;dl<8;dl++)
                                    {
                                        if (density[dl]<=mean_density)
                                        {
                                            new_direction_for_migration[new_direction_number_for_migration]=direction[dl];
                                            new_direction_number_for_migration=new_direction_number_for_migration+1;
                                        }
                                    }
                                    int new_loci_for_migratio_number=0;
                                    for (int density_nonzero=0; density_nonzero<new_direction_number_for_migration; density_nonzero++)
                                    {
                                        if (new_direction_for_migration[density_nonzero]>0)
                                        {
                                            new_loci_for_migratio_number=new_loci_for_migratio_number+1;
                                        }
                                    }
                                    int *new_loci_for_migration=new int[new_loci_for_migratio_number];
                                    int locinumber=0;
                                    for (int loci_mig=0;loci_mig<new_direction_number_for_migration;loci_mig++)
                                    {
                                        if (new_direction_for_migration[loci_mig]>0)
                                        {
                                            new_loci_for_migration[locinumber]=new_direction_for_migration[loci_mig];
                                            locinumber=locinumber+1;
                                        }
                                    }
                                    if (locinumber>0)
                                    {
                                        stateless_shuffle(new_loci_for_migration, new_loci_for_migration + locinumber, cell_rng_id, rng_time_step, rng_event++);
                                        
                                        order=new_loci_for_migration[0];
                                    }
                                    delete[] new_direction_for_migration;
                                    new_direction_for_migration = NULL;
                                    delete[] new_loci_for_migration;
                                    new_loci_for_migration = NULL;
                                    cell_array(i,24)=1;
                                }
                                break;
                            }
                            case 1://////////////////////////////////////////////////////following migration direction dudgement////////////////////////////////////////////////////////////
                            {
                                int pro_loci_left=0;
                                int pro_loci_right=0;
                                int pro_loci_mid=0;
//                                int pro_loci_right_riht=0;
//                                int pro_loci_left_left=0;
                                for (int x=0; x<new_loci; x++)
                                {
                                    if (direction1[x]==8)
                                    {
                                        pro_loci_left++;
                                    }
                                    else if (direction1[x]==2)
                                    {
                                        pro_loci_right++;
                                    }
                                    else if (direction1[x]==1)
                                    {
                                        pro_loci_mid++;
                                    }
//                                    else if (direction1[x]==3)
//                                    {
//                                        pro_loci_right_riht++;
//                                    }
//                                    else if (direction1[x]==7)
//                                    {
//                                        pro_loci_left_left++;
//                                    }
                                }
                                if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                                {
                                    int order_loci[20]={8,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
                                    int n=20;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    int order_loci[2]={2,8};
                                    int n=2;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                                {
                                    int order_loci[10]={1,8,1,1,1,1,1,1,1,1};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                                {
                                    int order_loci[10]={1,2,1,1,1,1,1,1,1,1};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                                {
                                    order=8;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                                {
                                    order=1;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    order=2;
                                }
//                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
//                                {
//                                    if(pro_loci_right_riht==1 && pro_loci_left_left==1)
//                                    {
//                                        int order_loci[2]={3,7};
//                                        int n=2;
//                                        stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
//                                        order=order_loci[0];
//                                    }
//                                    else if(pro_loci_right_riht==1 && pro_loci_left_left==0)
//                                    {
//                                        order=3;
//                                    }
//                                    else if(pro_loci_right_riht==0 && pro_loci_left_left==1)
//                                    {
//                                        order=1;
//                                    }
//                                    //                                    cell_array(i,20)=0;
//                                    //                                    cell_array(i,23)=0;
//                                }
                                break;
                            }
                            case 8:
                            {
                                int pro_loci_left=0;
                                int pro_loci_right=0;
                                int pro_loci_mid=0;
//                                int pro_loci_right_riht=0;
//                                int pro_loci_left_left=0;
                                for (int x=0; x<new_loci; x++)
                                {
                                    if (direction1[x]==7)
                                    {
                                        pro_loci_left++;
                                    }
                                    else if (direction1[x]==1)
                                    {
                                        pro_loci_right++;
                                    }
                                    else if (direction1[x]==8)
                                    {
                                        pro_loci_mid++;
                                    }
//                                    else if (direction1[x]==2)
//                                    {
//                                        pro_loci_right_riht++;
//                                    }
//                                    else if (direction1[x]==6)
//                                    {
//                                        pro_loci_left_left++;
//                                    }
                                }
                                if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                                {
                                    int order_loci[20]={7,8,1,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8};
                                    int n=20;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    int order_loci[2]={7,1};
                                    int n=2;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                                {
                                    int order_loci[10]={7,8,8,8,8,8,8,8,8,8};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                                {
                                    int order_loci[10]={1,8,8,8,8,8,8,8,8,8};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                                {
                                    order=7;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                                {
                                    order=8;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    order=1;
                                }
//                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
//                                {
//                                    if(pro_loci_right_riht==1 && pro_loci_left_left==1)
//                                    {
//                                        int order_loci[2]={2,6};
//                                        int n=2;
//                                        stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
//                                        order=order_loci[0];
//                                    }
//                                    else if(pro_loci_right_riht==1 && pro_loci_left_left==0)
//                                    {
//                                        order=2;
//                                    }
//                                    else if(pro_loci_right_riht==0 && pro_loci_left_left==1)
//                                    {
//                                        order=6;
//                                    }
//                                    //                                    cell_array(i,20)=0;
//                                    //                                    cell_array(i,23)=0;
//                                }
                                break;
                            }
                            case 2:
                            {
                                int pro_loci_left=0;
                                int pro_loci_right=0;
                                int pro_loci_mid=0;
//                                int pro_loci_right_riht=0;
//                                int pro_loci_left_left=0;
                                for (int x=0; x<new_loci; x++)
                                {
                                    if (direction1[x]==1)
                                    {
                                        pro_loci_left++;
                                    }
                                    else if (direction1[x]==3)
                                    {
                                        pro_loci_right++;
                                    }
                                    else if (direction1[x]==2)
                                    {
                                        pro_loci_mid++;
                                    }
//                                    else if (direction1[x]==4)
//                                    {
//                                        pro_loci_right_riht++;
//                                    }
//                                    else if (direction1[x]==8)
//                                    {
//                                        pro_loci_left_left++;
//                                    }
                                }
                                if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                                {
                                    int order_loci[20]={1,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
                                    int n=20;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    int order_loci[2]={1,3};
                                    int n=2;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                                {
                                    int order_loci[10]={1,2,2,2,2,2,2,2,2,2};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                                {
                                    int order_loci[10]={3,2,2,2,2,2,2,2,2,2};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                                {
                                    order=1;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                                {
                                    order=2;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    order=3;
                                }
//                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
//                                {
//                                    if(pro_loci_right_riht==1 && pro_loci_left_left==1)
//                                    {
//                                        int order_loci[2]={8,4};
//                                        int n=2;
//                                        stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
//                                        order=order_loci[0];
//                                    }
//                                    else if(pro_loci_right_riht==1 && pro_loci_left_left==0)
//                                    {
//                                        order=4;
//                                    }
//                                    else if(pro_loci_right_riht==0 && pro_loci_left_left==1)
//                                    {
//                                        order=8;
//                                    }
//                                    //                                    cell_array(i,20)=0;
//                                    //                                    cell_array(i,23)=0;
//                                }
                                break;
                            }
                            case 3:
                            {
                                int pro_loci_left=0;
                                int pro_loci_right=0;
                                int pro_loci_mid=0;
                                int pro_loci_right_riht=0;
                                int pro_loci_left_left=0;
                                for (int x=0; x<new_loci; x++)
                                {
                                    if (direction1[x]==2)
                                    {
                                        pro_loci_left++;
                                    }
                                    else if (direction1[x]==4)
                                    {
                                        pro_loci_right++;
                                    }
                                    else if (direction1[x]==3)
                                    {
                                        pro_loci_mid++;
                                    }
                                    else if (direction1[x]==5)
                                    {
                                        pro_loci_right_riht++;
                                    }
                                    else if (direction1[x]==1)
                                    {
                                        pro_loci_left_left++;
                                    }
                                }
                                if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                                {
                                    int order_loci[20]={2,3,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
                                    int n=20;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    int order_loci[2]={2,4};
                                    int n=2;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                                {
                                    int order_loci[10]={2,3,3,3,3,3,3,3,3,3};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                                {
                                    int order_loci[10]={4,3,3,3,3,3,3,3,3,3};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                                {
                                    order=2;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                                {
                                    order=3;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    order=4;
                                }
//                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
//                                {
//
//                                    if(pro_loci_right_riht==1 && pro_loci_left_left==1)
//                                    {
//                                        int order_loci[2]={1,5};
//                                        int n=2;
//                                        stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
//                                        order=order_loci[0];
//                                    }
//                                    else if(pro_loci_right_riht==1 && pro_loci_left_left==0)
//                                    {
//                                        order=5;
//                                    }
//                                    else if(pro_loci_right_riht==0 && pro_loci_left_left==1)
//                                    {
//                                        order=1;
//                                    }
//                                    //                                    cell_array(i,20)=0;
//                                    //                                    cell_array(i,23)=0;
//                                }
                                break;
                            }
                            case 4:
                            {
                                int pro_loci_left=0;
                                int pro_loci_right=0;
                                int pro_loci_mid=0;
//                                int pro_loci_right_riht=0;
//                                int pro_loci_left_left=0;
                                for (int x=0; x<new_loci; x++)
                                {
                                    if (direction1[x]==3)
                                    {
                                        pro_loci_left++;
                                    }
                                    else if (direction1[x]==5)
                                    {
                                        pro_loci_right++;
                                    }
                                    else if (direction1[x]==4)
                                    {
                                        pro_loci_mid++;
                                    }
//                                    else if (direction1[x]==6)
//                                    {
//                                        pro_loci_right_riht++;
//                                    }
//                                    else if (direction1[x]==2)
//                                    {
//                                        pro_loci_left_left++;
//                                    }
                                }
                                if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                                {
                                    int order_loci[20]={3,4,5,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};
                                    int n=20;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    int order_loci[2]={3,5};
                                    int n=2;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                                {
                                    int order_loci[10]={3,4,4,4,4,4,4,4,4,4};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                                {
                                    int order_loci[10]={5,4,4,4,4,4,4,4,4,4};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                                {
                                    order=3;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                                {
                                    order=4;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    order=5;
                                }
//                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
//                                {
//                                    if(pro_loci_right_riht==1 && pro_loci_left_left==1)
//                                    {
//                                        int order_loci[2]={2,6};
//                                        int n=2;
//                                        stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
//                                        order=order_loci[0];
//                                    }
//                                    else if(pro_loci_right_riht==1 && pro_loci_left_left==0)
//                                    {
//                                        order=6;
//                                    }
//                                    else if(pro_loci_right_riht==0 && pro_loci_left_left==1)
//                                    {
//                                        order=2;
//                                    }
//                                    //                                    cell_array(i,20)=0;
//                                    //                                    cell_array(i,23)=0;
//                                }
                                break;
                            }
                            case 5:
                            {
                                int pro_loci_left=0;
                                int pro_loci_right=0;
                                int pro_loci_mid=0;
//                                int pro_loci_right_riht=0;
//                                int pro_loci_left_left=0;
                                for (int x=0; x<new_loci; x++)
                                {
                                    if (direction1[x]==4)
                                    {
                                        pro_loci_left++;
                                    }
                                    else if (direction1[x]==6)
                                    {
                                        pro_loci_right++;
                                    }
                                    else if (direction1[x]==5)
                                    {
                                        pro_loci_mid++;
                                    }
//                                    else if (direction1[x]==7)
//                                    {
//                                        pro_loci_right_riht++;
//                                    }
//                                    else if (direction1[x]==3)
//                                    {
//                                        pro_loci_left_left++;
//                                    }
                                }
                                if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                                {
                                    int order_loci[20]={4,5,6,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5};
                                    int n=20;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    int order_loci[2]={4,6};
                                    int n=2;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                                {
                                    int order_loci[10]={4,5,5,5,5,5,5,5,5,5};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                                {
                                    int order_loci[10]={6,5,5,5,5,5,5,5,5,5};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                                {
                                    order=4;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                                {
                                    order=5;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    order=6;
                                }
//                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
//                                {
//
//                                    if(pro_loci_right_riht==1 && pro_loci_left_left==1)
//                                    {
//                                        int order_loci[2]={3,7};
//                                        int n=2;
//                                        stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
//                                        order=order_loci[0];
//                                    }
//                                    else if(pro_loci_right_riht==1 && pro_loci_left_left==0)
//                                    {
//                                        order=7;
//                                    }
//                                    else if(pro_loci_right_riht==0 && pro_loci_left_left==1)
//                                    {
//                                        order=3;
//                                    }
//                                    //                                    cell_array(i,20)=0;
//                                    //                                    cell_array(i,23)=0;
//                                }
                                break;
                            }
                            case 6:
                            {
                                int pro_loci_left=0;
                                int pro_loci_right=0;
                                int pro_loci_mid=0;
//                                int pro_loci_right_riht=0;
//                                int pro_loci_left_left=0;
                                for (int x=0; x<new_loci; x++)
                                {
                                    if (direction1[x]==5)
                                    {
                                        pro_loci_left++;
                                    }
                                    else if (direction1[x]==7)
                                    {
                                        pro_loci_right++;
                                    }
                                    else if (direction1[x]==6)
                                    {
                                        pro_loci_mid++;
                                    }
//                                    else if (direction1[x]==8)
//                                    {
//                                        pro_loci_right_riht++;
//                                    }
//                                    else if (direction1[x]==4)
//                                    {
//                                        pro_loci_left_left++;
//                                    }
                                }
                                if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                                {
                                    int order_loci[20]={5,6,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6};
                                    int n=20;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    int order_loci[2]={5,7};
                                    int n=2;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                                {
                                    int order_loci[10]={5,6,6,6,6,6,6,6,6,6};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                                {
                                    int order_loci[10]={7,6,6,6,6,6,6,6,6,6};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                                {
                                    order=5;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                                {
                                    order=6;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    order=7;
                                }
//                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
//                                {
//                                    if(pro_loci_right_riht==1 && pro_loci_left_left==1)
//                                    {
//                                        int order_loci[2]={4,8};
//                                        int n=2;
//                                        stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
//                                        order=order_loci[0];
//                                    }
//                                    else if(pro_loci_right_riht==1 && pro_loci_left_left==0)
//                                    {
//                                        order=8;
//                                    }
//                                    else if(pro_loci_right_riht==0 && pro_loci_left_left==1)
//                                    {
//                                        order=4;
//                                    }
//                                    //                                    cell_array(i,20)=0;
//                                    //                                    cell_array(i,23)=0;
//                                }
                                break;
                            }
                            case 7:
                            {
                                int pro_loci_left=0;
                                int pro_loci_right=0;
                                int pro_loci_mid=0;
//                                int pro_loci_right_riht=0;
//                                int pro_loci_left_left=0;
                                for (int x=0; x<new_loci; x++)
                                {
                                    if (direction1[x]==6)
                                    {
                                        pro_loci_left++;
                                    }
                                    else if (direction1[x]==8)
                                    {
                                        pro_loci_right++;
                                    }
                                    else if (direction1[x]==7)
                                    {
                                        pro_loci_mid++;
                                    }
//                                    else if (direction1[x]==1)
//                                    {
//                                        pro_loci_right_riht++;
//                                    }
//                                    else if (direction1[x]==5)
//                                    {
//                                        pro_loci_left_left++;
//                                    }
                                }
                                if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                                {
                                    int order_loci[20]={6,7,8,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7};
                                    int n=20;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    int order_loci[2]={6,8};
                                    int n=2;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                                {
                                    int order_loci[10]={6,7,7,7,7,7,7,7,7,7};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                                {
                                    int order_loci[10]={8,7,7,7,7,7,7,7,7,7};
                                    int n=10;
                                    stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
                                    order=order_loci[0];
                                }
                                else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                                {
                                    order=6;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                                {
                                    order=7;
                                }
                                else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                                {
                                    order=8;
                                }
//                                else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
//                                {
//                                    if(pro_loci_right_riht==1 && pro_loci_left_left==1)
//                                    {
//                                        int order_loci[2]={1,5};
//                                        int n=2;
//                                        stateless_shuffle(order_loci, order_loci + n, cell_rng_id, rng_time_step, rng_event++);
//                                        order=order_loci[0];
//                                    }
//                                    else if(pro_loci_right_riht==1 && pro_loci_left_left==0)
//                                    {
//                                        order=1;
//                                    }
//                                    else if(pro_loci_right_riht==0 && pro_loci_left_left==1)
//                                    {
//                                        order=5;
//                                    }
//                                    //                                    cell_array(i,20)=0;
//                                    //                                    cell_array(i,23)=0;
//                                }
                                break;
                            }
                        }
                        /////////////////////////////////////////////////////migration///////////////////////////////////////
                                                if (order>=1 && order<=8)
                        {
                            apply_small_migration_move(i, cell_array, Visual_range, x1, y1, order);
                        }
                        else
                        {
                            cell_array(i,23)=0;
                            cell_array(i,20)=0;
                        }
                        delete[] direction1;
                        direction1 = NULL;
                        //                        cell_array(i,20)=0;
                    }
                    break;
                }
            }
            break;
        }
        case 2://K cells
        {
            int cell_shape=(int)cell_array(i,14);
            switch (cell_shape)
            {
                case 0:
                {
                    int direction[8]={0};
                    fill_big_migration_directions(x1, y1, Visual_range, direction);
                    int mloci=0;
                    for (int mlo=0; mlo<8; mlo++)
                    {
                        if (direction[mlo]!=0)
                        {
                            mloci++;
                        }
                    }
                    ////////////////////////////////////////////////migration direction dudgement////////////////////////////////////////////
                    if (mloci>0)
                    {
                        int *direction1=new int[mloci];
                        int new_loci=0;
                        for (int loci=0; loci<8; loci++)
                        {
                            if (direction[loci]!=0)
                            {
                                direction1[new_loci]=direction[loci];
                                new_loci++;
                            }
                        }
                        int order=0;
                        long length_dir=new_loci;
                        if (length_dir==0)
                        {
                            order=direction1[0];
                        }
                        else
                        {
                            stateless_shuffle(direction1, direction1 + new_loci, cell_rng_id, rng_time_step, rng_event++);
                            
                            order=direction1[0];
                        }
                        ////////////////////////////////////////////////////////////////////migration//////////////////////////////////////////////////////
                                                if (order>=1 && order<=8)
                        {
                            apply_big_migration_move(i, cell_array, Visual_range, x1, y1, order);
                        }
                        else
                        {
                            cell_array(i,23)=0;
                            cell_array(i,20)=0;
                        }
                        delete[] direction1;
                        direction1 = NULL;
                        //                        cell_array(i,20)=0;
                        cell_array(i,24)=1;
                    }
                    break;
                }
                default :
                {
                    int direction[8]={0};
                    fill_small_migration_directions(x1, y1, Visual_range, direction);
                    int mloci=0;
                    for (int mlo=0; mlo<8; mlo++)
                    {
                        if (direction[mlo]!=0)
                        {
                            mloci++;
                        }
                    }
                    if (mloci>0)
                    {
                        int *direction1=new int[mloci];
                        int new_loci=0;
                        for (int loci=0; loci<8; loci++)
                        {
                            if (direction[loci]!=0)
                            {
                                direction1[new_loci]=direction[loci];
                                new_loci++;
                            }
                        }
                        int order=0;
                        long length_dir=new_loci;
                        if (length_dir==0)
                        {
                            order=direction1[0];
                        }
                        else if (length_dir>0)
                        {
                            stateless_shuffle(direction1, direction1 + new_loci, cell_rng_id, rng_time_step, rng_event++);
                            
                            order=direction1[0];
                        }
                                                if (order>=1 && order<=8)
                        {
                            apply_small_migration_move(i, cell_array, Visual_range, x1, y1, order);
                        }
                        else
                        {
                            cell_array(i,23)=0;
                            cell_array(i,20)=0;
                        }
                        delete[] direction1;
                        direction1 = NULL;
                        //                        cell_array(i,20)=0;
                        cell_array(i,24)=1;
                    }
                    break;
                }
            }
            break;
        }
    }
    migration_judgement=migration_judgement+1;
}
#endif /* migration_hpp */
