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



#include "migration.hpp"

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
#include "visual_range.hpp"
#include "cell_columns.hpp"
#include "deltah_calculation.hpp"
#include "cell_motion.hpp"
#include "stateless_rng.hpp"
#include <chrono>

using std::chrono::high_resolution_clock;

bool migration_visual_site_empty(const VisualRange &Visual_range, int x, int y)
{
    return Visual_range.occupied(x, y) == 0;
}

void fill_big_migration_directions(int x1, int y1, const VisualRange &Visual_range, int direction[8])
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

void fill_small_migration_directions(int x1, int y1, const VisualRange &Visual_range, int direction[8])
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

void add_unique_migration_label(int labels[100], int &count, int label)
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

void add_migration_density_site(const VisualRange &Visual_range, int x1, int y1, int local_x, int local_y, int labels[100], int &count)
{
    add_unique_migration_label(labels, count, (int)Visual_range.cell_label(x1 - 5 + local_x, y1 - 5 + local_y));
}

double big_migration_density(int x1, int y1, const VisualRange &Visual_range, int direction_index)
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

double small_migration_density(int x1, int y1, const VisualRange &Visual_range, int direction_index)
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

void apply_big_migration_move(int i, CellStore &cell_array, VisualRange &Visual_range, int x1, int y1, int order)
{
    int row = i - 1;
    auto &x1_values = cell_array.x1();
    auto &x2_values = cell_array.x2();
    auto &x3_values = cell_array.x3();
    auto &x4_values = cell_array.x4();
    auto &y1_values = cell_array.y1();
    auto &y2_values = cell_array.y2();
    auto &y3_values = cell_array.y3();
    auto &y4_values = cell_array.y4();
    auto &ids = cell_array.id();
    auto &migration_directions = cell_array.migration_direction();
    auto &migration_elapsed = cell_array.migration_elapsed();
    long cell_label_1=Visual_range.cell_label(x1,y1);
    long cellstage=Visual_range.stage(x1,y1);
    Visual_range.clear_square(x1,y1);
    int dx = migration_direction_dx(order);
    int dy = migration_direction_dy(order);
    x1_values[row] += dx;
    x2_values[row] += dx;
    x3_values[row] += dx;
    x4_values[row] += dx;
    y1_values[row] += dy;
    y2_values[row] += dy;
    y3_values[row] += dy;
    y4_values[row] += dy;
    Visual_range.write_square((int)x1_values[row], (int)y1_values[row], (long)ids[row], cellstage, cell_label_1);
    migration_directions[row]=order;
    migration_elapsed[row]=0;
}

void apply_small_migration_move(int i, CellStore &cell_array, VisualRange &Visual_range, int x1, int y1, int order)
{
    int row = i - 1;
    auto &x1_values = cell_array.x1();
    auto &y1_values = cell_array.y1();
    auto &ids = cell_array.id();
    auto &migration_directions = cell_array.migration_direction();
    auto &migration_elapsed = cell_array.migration_elapsed();
    long cell_label_1=Visual_range.cell_label(x1,y1);
    long cellstage=Visual_range.stage(x1,y1);
    Visual_range.clear_site(x1,y1);
    int dx = migration_direction_dx(order);
    int dy = migration_direction_dy(order);
    x1_values[row] += dx;
    y1_values[row] += dy;
    Visual_range.write_site((int)x1_values[row], (int)y1_values[row], (long)ids[row], cellstage, cell_label_1);
    migration_directions[row]=order;
    migration_elapsed[row]=0;
}

void migration(int i, double deltah, CellStore &cell_array, VisualRange &Visual_range, double &migration_judgement, long rng_time_step, long rng_event_base)
{
    (void)deltah;
    int row = i - 1;
    auto &x1_values = cell_array.x1();
    auto &y1_values = cell_array.y1();
    auto &types = cell_array.type();
    auto &stages = cell_array.stage();
    auto &ids = cell_array.id();
    auto &migration_elapsed = cell_array.migration_elapsed();
    auto &migration_directions = cell_array.migration_direction();
    auto &migration_follow_flags = cell_array.migration_follow_flag();
    int x1=(int)x1_values[row];
    int y1=(int)y1_values[row];
    long cell_rng_id = (long)ids[row];
    if (cell_rng_id == 0)
    {
        cell_rng_id = i;
    }
    long rng_event = rng_event_base + ((long)types[row] * 10000) + ((long)stages[row] * 1000);
    int cell_type=(int)types[row];
    switch (cell_type)
    {
        case 1://r cells
        {
            int cell_shape=stages[row];
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
                        int migration_direction=migration_directions[row];
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
                                    migration_follow_flags[row]=1;
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
//                                    //                                    migration_elapsed[row]=0;
//                                    //                                    migration_directions[row]=0;
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
//                                    //                                    migration_elapsed[row]=0;
//                                    //                                    migration_directions[row]=0;
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
//                                    //                                    migration_elapsed[row]=0;
//                                    //                                    migration_directions[row]=0;
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
//                                    //                                    migration_elapsed[row]=0;
//                                    //                                    migration_directions[row]=0;
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
//                                    //                                    migration_elapsed[row]=0;
//                                    //                                    migration_directions[row]=0;
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
//                                    //                                    migration_elapsed[row]=0;
//                                    //                                    migration_directions[row]=0;
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
//                                    //                                    migration_elapsed[row]=0;
//                                    //                                    migration_directions[row]=0;
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
//                                    //                                    migration_elapsed[row]=0;
//                                    //                                    migration_directions[row]=0;
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
                            migration_directions[row]=0;
                            migration_elapsed[row]=0;
                        }
                        delete[] direction1;
                        direction1 = NULL;
                        //                        migration_elapsed[row]=0;
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
                        int migration_direction=(int)migration_directions[row];
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
                                    migration_follow_flags[row]=1;
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
//                                    //                                    migration_elapsed[row]=0;
//                                    //                                    migration_directions[row]=0;
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
//                                    //                                    migration_elapsed[row]=0;
//                                    //                                    migration_directions[row]=0;
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
//                                    //                                    migration_elapsed[row]=0;
//                                    //                                    migration_directions[row]=0;
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
//                                    //                                    migration_elapsed[row]=0;
//                                    //                                    migration_directions[row]=0;
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
//                                    //                                    migration_elapsed[row]=0;
//                                    //                                    migration_directions[row]=0;
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
//                                    //                                    migration_elapsed[row]=0;
//                                    //                                    migration_directions[row]=0;
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
//                                    //                                    migration_elapsed[row]=0;
//                                    //                                    migration_directions[row]=0;
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
//                                    //                                    migration_elapsed[row]=0;
//                                    //                                    migration_directions[row]=0;
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
                            migration_directions[row]=0;
                            migration_elapsed[row]=0;
                        }
                        delete[] direction1;
                        direction1 = NULL;
                        //                        migration_elapsed[row]=0;
                    }
                    break;
                }
            }
            break;
        }
        case 2://K cells
        {
            int cell_shape=(int)stages[row];
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
                            migration_directions[row]=0;
                            migration_elapsed[row]=0;
                        }
                        delete[] direction1;
                        direction1 = NULL;
                        //                        migration_elapsed[row]=0;
                        migration_follow_flags[row]=1;
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
                            migration_directions[row]=0;
                            migration_elapsed[row]=0;
                        }
                        delete[] direction1;
                        direction1 = NULL;
                        //                        migration_elapsed[row]=0;
                        migration_follow_flags[row]=1;
                    }
                    break;
                }
            }
            break;
        }
    }
    migration_judgement=migration_judgement+1;
}
