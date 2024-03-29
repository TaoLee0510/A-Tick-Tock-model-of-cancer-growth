//
//  CellMigrationDivision.hpp
//  CCDS
//
//  Created by Tao Lee on 12/7/22.
//  Copyright © 2022 Tao Lee. All rights reserved.
//

#ifndef CellMigrationDivision_hpp
#define CellMigrationDivision_hpp

#include <stdio.h>
#include <omp.h>
#include <iostream>
#include <time.h>
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
#include <getopt.h>
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
#include "random_uniform.hpp"
#include "outer_corr.hpp"
#include "outer_cell_count.hpp"
#include "deltah_calculation.hpp"
#include "outer_initiation_array.hpp"
#include "out_initiation_visualrange.hpp"
#include "inner_count.hpp"
#include "inner_initiation_array.hpp"
#include "density_growth_rate_calculation_1.hpp"
#include "stage_convert.hpp"
#include "death_judgement.hpp"
#include "save_data_free_living.hpp"
#include "random_migration.hpp"
#include "migration.hpp"
#include "division.hpp"
#include "migrate_activation.hpp"
#include "density_calculation.hpp"
#include "deltah_recalculation.hpp"
#include "sortRow.hpp"
#include <omp.h>
#include <chrono>
using std::chrono::high_resolution_clock;



void CellMigrationDivision(int &DDM, int i, double &deltah,Array<double, 2> &cell_array, Array<long, 3> &Visual_range, Array<int,2> cor_big, Array<int, 2> area_square, Array<int, 2> sub_area_square, Array<int, 2> cor_small, Array<int, 2> area_square_s, Array<int, 2>  sub_area_square_s,double &migration_judgement,double deathjudge, double beta_distribution_alpha_mig_time,double beta_distribution_beta_mig_time,int chemotaxis,double bunderD,Array<long, 3> sub_visual,int borderx,int bordery,double beta_distribution_alpha_for_normal_migration,double migration_rate_r_mean_quia,double beta_distribution_beta_for_normal_migration, double max_growth_rate_r, double max_growth_rate_K,  Array<double,2> cell_array_temp,  Array<int,2> cor_big_1, Array<int, 2> cor_big_1_change_shape, Array<int, 2> cor_small_1, Array<int, 2> proliferation_loci, Array<double, 2> cell_temp,int &cell_label,int utralsmall,int Col)
{
    if (cell_array(i,1)>=100 && cell_array(i,5) >=100 && cell_array(i,1)<=borderx && cell_array(i,5)<=bordery)
    {
//        if (cell_array(i,11)>deathjudge)
        if (cell_array(i,17)>0)
        {
            if (cell_array(i,16)<cell_array(i,17))
            {
                double expected_division_time=24/cell_array(i,11);
                double undividing_time=0.9*expected_division_time;
                if (cell_array(i,16)<=undividing_time)
                {
                    auto start = std::chrono::high_resolution_clock::now();
                    const gsl_rng_type *T10;
                    gsl_rng *r10;
                    gsl_rng_env_setup();
                    T10 = gsl_rng_ranlxs0;
                    auto end = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start);
                    gsl_rng_default_seed = (duration.count());
                    r10 = gsl_rng_alloc(T10);
                    if (cell_array(i,25)==0)
                    {
                        switch (DDM)
                        {
                            case 1:
                            {
                                double Dr=density_calculation(i, sub_visual, Visual_range, cell_array);
                                if (Dr>=bunderD)
                                {
                                    cell_array(i,25)=1;
                                    cell_array(i,28)=cell_array(i,12);
                                    cell_array(i,26)=gsl_ran_beta(r10,beta_distribution_alpha_mig_time,beta_distribution_beta_mig_time)*(cell_array(i,17)-cell_array(i,16));
                                }
                                cell_array(i,21)=1/cell_array(i,28);
                                break;
                            }
                            case 0:
                            {
                                cell_array(i,21)=1/cell_array(i,28);
                                break;
                            }
                        }
                        if (cell_array(i,20)>=cell_array(i,21))
                        {
                            random_migration(i, deltah, cell_array, Visual_range, cor_big, area_square, sub_area_square, cor_small, area_square_s, sub_area_square_s,migration_judgement);
                        }
                        else
                        {
                            cell_array(i,20)=cell_array(i,20)+deltah;
                        }
                    }
                    else
                    {
                        if (cell_array(i,27)>=cell_array(i,26))
                        {
                            cell_array(i,25)=0;
                            cell_array(i,26)=0;
                            cell_array(i,27)=0;
                            cell_array(i,28)=gsl_ran_beta(r10,beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_r_mean_quia;
                            cell_array(i,21)=1/cell_array(i,28);
                            if (cell_array(i,20)>=cell_array(i,21))
                            {
                                switch (chemotaxis)
                                {
                                    case 0:
                                    {
                                        random_migration(i, deltah, cell_array, Visual_range, cor_big, area_square, sub_area_square, cor_small, area_square_s, sub_area_square_s,migration_judgement);
                                        break;
                                    }
                                    case 1:
                                    {
                                        migration(i, deltah,cell_array, Visual_range, cor_big, area_square, sub_area_square, cor_small, area_square_s, sub_area_square_s,migration_judgement);
                                    }
                                }
                            }
                            else
                            {
                                cell_array(i,20)=cell_array(i,20)+deltah;
                            }
                        }
                        else
                        {
                            if (cell_array(i,20)>=cell_array(i,21))
                            {
                                switch (chemotaxis)
                                {
                                    case 0:
                                    {
                                        random_migration(i, deltah, cell_array, Visual_range, cor_big, area_square, sub_area_square, cor_small, area_square_s, sub_area_square_s,migration_judgement);
                                        break;
                                    }
                                    case 1:
                                    {
                                        migration(i, deltah,cell_array, Visual_range, cor_big, area_square, sub_area_square, cor_small, area_square_s, sub_area_square_s,migration_judgement);
                                        break;
                                    }
                                }
                            }
                            else
                            {
                                cell_array(i,20)=cell_array(i,20)+deltah;
                                cell_array(i,27)=cell_array(i,27)+deltah;
                                cell_array(i,28)=cell_array(i,12);
                            }
                        }
                    }
                    gsl_rng_free(r10);
                }
                cell_array(i,16)=cell_array(i,16)+deltah;
            }
            else
            {
                division(i, max_growth_rate_r, max_growth_rate_K, cell_array, cell_array_temp, Visual_range, cor_big_1, cor_big_1_change_shape, cor_small_1, proliferation_loci, cell_temp,cell_label,deltah,utralsmall, Col);
            }
        }
        else
        {
            double D_time_1=1.5*(24/cell_array(i,10));
            double D_time_2=0.9*cell_array(i,18);
            double D_time = 0;
            if (D_time_1<=D_time_2)
            {
                D_time = D_time_1;
            }
            else
            {
                D_time = D_time_2;
            }
            if (cell_array(i,19)<=D_time)
            {
                if (cell_array(i,20)>=cell_array(i,21))
                {
                    switch (chemotaxis)
                    {
                        case 0:
                        {
                            random_migration(i, deltah, cell_array, Visual_range, cor_big, area_square, sub_area_square, cor_small, area_square_s, sub_area_square_s,migration_judgement);
                            break;
                        }
                        case 1:
                        {
                            if(cell_array(i,25)==0)
                            {
                                random_migration(i, deltah, cell_array, Visual_range, cor_big, area_square, sub_area_square, cor_small, area_square_s, sub_area_square_s,migration_judgement);
                            }
                            else
                            {
                                migration(i, deltah,cell_array, Visual_range, cor_big, area_square, sub_area_square, cor_small, area_square_s, sub_area_square_s,migration_judgement);
                            }
                            break;
                        }
                    }
                }
                else
                {
                    cell_array(i,20)=cell_array(i,20)+deltah;
                }
            }
            else
            {
                cell_array(i,20)=cell_array(i,20)+deltah;
            }
        }
    }
}




#endif /* CellMigrationDivision_hpp */
