//
//  CellMigration.hpp
//  CCDS
//
//  Created by Tao Lee on 12/7/22.
//  Copyright © 2022 Tao Lee. All rights reserved.
//

#include "CellMigration.hpp"

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
#include "visual_range.hpp"
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
#include "stateless_rng.hpp"
#include "free_living_division.hpp"
#include "migrate_activation.hpp"
#include "density_calculation.hpp"
#include "deltah_recalculation.hpp"
#include <omp.h>
#include <chrono>

using std::chrono::high_resolution_clock;

void CellMigration(int &DDM, int i, double &deltah,CellStore &cell_array, VisualRange &Visual_range, double &migration_judgement,double deathjudge, double beta_distribution_alpha_mig_time,double beta_distribution_beta_mig_time,int chemotaxis,double bunderD,int borderx,int bordery,double beta_distribution_alpha_for_normal_migration,double migration_rate_r_mean_quia,double beta_distribution_beta_for_normal_migration,long rng_time_step)
{
    long cell_rng_id = (long)cell_array.id()[i - 1];
    if (cell_rng_id == 0)
    {
        cell_rng_id = i;
    }
    long rng_event = 100;
    if (cell_array.x1()[i - 1]>=100 && cell_array.y1()[i - 1] >=100 && cell_array.x1()[i - 1]<=borderx && cell_array.y1()[i - 1]<=bordery)
    {
//        if (cell_array.density_growth_rate()[i - 1]>deathjudge)
        if (cell_array.division_time()[i - 1]>0)
        {
            if (cell_array.division_elapsed()[i - 1]<cell_array.division_time()[i - 1])
            {
//                double undividing_time=21.6/cell_array.density_growth_rate()[i - 1];//24*0.9=21.6
                double undividing_time=cell_array.division_time()[i - 1]*0.9;
                if (cell_array.division_elapsed()[i - 1]<=undividing_time)
                {
                    int migration_label=(int)cell_array.migration_active()[i - 1];
                    switch (migration_label)
                    {
                        case 0:
                        {
                            switch (DDM)
                            {
                                case 1:
                                {
                                    double Dr=density_calculation(i, Visual_range, cell_array);
                                    if (Dr>=bunderD)
                                    {
                                        cell_array.migration_active()[i - 1]=1;
                                        cell_array.migration_rate()[i - 1]=cell_array.migration_rate_base()[i - 1];

                                        cell_array.migration_duration()[i - 1]=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha_mig_time, beta_distribution_beta_mig_time)*(cell_array.division_time()[i - 1]-cell_array.division_elapsed()[i - 1]);
                                    }
                                    cell_array.migration_interval()[i - 1]=1/cell_array.migration_rate()[i - 1];
                                    break;
                                }
                                case 0:
                                {
                                    cell_array.migration_interval()[i - 1]=1/cell_array.migration_rate()[i - 1];
                                    break;
                                }
                            }
                            if (cell_array.migration_elapsed()[i - 1]>=cell_array.migration_interval()[i - 1])
                            {
                                random_migration(i, deltah, cell_array, Visual_range, migration_judgement, rng_time_step, 1000 + rng_event++);
                            }
                            else
                            {
                                cell_array.migration_elapsed()[i - 1]=cell_array.migration_elapsed()[i - 1]+deltah;
                            }
                            break;
                        }
                        case 1:
                        {
                            if (cell_array.migration_passed()[i - 1]>=cell_array.migration_duration()[i - 1])
                            {
                                cell_array.migration_active()[i - 1]=0;
                                cell_array.migration_duration()[i - 1]=0;
                                cell_array.migration_passed()[i - 1]=0;
                                cell_array.migration_rate()[i - 1]=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha_for_normal_migration, beta_distribution_beta_for_normal_migration)*migration_rate_r_mean_quia;
                                cell_array.migration_interval()[i - 1]=1/cell_array.migration_rate()[i - 1];
                                if (cell_array.migration_elapsed()[i - 1]>=cell_array.migration_interval()[i - 1])
                                {
                                    switch (chemotaxis)
                                    {
                                        case 0:
                                        {
                                            random_migration(i, deltah, cell_array, Visual_range, migration_judgement, rng_time_step, 1000 + rng_event++);
                                            break;
                                        }
                                        case 1:
                                        {
                                            migration(i, deltah,cell_array, Visual_range, migration_judgement, rng_time_step, 2000 + rng_event++);
                                            break;
                                        }
                                    }
                                }
                                else
                                {
                                    cell_array.migration_elapsed()[i - 1]=cell_array.migration_elapsed()[i - 1]+deltah;
                                }
                            }
                            else //cell_array.migration_passed()[i - 1]<cell_array.migration_duration()[i - 1]
                            {
                                if (cell_array.migration_elapsed()[i - 1]>=cell_array.migration_interval()[i - 1])
                                {
                                    switch (chemotaxis)
                                    {
                                        case 0:
                                        {
                                            random_migration(i, deltah, cell_array, Visual_range, migration_judgement, rng_time_step, 1000 + rng_event++);
                                            break;
                                        }
                                        case 1:
                                        {
                                            migration(i, deltah,cell_array, Visual_range, migration_judgement, rng_time_step, 2000 + rng_event++);
                                            break;
                                        }
                                    }
                                }
                                else//cell_array.migration_elapsed()[i - 1]<cell_array.migration_interval()[i - 1]
                                {
                                    cell_array.migration_elapsed()[i - 1]=cell_array.migration_elapsed()[i - 1]+deltah;
                                    cell_array.migration_passed()[i - 1]=cell_array.migration_passed()[i - 1]+deltah;
                                    cell_array.migration_rate()[i - 1]=cell_array.migration_rate_base()[i - 1];
                                }
                            }

                            break;
                        }
                    }
                }
//                cell_array.division_elapsed()[i - 1]=cell_array.division_elapsed()[i - 1]+deltah;
            }
        }
        else
        {
            double D_time_1=1.5*(24/cell_array.growth_rate()[i - 1]);
            double D_time_2=0.9*cell_array.death_time()[i - 1];
            double D_time = 0;
            if (D_time_1<=D_time_2)
            {
                D_time = D_time_1;
            }
            else
            {
                D_time = D_time_2;
            }
            if (cell_array.death_elapsed()[i - 1]<=D_time)
            {
                if (cell_array.migration_elapsed()[i - 1]>=cell_array.migration_interval()[i - 1])
                {
                    switch (chemotaxis)
                    {
                        case 0:
                        {
                            random_migration(i, deltah, cell_array, Visual_range, migration_judgement, rng_time_step, 1000 + rng_event++);
                            break;
                        }
                        case 1:
                        {
                            if(cell_array.migration_active()[i - 1]==0)
                            {
                                random_migration(i, deltah, cell_array, Visual_range, migration_judgement, rng_time_step, 1000 + rng_event++);
                            }
                            else
                            {
                                migration(i, deltah,cell_array, Visual_range, migration_judgement, rng_time_step, 2000 + rng_event++);
                            }
                            break;
                        }
                    }
                }
                else
                {
                    cell_array.migration_elapsed()[i - 1]=cell_array.migration_elapsed()[i - 1]+deltah;
                }
            }
            else
            {
                cell_array.migration_elapsed()[i - 1]=cell_array.migration_elapsed()[i - 1]+deltah;
            }
        }
    }
}
