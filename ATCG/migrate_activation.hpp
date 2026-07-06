//
//  migrate_activation.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef migrate_activation_hpp
#define migrate_activation_hpp

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
#include "density_calculation.hpp"
#include "deltah_calculation.hpp"
#include "stateless_rng.hpp"

inline void migrate_activation(Array<double, 2> &cell_array, double bunderD, Array<long, 3> &sub_visual, const Array<long,3> &Visual_range,double migration_time_range, double migration_rate_r_mean_quia, double beta_distribution_alpha_for_normal_migration, double beta_distribution_beta_for_normal_migration, double beta_distribution_alpha_mig_time, double beta_distribution_beta_mig_time, int DDM, long rng_time_step)
{
    (void)migration_time_range;
    int i;
    int C= cell_array.rows();
    for(i=1;i<=C;i++)
    {
        long cell_rng_id = (long)cell_array(i,cell_col::kId);
        if (cell_rng_id == 0)
        {
            cell_rng_id = i;
        }
        long rng_event = 100;
        if (DDM==1)
        {
            double Dr=density_calculation(i, sub_visual, Visual_range, cell_array);
            int dudgement=0;
            if (Dr<bunderD)
            {
                dudgement=1;
            }
            switch (dudgement)
            {
                case 0:{
//                    double probability_of_time=1/migration_time_range;
                    cell_array(i,cell_col::kMigrationActive)=1;
                    double inherent_migration_speed=(double)cell_array(i,cell_col::kMigrationRateBase);
                    cell_array(i,cell_col::kMigrationRate)=inherent_migration_speed;
                    cell_array(i,cell_col::kMigrationDuration)=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha_mig_time,beta_distribution_beta_mig_time)*(cell_array(i,cell_col::kDivisionTime)-cell_array(i,cell_col::kDivisionElapsed));
//                    double mig_low=migration_rate_r_mean_quia*10;
////                    double mig_threshold=migration_time_range*3/4;
//                    if (mig<=mig_low)
//                    {
//                    }
//                    else
//                    {
//                        cell_array(i,26)=mig;
//                    }
                    break;
                }
                case 1:
                {
                    int cell_type=(int)cell_array(i,cell_col::kType);
                    switch (cell_type)
                    {
                        case 1:
                        {
                            cell_array(i,cell_col::kMigrationRate)=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_r_mean_quia;
                            break;
                        }
                        default:
                        {
                            double inherent_migration_speed=cell_array(i,cell_col::kMigrationRateBase);
                            cell_array(i,cell_col::kMigrationRate)=inherent_migration_speed;
                            break;
                        }
                    }
                    break;
                }
            }
        }
        else
        {
            int cell_type=(int)cell_array(i,cell_col::kType);
            switch (cell_type)
            {
                case 1:
                {
                    cell_array(i,cell_col::kMigrationRate)=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_r_mean_quia;
                    break;
                }
                default:
                {
                    double inherent_migration_speed=cell_array(i,cell_col::kMigrationRateBase);
                    cell_array(i,cell_col::kMigrationRate)=inherent_migration_speed;
                    break;
                }
            }
        }
    }
}

inline void migrate_activation(CellStore &cells, double bunderD, Array<long, 3> &sub_visual, const Array<long,3> &Visual_range,double migration_time_range, double migration_rate_r_mean_quia, double beta_distribution_alpha_for_normal_migration, double beta_distribution_beta_for_normal_migration, double beta_distribution_alpha_mig_time, double beta_distribution_beta_mig_time, int DDM, long rng_time_step)
{
    (void)migration_time_range;
    CellStore::Column &ids = cells.id();
    CellStore::Column &types = cells.type();
    CellStore::Column &migration_rate_base = cells.migration_rate_base();
    CellStore::Column &division_elapsed = cells.division_elapsed();
    CellStore::Column &division_time = cells.division_time();
    CellStore::Column &migration_active = cells.migration_active();
    CellStore::Column &migration_duration = cells.migration_duration();
    CellStore::Column &migration_rate = cells.migration_rate();

    int C= cells.rows();
    for(int i=1;i<=C;i++)
    {
        int row = i - 1;
        long cell_rng_id = (long)ids[row];
        if (cell_rng_id == 0)
        {
            cell_rng_id = i;
        }
        long rng_event = 100;
        if (DDM==1)
        {
            double Dr=density_calculation(i, sub_visual, Visual_range, cells);
            int dudgement=0;
            if (Dr<bunderD)
            {
                dudgement=1;
            }
            switch (dudgement)
            {
                case 0:{
                    migration_active[row]=1;
                    double inherent_migration_speed=(double)migration_rate_base[row];
                    migration_rate[row]=inherent_migration_speed;
                    migration_duration[row]=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha_mig_time,beta_distribution_beta_mig_time)*(division_time[row]-division_elapsed[row]);
                    break;
                }
                case 1:
                {
                    int cell_type=(int)types[row];
                    switch (cell_type)
                    {
                        case 1:
                        {
                            migration_rate[row]=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_r_mean_quia;
                            break;
                        }
                        default:
                        {
                            double inherent_migration_speed=migration_rate_base[row];
                            migration_rate[row]=inherent_migration_speed;
                            break;
                        }
                    }
                    break;
                }
            }
        }
        else
        {
            int cell_type=(int)types[row];
            switch (cell_type)
            {
                case 1:
                {
                    migration_rate[row]=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_r_mean_quia;
                    break;
                }
                default:
                {
                    double inherent_migration_speed=migration_rate_base[row];
                    migration_rate[row]=inherent_migration_speed;
                    break;
                }
            }
        }
    }
}
#endif /* migrate_activation_hpp */
