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
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "density_calculation.hpp"
#include "deltah_calculation.hpp"

void migrate_activation(Array<float, 2> &cell_array, double bunderD, Array<int, 3> sub_visual, Array<int,3> Visual_range,double migration_time_range, double migration_rate_r_mean_quia, double beta_distribution_alpha_for_normal_migration, double beta_distribution_beta_for_normal_migration, int DDM)
{
    const gsl_rng_type *T10;
    gsl_rng *r10;
    gsl_rng_env_setup();
    T10 = gsl_rng_ranlxs0;
    gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    r10 = gsl_rng_alloc(T10);
    int i;
    int C= cell_array.rows();
    for(i=1;i<=C;i++)
    {
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
                    double probability_of_time=1/migration_time_range;
                    cell_array(i,25)=1;
                    double inherent_migration_speed=(double)cell_array(i,12);
                    cell_array(i,28)=inherent_migration_speed;
                    cell_array(i,26)=gsl_ran_geometric(r10,probability_of_time);
                    break;
                }
                case 1:
                {
                    int cell_type=(int)cell_array(i,9);
                    switch (cell_type)
                    {
                        case 1:
                        {
                            cell_array(i,28)=gsl_ran_beta(r10,beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_r_mean_quia;
                            break;
                        }
                        default:
                        {
                            double inherent_migration_speed=cell_array(i,12);
                            cell_array(i,28)=inherent_migration_speed;
                            break;
                        }
                    }
                    break;
                }
            }
        }
        else
        {
            int cell_type=(int)cell_array(i,9);
            switch (cell_type)
            {
                case 1:
                {
                    cell_array(i,28)=gsl_ran_beta(r10,beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_r_mean_quia;
                    break;
                }
                default:
                {
                    double inherent_migration_speed=cell_array(i,12);
                    cell_array(i,28)=inherent_migration_speed;
                    break;
                }
            }
        }
    }
    gsl_rng_free(r10);
}
#endif /* migrate_activation_hpp */
