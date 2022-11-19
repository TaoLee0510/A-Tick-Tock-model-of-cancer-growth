//
//  cell_type_transform.hpp
//  CCDS
//
//  Created by Tao Lee on 11/6/22.
//  Copyright © 2022 Tao Lee. All rights reserved.
//

#ifndef cell_type_transform_hpp
#define cell_type_transform_hpp

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


void cell_type_transform(Array<float, 2> &cell_temp, double beta_distribution_alpha_for_normal_migration,double beta_distribution_beta_for_normal_migration,double migration_rate_K_mean,double uniup_K, double unilow_K,double sigmahatK,double muhatK,int &K_label,int i,Array<int, 3> sub_visual, Array<int,3> Visual_range,Array<float, 2> &cell_array, double beta_distribution_alpha, double beta_distribution_beta, double migration_rate_r_mean,double migration_rate_r_mean_quia,double beta_distribution_expected_for_normal_migration,int &r_label,double K_formation_rate)
{
    double Dr1=density_calculation(i, sub_visual, Visual_range, cell_array);
    double initial_K_growth_rate1;
    double migration_rate_K2;
    const gsl_rng_type *T78;
    gsl_rng *r78;
    gsl_rng_env_setup();
    T78 = gsl_rng_ranlxs0;
    gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    r78 = gsl_rng_alloc(T78);
    if(cell_array(i,14)==0)
    {
        if (cell_array(i,9)==1)
        {
            if (Dr1>=0.5)
            {
                std::random_device r;
                std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
                std::mt19937 RNG(seed);
                const gsl_rng_type *T77;
                gsl_rng *r77;
                gsl_rng_env_setup();
                T77 = gsl_rng_ranlxs0;
                gsl_rng_default_seed = ((unsigned long)(time(NULL)));
                r77 = gsl_rng_alloc(T77);
                double random_uni=gsl_rng_uniform(r77);
                if(random_uni<=K_formation_rate)
                {
                    K_label=K_label+1;
                    cell_temp(1,9)=2;//cell_array type
                    Array<float,2> radom_number(1,100,FortranArray<2>());
                    radom_number=random_uniform(100);
                    double rangK2 = uniup_K - unilow_K;
                    double rand2 = (radom_number(1,50)*rangK2)+unilow_K;
                    initial_K_growth_rate1=gsl_cdf_gaussian_Pinv(rand2, sigmahatK) + muhatK;
                    cell_temp(1,10)=initial_K_growth_rate1;//    $10: inherent growth rate
                    cell_temp(1,11)=initial_K_growth_rate1;// $11: density growth rate
                    migration_rate_K2=gsl_ran_beta(r77,beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_K_mean;
                    cell_temp(1,12)=migration_rate_K2;// $12: inherent migration rate
                    cell_temp(1,13)=cell_array(i,13);// $13: mass absorb rate
                    cell_temp(1,14)=0;// $14: cell_array stage
                    cell_temp(1,15)=K_label;// $15: cell_array index
                    cell_temp(1,22)=1;//    $22: cell_array viability
                    cell_temp(1,23)=0;
                    cell_temp(1,24)=0;
                    cell_array(i,16)=0;//    $16: pass time to next division
                    cell_temp(1,16)=0;//    $16: pass time to next division
                    cell_temp(1,25)=cell_array(i,25);//    $25: migration judgement lables:  0: non_migration  1: migration
                    cell_temp(1,26)=cell_array(i,26);//    $26: migration lasted time
                    cell_temp(1,27)=cell_array(i,27);//    $27: passed time of migration
                    cell_temp(1,28)=migration_rate_K2;//    $28: migration rate
                    
                }
                else
                {
                    cell_temp(1,9)=cell_array(i,9);
                    cell_temp(1,10)=cell_array(i,10);
                    cell_temp(1,11)=cell_array(i,11);
                    
                    
                    double mig=gsl_ran_beta(r77,beta_distribution_alpha,beta_distribution_beta)*migration_rate_r_mean;
                    if (mig<=migration_rate_r_mean_quia)
                    {
                        cell_array(i,12)=migration_rate_r_mean_quia*beta_distribution_expected_for_normal_migration;
                    }
                    else
                    {
                        cell_array(i,12)=mig;
                    }
                    
                    double mig1=gsl_ran_beta(r77,beta_distribution_alpha,beta_distribution_beta)*migration_rate_r_mean;
                    if (mig1<=migration_rate_r_mean_quia)
                    {
                        cell_temp(1,12)=migration_rate_r_mean_quia*beta_distribution_expected_for_normal_migration;
                    }
                    else
                    {
                        cell_temp(1,12)=mig1;
                    }
                    
                    cell_temp(1,13)=cell_array(i,13);
                    cell_temp(1,15)=r_label+1;
                    cell_temp(1,18)=0;
                    cell_temp(1,19)=cell_array(i,19);
                    cell_temp(1,21)=cell_array(i,21);
                    cell_temp(1,22)=cell_array(i,22);
                    cell_temp(1,23)=0;
                    cell_temp(1,24)=0;
                    cell_array(i,16)=0;
                    cell_temp(1,16)=0;
                    cell_temp(1,25)=cell_array(i,25);
                    cell_temp(1,26)=cell_array(i,26);
                    cell_temp(1,27)=cell_array(i,27);
                    cell_temp(1,28)=gsl_ran_beta(r77,beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_r_mean_quia;
                }
            }
            else
            {
                cell_temp(1,9)=cell_array(i,9);
                cell_temp(1,10)=cell_array(i,10);
                cell_temp(1,11)=cell_array(i,11);
                
                
                double mig=gsl_ran_beta(r78,beta_distribution_alpha,beta_distribution_beta)*migration_rate_r_mean;
                if (mig<=migration_rate_r_mean_quia)
                {
                    cell_array(i,12)=migration_rate_r_mean_quia*beta_distribution_expected_for_normal_migration;
                }
                else
                {
                    cell_array(i,12)=mig;
                }
                
                
                cell_temp(1,12)=cell_array(i,12);
                cell_temp(1,13)=cell_array(i,13);
                cell_temp(1,15)=r_label+1;
                cell_temp(1,18)=0;
                cell_temp(1,19)=cell_array(i,19);
                cell_temp(1,21)=cell_array(i,21);
                cell_temp(1,22)=cell_array(i,22);
                cell_temp(1,23)=0;
                cell_temp(1,24)=0;
                cell_array(i,16)=0;
                cell_temp(1,16)=0;
                cell_temp(1,25)=cell_array(i,25);
                cell_temp(1,26)=cell_array(i,26);
                cell_temp(1,27)=cell_array(i,27);
                cell_temp(1,28)=gsl_ran_beta(r78,beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_r_mean_quia;
            }
        }
        else if (cell_array(i,9)==2)
        {
            K_label=K_label+1;
           
            cell_array(i,12)=gsl_ran_beta(r78,beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_K_mean;
            
            cell_temp(1,9)=cell_array(i,9);
            cell_temp(1,10)=cell_array(i,10);
            cell_temp(1,11)=cell_array(i,11);
            
            cell_temp(1,12)=cell_array(i,12)=gsl_ran_beta(r78,beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_K_mean;
            
            cell_temp(1,13)=cell_array(i,13);
            cell_temp(1,15)=K_label;
            cell_temp(1,18)=0;
            cell_temp(1,19)=cell_array(i,19);
            cell_temp(1,21)=cell_array(i,21);
            cell_temp(1,22)=cell_array(i,22);
            cell_temp(1,23)=0;
            cell_temp(1,24)=0;
            cell_array(i,16)=0;
            cell_temp(1,16)=0;
            cell_temp(1,25)=cell_array(i,25);
            cell_temp(1,26)=cell_array(i,26);
            cell_temp(1,27)=cell_array(i,27);
            cell_temp(1,28)=cell_array(i,12);
        }
    }
    else if(cell_array(i,14)==1)
    {
        if (cell_array(i,9)==1)
        {
            if(Dr1>=0.5)
            {
                std::random_device r;
                std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
                std::mt19937 RNG(seed);
                const gsl_rng_type *T77;
                gsl_rng *r77;
                gsl_rng_env_setup();
                T77 = gsl_rng_ranlxs0;
                gsl_rng_default_seed = ((unsigned long)(time(NULL)));
                r77 = gsl_rng_alloc(T77);
                double random_uni=gsl_rng_uniform(r77);
                if(random_uni<=K_formation_rate)
                {
                    K_label=K_label+1;
                    cell_temp(1,9)=2;//cell_array type
                    Array<float,2> radom_number(1,100,FortranArray<2>());
                    radom_number=random_uniform(100);
                    double rangK2 = uniup_K - unilow_K;
                    double rand2 = (radom_number(1,50)*rangK2)+unilow_K;
                    initial_K_growth_rate1=gsl_cdf_gaussian_Pinv(rand2, sigmahatK) + muhatK;
                    cell_temp(1,10)=initial_K_growth_rate1;//    $10: inherent growth rate
                    cell_temp(1,11)=initial_K_growth_rate1;// $11: density growth rate
                    migration_rate_K2=gsl_ran_beta(r77,beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_K_mean;
                    cell_temp(1,12)=migration_rate_K2;// $12: inherent migration rate
                    cell_temp(1,13)=cell_array(i,13);// $13: mass absorb rate
                    cell_temp(1,15)=K_label;// $15: cell_array index
                    cell_temp(1,22)=1;//    $22: cell_array viability
                    cell_temp(1,23)=0;
                    cell_temp(1,24)=0;
                    cell_array(i,16)=0;//    $16: pass time to next division
                    cell_temp(1,16)=0;//    $16: pass time to next division
                    cell_temp(1,25)=cell_array(i,25);//    $25: migration judgement lables:  0: non_migration  1: migration
                    cell_temp(1,26)=cell_array(i,26);//    $26: migration lasted time
                    cell_temp(1,27)=cell_array(i,27);//    $27: passed time of migration
                    cell_temp(1,28)=migration_rate_K2;//    $28: migration rate
                }
                else
                {
                    cell_temp(1,9)=cell_array(i,9);
                    cell_temp(1,10)=cell_array(i,10);
                    cell_temp(1,11)=cell_array(i,11);
                    
                    double mig=gsl_ran_beta(r77,beta_distribution_alpha,beta_distribution_beta)*migration_rate_r_mean;
                    if (mig<=migration_rate_r_mean_quia)
                    {
                        cell_array(i,12)=migration_rate_r_mean_quia*beta_distribution_expected_for_normal_migration;
                    }
                    else
                    {
                        cell_array(i,12)=mig;
                    }
                    
                    double mig1=gsl_ran_beta(r77,beta_distribution_alpha,beta_distribution_beta)*migration_rate_r_mean;
                    if (mig1<=migration_rate_r_mean_quia)
                    {
                        cell_temp(1,12)=migration_rate_r_mean_quia*beta_distribution_expected_for_normal_migration;
                    }
                    else
                    {
                        cell_temp(1,12)=mig1;
                    }
                    
                    cell_temp(1,13)=cell_array(i,13);
                    cell_temp(1,15)=r_label+1;
                    cell_temp(1,18)=cell_array(i,18);
                    cell_temp(1,19)=cell_array(i,19);
                    cell_temp(1,21)=cell_array(i,21);
                    cell_temp(1,22)=cell_array(i,22);
                    cell_temp(1,25)=cell_array(i,25);
                    cell_temp(1,26)=cell_array(i,26);
                    cell_temp(1,27)=cell_array(i,27);
                    cell_temp(1,28)=gsl_ran_beta(r77,beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_r_mean_quia;
                }
            }
            else
            {
                cell_temp(1,9)=cell_array(i,9);
                cell_temp(1,10)=cell_array(i,10);
                cell_temp(1,11)=cell_array(i,11);
                
                double mig=gsl_ran_beta(r78,beta_distribution_alpha,beta_distribution_beta)*migration_rate_r_mean;
                if (mig<=migration_rate_r_mean_quia)
                {
                    cell_array(i,12)=migration_rate_r_mean_quia*beta_distribution_expected_for_normal_migration;
                }
                else
                {
                    cell_array(i,12)=mig;
                }
                
                double mig1=gsl_ran_beta(r78,beta_distribution_alpha,beta_distribution_beta)*migration_rate_r_mean;
                if (mig1<=migration_rate_r_mean_quia)
                {
                    cell_temp(1,12)=migration_rate_r_mean_quia*beta_distribution_expected_for_normal_migration;
                }
                else
                {
                    cell_temp(1,12)=mig1;
                }
                
                cell_temp(1,13)=cell_array(i,13);
                cell_temp(1,15)=r_label+1;
                cell_temp(1,18)=cell_array(i,18);
                cell_temp(1,19)=cell_array(i,19);
                cell_temp(1,21)=cell_array(i,21);
                cell_temp(1,22)=cell_array(i,22);
                cell_temp(1,25)=cell_array(i,25);
                cell_temp(1,26)=cell_array(i,26);
                cell_temp(1,27)=cell_array(i,27);
                cell_temp(1,28)=gsl_ran_beta(r78,beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_r_mean_quia;
            }
        }
        else if (cell_array(i,9)==2)
        {
            K_label=K_label+1;
         
            cell_array(i,12)=gsl_ran_beta(r78,beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_K_mean;
            
            cell_temp(1,9)=cell_array(i,9);
            cell_temp(1,10)=cell_array(i,10);
            cell_temp(1,11)=cell_array(i,11);
            
            cell_temp(1,12)=gsl_ran_beta(r78,beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_K_mean;
            cell_temp(1,13)=cell_array(i,13);
            cell_temp(1,15)=K_label;
            cell_temp(1,18)=0;
            cell_temp(1,19)=cell_array(i,19);
            cell_temp(1,21)=cell_array(i,21);
            cell_temp(1,22)=cell_array(i,22);
            cell_temp(1,23)=0;
            cell_temp(1,24)=0;
            cell_array(i,16)=0;
            cell_temp(1,16)=0;
            cell_temp(1,25)=cell_array(i,25);
            cell_temp(1,26)=cell_array(i,26);
            cell_temp(1,27)=cell_array(i,27);
            cell_temp(1,28)=cell_array(i,12);
        }
    }
}


#endif /* cell_type_transform_hpp */
