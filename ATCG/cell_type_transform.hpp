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
#define BZ_THREADSAFE
#define BZ_THREADSAFE_USE_OPENMP
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "density_calculation.hpp"
#include "cell_store.hpp"
#include "visual_range.hpp"
#include "stateless_rng.hpp"
#include <chrono>

using std::chrono::high_resolution_clock;

template <typename CellArray>
inline void cell_type_transform(CellRowBuffer &cell_temp, double beta_distribution_alpha_for_normal_migration,double beta_distribution_beta_for_normal_migration,double migration_rate_K_mean,double uniup_K, double unilow_K,double sigmahatK,double muhatK,long &K_label,int i,const VisualRange &Visual_range,CellArray &cell_array, double beta_distribution_alpha, double beta_distribution_beta, double migration_rate_r_mean,double migration_rate_r_mean_quia,double beta_distribution_expected_for_normal_migration,long &r_label,double K_formation_rate, long cell_rng_id, long rng_time_step, long &rng_event)
{
    int row = i - 1;
    auto &types = cell_array.type();
    auto &growth_rates = cell_array.growth_rate();
    auto &density_growth_rates = cell_array.density_growth_rate();
    auto &migration_rate_bases = cell_array.migration_rate_base();
    auto &random_labels = cell_array.random_label();
    auto &stages = cell_array.stage();
    auto &division_elapsed = cell_array.division_elapsed();
    auto &death_times = cell_array.death_time();
    auto &death_elapsed = cell_array.death_elapsed();
    auto &migration_intervals = cell_array.migration_interval();
    auto &viability = cell_array.viability();
    auto &migration_active = cell_array.migration_active();
    auto &migration_duration = cell_array.migration_duration();
    auto &migration_passed = cell_array.migration_passed();
    double Dr1=density_calculation(i, Visual_range, cell_array);
    double initial_K_growth_rate1;
    double migration_rate_K2;
    if (cell_rng_id == 0)
    {
        cell_rng_id = i;
    }
    if(stages[row]==0)
    {
        if (types[row]==1)
        {
            if (Dr1>=0.5)
            {
                double random_uni=stateless_uniform(cell_rng_id, rng_time_step, rng_event++);
                if(random_uni<=K_formation_rate)
                {
                    K_label=K_label+1;
                    cell_temp(1,9)=2;//cell_array type
                    long uniform_event_base = rng_event;
                    rng_event += 101;
                    double rangK2 = uniup_K - unilow_K;
                    double rand2 = (stateless_uniform(cell_rng_id, rng_time_step, uniform_event_base + 50)*rangK2)+unilow_K;
                    initial_K_growth_rate1=gsl_cdf_gaussian_Pinv(rand2, sigmahatK) + muhatK;
                    cell_temp(1,10)=initial_K_growth_rate1;//    $10: inherent growth rate
                    cell_temp(1,11)=initial_K_growth_rate1;// $11: density growth rate
                    migration_rate_K2=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_K_mean;
                    cell_temp(1,12)=migration_rate_K2;// $12: inherent migration rate
//                    cell_temp(1,13)=random_labels[row];// $13: mass absorb rate
                    cell_temp(1,14)=0;// $14: cell_array stage
                    cell_temp(1,15)=K_label;// $15: cell_array index
                    cell_temp(1,22)=1;//    $22: cell_array viability
                    cell_temp(1,23)=0;
                    cell_temp(1,24)=0;
                    division_elapsed[row]=0;//    $16: pass time to next division
                    cell_temp(1,16)=0;//    $16: pass time to next division
                    cell_temp(1,25)=migration_active[row];//    $25: migration judgement lables:  0: non_migration  1: migration
                    cell_temp(1,26)=migration_duration[row];//    $26: migration lasted time
                    cell_temp(1,27)=migration_passed[row];//    $27: passed time of migration
                    cell_temp(1,28)=migration_rate_K2;//    $28: migration rate
                    
                }
                else
                {
                    cell_temp(1,9)=types[row];
                    cell_temp(1,10)=growth_rates[row];
                    cell_temp(1,11)=density_growth_rates[row];
                    
                    
                    double mig=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha,beta_distribution_beta)*migration_rate_r_mean;
                    if (mig<=migration_rate_r_mean_quia)
                    {
                        migration_rate_bases[row]=migration_rate_r_mean_quia*beta_distribution_expected_for_normal_migration;
                    }
                    else
                    {
                        migration_rate_bases[row]=mig;
                    }
                    
                    double mig1=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha,beta_distribution_beta)*migration_rate_r_mean;
                    if (mig1<=migration_rate_r_mean_quia)
                    {
                        cell_temp(1,12)=migration_rate_r_mean_quia*beta_distribution_expected_for_normal_migration;
                    }
                    else
                    {
                        cell_temp(1,12)=mig1;
                    }
                    
//                    cell_temp(1,13)=random_labels[row];
                    cell_temp(1,15)=r_label+1;
                    cell_temp(1,18)=0;
                    cell_temp(1,19)=death_elapsed[row];
                    cell_temp(1,21)=migration_intervals[row];
                    cell_temp(1,22)=viability[row];
                    cell_temp(1,23)=0;
                    cell_temp(1,24)=0;
                    division_elapsed[row]=0;
                    cell_temp(1,16)=0;
                    cell_temp(1,25)=migration_active[row];
                    cell_temp(1,26)=migration_duration[row];
                    cell_temp(1,27)=migration_passed[row];
                    cell_temp(1,28)=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_r_mean_quia;
                }
            }
            else
            {
                cell_temp(1,9)=types[row];
                cell_temp(1,10)=growth_rates[row];
                cell_temp(1,11)=density_growth_rates[row];
                
                
                double mig=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha,beta_distribution_beta)*migration_rate_r_mean;
                if (mig<=migration_rate_r_mean_quia)
                {
                    migration_rate_bases[row]=migration_rate_r_mean_quia*beta_distribution_expected_for_normal_migration;
                }
                else
                {
                    migration_rate_bases[row]=mig;
                }
                
                
                cell_temp(1,12)=migration_rate_bases[row];
//                cell_temp(1,13)=random_labels[row];
                cell_temp(1,15)=r_label+1;
                cell_temp(1,18)=0;
                cell_temp(1,19)=death_elapsed[row];
                cell_temp(1,21)=migration_intervals[row];
                cell_temp(1,22)=viability[row];
                cell_temp(1,23)=0;
                cell_temp(1,24)=0;
                division_elapsed[row]=0;
                cell_temp(1,16)=0;
                cell_temp(1,25)=migration_active[row];
                cell_temp(1,26)=migration_duration[row];
                cell_temp(1,27)=migration_passed[row];
                cell_temp(1,28)=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_r_mean_quia;
            }
        }
        else if (types[row]==2)
        {
            K_label=K_label+1;
           
            migration_rate_bases[row]=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_K_mean;
            
            cell_temp(1,9)=types[row];
            cell_temp(1,10)=growth_rates[row];
            cell_temp(1,11)=density_growth_rates[row];
            
            cell_temp(1,12)=migration_rate_bases[row]=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_K_mean;
            
//            cell_temp(1,13)=random_labels[row];
            cell_temp(1,15)=K_label;
            cell_temp(1,18)=0;
            cell_temp(1,19)=death_elapsed[row];
            cell_temp(1,21)=migration_intervals[row];
            cell_temp(1,22)=viability[row];
            cell_temp(1,23)=0;
            cell_temp(1,24)=0;
            division_elapsed[row]=0;
            cell_temp(1,16)=0;
            cell_temp(1,25)=migration_active[row];
            cell_temp(1,26)=migration_duration[row];
            cell_temp(1,27)=migration_passed[row];
            cell_temp(1,28)=migration_rate_bases[row];
        }
    }
    else if(stages[row]==1)
    {
        if (types[row]==1)
        {
            if(Dr1>=0.5)
            {
                double random_uni=stateless_uniform(cell_rng_id, rng_time_step, rng_event++);
                if(random_uni<=K_formation_rate)
                {
                    K_label=K_label+1;
                    cell_temp(1,9)=2;//cell_array type
                    long uniform_event_base = rng_event;
                    rng_event += 101;
                    double rangK2 = uniup_K - unilow_K;
                    double rand2 = (stateless_uniform(cell_rng_id, rng_time_step, uniform_event_base + 50)*rangK2)+unilow_K;
                    initial_K_growth_rate1=gsl_cdf_gaussian_Pinv(rand2, sigmahatK) + muhatK;
                    cell_temp(1,10)=initial_K_growth_rate1;//    $10: inherent growth rate
                    cell_temp(1,11)=initial_K_growth_rate1;// $11: density growth rate
                    migration_rate_K2=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_K_mean;
                    cell_temp(1,12)=migration_rate_K2;// $12: inherent migration rate
//                    cell_temp(1,13)=random_labels[row];// $13: mass absorb rate
                    cell_temp(1,15)=K_label;// $15: cell_array index
                    cell_temp(1,22)=1;//    $22: cell_array viability
                    cell_temp(1,23)=0;
                    cell_temp(1,24)=0;
                    division_elapsed[row]=0;//    $16: pass time to next division
                    cell_temp(1,16)=0;//    $16: pass time to next division
                    cell_temp(1,25)=migration_active[row];//    $25: migration judgement lables:  0: non_migration  1: migration
                    cell_temp(1,26)=migration_duration[row];//    $26: migration lasted time
                    cell_temp(1,27)=migration_passed[row];//    $27: passed time of migration
                    cell_temp(1,28)=migration_rate_K2;//    $28: migration rate
                }
                else
                {
                    cell_temp(1,9)=types[row];
                    cell_temp(1,10)=growth_rates[row];
                    cell_temp(1,11)=density_growth_rates[row];
                    
                    double mig=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha,beta_distribution_beta)*migration_rate_r_mean;
                    if (mig<=migration_rate_r_mean_quia)
                    {
                        migration_rate_bases[row]=migration_rate_r_mean_quia*beta_distribution_expected_for_normal_migration;
                    }
                    else
                    {
                        migration_rate_bases[row]=mig;
                    }
                    
                    double mig1=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha,beta_distribution_beta)*migration_rate_r_mean;
                    if (mig1<=migration_rate_r_mean_quia)
                    {
                        cell_temp(1,12)=migration_rate_r_mean_quia*beta_distribution_expected_for_normal_migration;
                    }
                    else
                    {
                        cell_temp(1,12)=mig1;
                    }
                    
//                    cell_temp(1,13)=random_labels[row];
                    cell_temp(1,15)=r_label+1;
                    cell_temp(1,18)=death_times[row];
                    cell_temp(1,19)=death_elapsed[row];
                    cell_temp(1,21)=migration_intervals[row];
                    cell_temp(1,22)=viability[row];
                    cell_temp(1,25)=migration_active[row];
                    cell_temp(1,26)=migration_duration[row];
                    cell_temp(1,27)=migration_passed[row];
                    cell_temp(1,28)=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_r_mean_quia;
                }
            }
            else
            {
                cell_temp(1,9)=types[row];
                cell_temp(1,10)=growth_rates[row];
                cell_temp(1,11)=density_growth_rates[row];
                
                double mig=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha,beta_distribution_beta)*migration_rate_r_mean;
                if (mig<=migration_rate_r_mean_quia)
                {
                    migration_rate_bases[row]=migration_rate_r_mean_quia*beta_distribution_expected_for_normal_migration;
                }
                else
                {
                    migration_rate_bases[row]=mig;
                }
                
                double mig1=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha,beta_distribution_beta)*migration_rate_r_mean;
                if (mig1<=migration_rate_r_mean_quia)
                {
                    cell_temp(1,12)=migration_rate_r_mean_quia*beta_distribution_expected_for_normal_migration;
                }
                else
                {
                    cell_temp(1,12)=mig1;
                }
                
//                cell_temp(1,13)=random_labels[row];
                cell_temp(1,15)=r_label+1;
                cell_temp(1,18)=death_times[row];
                cell_temp(1,19)=death_elapsed[row];
                cell_temp(1,21)=migration_intervals[row];
                cell_temp(1,22)=viability[row];
                cell_temp(1,25)=migration_active[row];
                cell_temp(1,26)=migration_duration[row];
                cell_temp(1,27)=migration_passed[row];
                cell_temp(1,28)=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_r_mean_quia;
            }
        }
        else if (types[row]==2)
        {
            K_label=K_label+1;
         
            migration_rate_bases[row]=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_K_mean;
            
            cell_temp(1,9)=types[row];
            cell_temp(1,10)=growth_rates[row];
            cell_temp(1,11)=density_growth_rates[row];
            
            cell_temp(1,12)=stateless_beta(cell_rng_id, rng_time_step, rng_event++, beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_K_mean;
//            cell_temp(1,13)=random_labels[row];
            cell_temp(1,15)=K_label;
            cell_temp(1,18)=0;
            cell_temp(1,19)=death_elapsed[row];
            cell_temp(1,21)=migration_intervals[row];
            cell_temp(1,22)=viability[row];
            cell_temp(1,23)=0;
            cell_temp(1,24)=0;
            division_elapsed[row]=0;
            cell_temp(1,16)=0;
            cell_temp(1,25)=migration_active[row];
            cell_temp(1,26)=migration_duration[row];
            cell_temp(1,27)=migration_passed[row];
            cell_temp(1,28)=migration_rate_bases[row];
        }
    }
}


#endif /* cell_type_transform_hpp */
