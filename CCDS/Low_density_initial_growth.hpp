//
//  Low_density_initial_growth.hpp
//  CCDS
//
//  Created by Taolee on 2019/4/13.
//  Copyright © 2019 Tao Lee. All rights reserved.
//

#ifndef Low_density_initial_growth_hpp
#define Low_density_initial_growth_hpp

#include <stdio.h>
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
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "random_uniform.hpp"
#include "outer_corr.hpp"
#include "outer_cell_count.hpp"
#include "deltah_calculation.hpp"
#include "outer_initiation_array_low_density.hpp"
#include "out_initiation_visualrange.hpp"
#include "inner_count_low_density.hpp"
#include "inner_initiation_array.hpp"
#include "density_growth_rate_calculation_1.hpp"
#include "sort_time_per_generation.hpp"
#include "sort_cell_type.hpp"
#include "stage_convert.hpp"
#include "death_judgement.hpp"
#include "save_data.hpp"
#include "sort_time_division.hpp"
#include "random_migration.hpp"
#include "migration.hpp"
#include "division_low.hpp"
#include "migrate_activation.hpp"
#include "density_calculation.hpp"
#include "deltah_recalculation.hpp"


void Low_density_initial_growth(int Visual_range_x, int Visual_range_y, double R0, double R1, double mix_ratio_initial, float alpha, float beta, int DDM, int chemotaxis, double migration_rate_r_mean, double migration_rate_r_mean_quia, double migration_rate_K_mean, double deathjudge, double time_interval, int utralsmall, int allpng,int free_living)
{
    ///////////////////////////////////////////////////////// parameters definition//////////////////////////////////////////////////////////////////////////////////////
    double r_limit=0;
    double K_limit=0;
    double carrying_capacity_r=0;
    double carrying_capacity_K=0;
    double beta_distribution_alpha=0.05;
    double beta_distribution_expected=0.1;
    double beta_distribution_beta=(beta_distribution_alpha*(1-beta_distribution_expected))/beta_distribution_expected;/////////////*************************
    double beta_distribution_alpha_for_normal_migration=5;
    double beta_distribution_expected_for_normal_migration=0.5;
    double beta_distribution_beta_for_normal_migration=(beta_distribution_alpha_for_normal_migration*(1-beta_distribution_expected_for_normal_migration))/beta_distribution_expected_for_normal_migration;
    double death_time_range_r=48;
    double death_time_range_K=120;
    double deltah;
    double migration_time_range=72;///*********************
    
    double muhatr=1.1832;
    double sigmahatr=0.2441;
    double muhatK=0.6832;
    double sigmahatK=0.3764;
    double min_growth_rate_r=1.0722619;
    double min_growth_rate_K=0.33963482;
    double max_growth_rate_r=1.3171805;
    double max_growth_rate_K=0.99505180;
    
    
    double bunderD=0.9;///0.7
    int N0=0;
    int N01=0;
    int N0r;
    int N0K;
    
    int Vx=Visual_range_x+200;
    int Vy=Visual_range_y+200;
    int MMR=0;
    int MMR1=0;
    int MMR2=0;
    int cell_label=(Visual_range_x+200)*(Visual_range_y+200)+1;
    
    int K_label=100;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    const gsl_rng_type *T00;
    gsl_rng *r00;
    gsl_rng_env_setup();
    T00 = gsl_rng_ranlxs0;
    gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    r00 = gsl_rng_alloc(T00);
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 RNG(seed);
    std::uniform_real_distribution<> dis(0.05, 0.95);
    //////////////////////////////////////////////////////////////array definition///////////////////////////////////////////////////////////////////////
    Range all = Range::all();
    Array<int,3> Visual_range(Vx,Vy,4,FortranArray<3>());
    Visual_range(all,all,all)=0;
    Array<float,2> cell_array(1,28,FortranArray<2>());
    cell_array=0;
    //    $9: cell_array type
    //    $10: inherent growth rate
    //    $11: density growth rate
    //    $12: inherent migration rate
    //    $13: mass absorb rate
    //    $14: cell_array stage
    //    $15: cell_array index
    //    $16: pass time to next division
    //    $17: time for a generation
    //    $18: death time
    //    $19: time passed to death
    //    $20: pass time to next migrate
    //    $21: ones migrate time
    //    $22: cell_array viability
    //    $23: migration lable: 1=follow  0=initial
    //    $24: last migration direction
    //    $25: migration judgement lables:  0: non_migration  1: migration
    //    $26: migration lasted time
    //    $27: passed time of migration
    //    $28: migration rate
    Array<float,2> cell_array1(1,28,FortranArray<2>());
    cell_array1=0;
    Array<float,2> cell_array_temp(1,28,FortranArray<2>());
    cell_array_temp=0;
    Array<float, 2> cell_array_temp1(1,28,FortranArray<2>());
    cell_array_temp1=0;
    Array<float,2> cell_array2(1,28,FortranArray<2>());
    cell_array2=0;
    Array<int,2> cor_big(1,4,FortranArray<2>());
    cor_big=0;
    Array<int, 2> area_square(1,10,FortranArray<2>());
    area_square=0;
    Array<int, 2> sub_area_square(1,10,FortranArray<2>());
    sub_area_square=0;
    Array<int, 2> cor_small(1,3,FortranArray<2>());
    cor_small=0;
    Array<int, 2> area_square_s(9,9,FortranArray<2>());
    area_square_s=0;
    Array<int, 2>  sub_area_square_s(9,9,FortranArray<2>());
    sub_area_square_s=0;
    Array<int,2> Visual_range_1(10,10,FortranArray<2>());
    Visual_range_1=0;
    Array<int,2> Visual_range_2(10,10,FortranArray<2>());
    Visual_range_2=0;
    Array<int, 2> cor_big_1(1,16,FortranArray<2>());
    cor_big_1=0;
    Array<int, 2> cor_big_1_change_shape(1,16,FortranArray<2>());
    cor_big_1_change_shape=0;
    Array<int, 2> cor_small_1(1,8,FortranArray<2>());
    cor_small_1=0;
    Array<int, 2> proliferation_loci(1,4,FortranArray<2>());
    proliferation_loci=0;
    Array<float, 2> cell_temp(1,28,FortranArray<2>());
    cell_temp=0;
    Array<int, 3> sub_visual(3,3,4,FortranArray<3>());
    sub_visual=0;
    Array<float,2> cell_array0(1,28,FortranArray<2>());
    cell_array0=0;
    Array<float,2> cell_array_out(1,28,FortranArray<2>());
    cell_array_out=0;
    Array<float,2> cell_array_inner(1,28,FortranArray<2>());
    cell_array_out=0;
    Array<int,2> A(Visual_range_x/2,Visual_range_y/2,FortranArray<2>());
    A=0;
    //int NNy=4*(int)R1*(int)R1;
    int NNy=2000*2000;
    Array<double,2> colorspace(NNy,4,FortranArray<2>());
    colorspace=0;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (int i=1;i<=NNy;i++)
    {
        for (int j=1;j<=4;j++)
        {
            if (j==1)
            {
                colorspace(i,j)=i;
            }
            else
            {
                colorspace(i,j)=dis(RNG);
            }
        }
    }
    //////////////////////////////////////////////////////////// parameters calculation///////////////////////////////////////////////////////////////////////////////////
    double si_r;
    double si_K;
    double limit_r;
    double limit_K;
    double lambda_r = 0.0;
    double lambda_K = 0.0;
    double Cr = 0.0;
    double CK = 0.0;
    if (utralsmall==1)
    {
        carrying_capacity_r=36;
        carrying_capacity_K=42;
        r_limit=carrying_capacity_r*0.5;
        K_limit=carrying_capacity_K*0.5;
    }
    else if (utralsmall==0)
    {
        carrying_capacity_r=31;
        carrying_capacity_K=36;
        r_limit=carrying_capacity_r*0.5;
        K_limit=carrying_capacity_K*0.5;
    }
    si_r=carrying_capacity_r*carrying_capacity_r*log(carrying_capacity_r);
    si_K=carrying_capacity_K*carrying_capacity_K*log(carrying_capacity_K);
    limit_r=r_limit*r_limit*log(r_limit);
    limit_K=K_limit*K_limit*log(K_limit);
    lambda_r=carrying_capacity_r/(si_r-limit_r);
    lambda_K=carrying_capacity_K/(si_K-limit_K);
    Cr=1-(lambda_r*carrying_capacity_r*log(carrying_capacity_r));
    CK=1-(lambda_K*carrying_capacity_K*log(carrying_capacity_K));
    
    int N00=0;
    double unilow_K=0;
    double uniup_K=0;
    if(free_living == 0)
    {
        ////////////////////////////////////////////////////////////////////////outer initiation////////////////////////////////////////////////////////////////
        outer_corr(Visual_range_x,Visual_range_y,R0,R1,A);
        outer_cell_count(Visual_range_x,Visual_range_y,N0,R0,R1);
        N0r=N0*mix_ratio_initial;
        N0K=N0*(1-mix_ratio_initial);
        double migration_rate_r[N0r];
        double migration_rate_K[N0K];
        double unilow_r=gsl_cdf_gaussian_P(min_growth_rate_r-muhatr, sigmahatr );
        double uniup_r=gsl_cdf_gaussian_P(max_growth_rate_r-muhatr, sigmahatr );
        unilow_K=gsl_cdf_gaussian_P(min_growth_rate_K-muhatK, sigmahatK );
        uniup_K=gsl_cdf_gaussian_P(max_growth_rate_K-muhatK, sigmahatK );
        Array<float,2> radom_number(1,N0,FortranArray<2>());
        for (int x=1;x<=N0r;x++)
        {
            double mig=gsl_ran_beta(r00,beta_distribution_alpha,beta_distribution_beta)*migration_rate_r_mean;
            if (mig<=migration_rate_r_mean_quia)
            {
                migration_rate_r[x-1]=migration_rate_r_mean_quia*beta_distribution_expected_for_normal_migration;
            }
            else
            {
                migration_rate_r[x-1]=mig;
            }
        }
        for (int x=1;x<=N0K;x++)
        {
            migration_rate_K[x-1]=gsl_ran_beta(r00,beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_K_mean;
        }
        //////////* Outer deltah calculation*///////////////////////////////////////////
        double deltah1=deltah_calculation(N0, migration_rate_r,N0r,MMR1,DDM);
        /////////////////////////Initiation////////////////////////////////
        cell_array_out.resize(N0,28);
        cell_array_out=0;
        cell_array_out=outer_initiation_array_low_density(N0, Visual_range_x, Visual_range_y, A, uniup_r, unilow_r, sigmahatr, muhatr, uniup_K, unilow_K, sigmahatK, muhatK, N0r, N0K, migration_rate_r, migration_rate_K);
        Visual_range=outer_initiation_visualrange(cell_array_out, N0, Vx, Vy, cell_label);
        ////////////////////////////////////////////////////////////////////////* Inner cells initiation*////////////////////////////////////////////////////////////////
        N01=inner_count_low_density(Visual_range_x, Visual_range_y, Visual_range, N01, R1);
        int NN=N0+N01;
        //////////////////////////*Parameters calculation*////////////////////
        int N0r1=N01*mix_ratio_initial;
        int N0K1=N01*(1-mix_ratio_initial);
        double migration_rate_r1[N0r1];
        double migration_rate_K1[N0K1];
        double unilow_r1=gsl_cdf_gaussian_P(min_growth_rate_r-muhatr, sigmahatr );
        double uniup_r1=gsl_cdf_gaussian_P(max_growth_rate_r-muhatr, sigmahatr );
        double unilow_K1=gsl_cdf_gaussian_P(min_growth_rate_K-muhatK, sigmahatK );
        double uniup_K1=gsl_cdf_gaussian_P(max_growth_rate_K-muhatK, sigmahatK );
        
        /////////////////////////*Migration seepd generation*/////////////////////
        Array<float,2> radom_number1(1,N01,FortranArray<2>());
        for (int x=1;x<=N0r1;x++)
        {
            double mig=gsl_ran_beta(r00,beta_distribution_alpha,beta_distribution_beta)*migration_rate_r_mean;
            if (mig<=migration_rate_r_mean_quia)
            {
                migration_rate_r1[x-1]=migration_rate_r_mean_quia*beta_distribution_expected_for_normal_migration;
            }
            else
            {
                migration_rate_r1[x-1]=mig;
            }
        }
        for (int x=1;x<=N0K1;x++)
        {
            migration_rate_K1[x-1]=gsl_ran_beta(r00,beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_K_mean;
        }
        //////////* Inner deltah calculation*///////////////////////////////////////////
        double deltah2=deltah_calculation(N01, migration_rate_r1, N0r1,MMR2,DDM);
        /////////////////////////*initiation*////////////////////////////////
        cell_array_inner.resize(N01,28);
        cell_array_inner=0;
        inner_initiation_array(N0, N01, R0,Visual_range_x, Visual_range_y,cell_array_inner, Visual_range, uniup_r1, unilow_r1, sigmahatr, muhatr, uniup_K1, unilow_K1, sigmahatK, muhatK, N0r1, N0K1, migration_rate_r1, migration_rate_K1);
        for (int x=1; x<=N01; x++)
        {
            int x1 = cell_array_inner(x,1);
            int y1 = cell_array_inner(x,5);
            int cell_array_index=cell_array_inner(x,15);
            int cell_array_stage=cell_array_inner(x,14);
            Visual_range(x1,y1,1)=1;
            Visual_range(x1,y1,2)=cell_array_index;
            Visual_range(x1,y1,3)=cell_array_stage;
            Visual_range(x1,y1,4)=cell_label;
            cell_label=cell_label+1;
        }
        //////////////////////////////////////////////////////////////////////////cell_array combine////////////////////////////////////////////////////////////
        cell_array0.resize(NN,28);
        cell_array0=0;
        for (int i=1;i<=NN;i++)
        {
            if (i<=N0)
            {
                cell_array0(i,all)=cell_array_out(i,all);
            }
            else
            {
                cell_array0(i,all)=cell_array_inner(i-N0,all);
            }
        }
        //////////////////////////////////////////////////////////////////////////parameters renew////////////////////////////////////////////////////////////
        if (deltah1<deltah2)
        {
            deltah=deltah1;
        }
        else
        {
            deltah=deltah2;
        }
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if (DDM==1 && deltah>0.1) //(0.008)
        {
            cout << "error: Initiation false, simulation aborted" <<endl;//////////////////////////////////////// error mesage ////////////
            exit(0);
        }
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if (MMR1<MMR2)
        {
            MMR=MMR2;
        }
        else
        {
            MMR=MMR1;
        }
        N00=N0;
        N0=NN;
        N0r=N0r+N0r1;
        N0K=N0K+N0K1;
    }
    else{
        A(Visual_range_x/4,Visual_range_y/4)=1;
        N0=1;
        N0r=N0*mix_ratio_initial;
        N0K=1;
        double migration_rate_r[N0r];
        double migration_rate_K[N0K];
        double unilow_r=gsl_cdf_gaussian_P(min_growth_rate_r-muhatr, sigmahatr );
        double uniup_r=gsl_cdf_gaussian_P(max_growth_rate_r-muhatr, sigmahatr );
        unilow_K=gsl_cdf_gaussian_P(min_growth_rate_K-muhatK, sigmahatK );
        uniup_K=gsl_cdf_gaussian_P(max_growth_rate_K-muhatK, sigmahatK );
        Array<float,2> radom_number(1,N0,FortranArray<2>());
        for (int x=1;x<=N0r;x++)
        {
            double mig=gsl_ran_beta(r00,beta_distribution_alpha,beta_distribution_beta)*migration_rate_r_mean;
            if (mig<=migration_rate_r_mean_quia)
            {
                migration_rate_r[x-1]=migration_rate_r_mean_quia*beta_distribution_expected_for_normal_migration;
            }
            else
            {
                migration_rate_r[x-1]=mig;
            }
            for (int x=1;x<=N0K;x++)
            {
                migration_rate_K[x-1]=gsl_ran_beta(r00,beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_K_mean;
            }
        }
        double deltah1=deltah_calculation(N0, migration_rate_r,N0r,MMR1,DDM);
        /////////////////////////Initiation////////////////////////////////
        cell_array_out.resize(N0,28);
        cell_array_out=0;
        cell_array_out=outer_initiation_array(N0, Visual_range_x, Visual_range_y, A, uniup_r, unilow_r, sigmahatr, muhatr, uniup_K, unilow_K, sigmahatK, muhatK, N0r, N0K, migration_rate_r, migration_rate_K);
        Visual_range=outer_initiation_visualrange(cell_array_out, N0, Vx, Vy, cell_label);
        
        cell_array0.resize(N0,28);
        cell_array0=0;
        cell_array0(all,all)=cell_array_out(all,all);
        
        deltah=deltah1;
        N00=N0;
        MMR=MMR1;
    }
    
    //////////////////////////////////////////////////////////////////////////output parameters/////////////////////////////////////////////////////
    char dirname [100] = {'\0'};
    sprintf(dirname, "mkdir ./a_%.1f_b_%.1f",alpha,beta);
    system(dirname);
    char dirname1 [100] = {'\0'};
    sprintf(dirname1, "mkdir ./a_%.1f_b_%.1f_pics",alpha,beta);
    system(dirname1);
    char dirname3 [100] = {'\0'};
    sprintf(dirname3, "mkdir ./a_%.1f_b_%.1f_clonepics",alpha,beta);
    system(dirname3);
    if (allpng==1)
    {
        char dirname2 [100] = {'\0'};
        sprintf(dirname2, "mkdir ./a_%.1f_b_%.1f_picsall",alpha,beta);
        system(dirname2);
        char dirname4 [100] = {'\0'};
        sprintf(dirname4, "mkdir ./a_%.1f_b_%.1f_clonepicsall",alpha,beta);
        system(dirname4);
    }
    char filedir [100] = {'\0'};
    sprintf(filedir, "./Parameters.txt");
    FILE * fid1;
    fid1=fopen (filedir,"w+");
    fprintf(fid1, "%s %s %lf\n" ,"R0", "=", R0);
    fprintf(fid1, "%s %s %lf\n" ,"R1", "=", R1);
    fprintf(fid1, "%s %s %d\n" ,"Visual_range_x", "=", Visual_range_x);
    fprintf(fid1, "%s %s %d\n" ,"Visual_range_y", "=", Visual_range_y);
    fprintf(fid1, "%s %s %lf\n" ,"mix_ratio_initial", "=", mix_ratio_initial);
    fprintf(fid1, "%s %s %lf\n" ,"time_interval", "=", time_interval);
    fprintf(fid1, "%s %s %lf\n" ,"r_limit", "=", r_limit);
    fprintf(fid1, "%s %s %lf\n" ,"K_limit", "=", K_limit);
    fprintf(fid1, "%s %s %lf\n" ,"carrying_capacity_r", "=", carrying_capacity_r);
    fprintf(fid1, "%s %s %lf\n" ,"carrying_capacity_K", "=", carrying_capacity_K);
    fprintf(fid1, "%s %s %lf\n" ,"alpha", "=", alpha);
    fprintf(fid1, "%s %s %lf\n" ,"beta", "=", beta);
    fprintf(fid1, "%s %s %lf\n" ,"beta_distribution_alpha", "=", beta_distribution_alpha);
    fprintf(fid1, "%s %s %lf\n" ,"beta_distribution_expected", "=", beta_distribution_expected);
    fprintf(fid1, "%s %s %lf\n" ,"migration_rate_r_mean", "=", migration_rate_r_mean);
    fprintf(fid1, "%s %s %lf\n" ,"beta_distribution_alpha_for_normal_migration", "=", beta_distribution_alpha_for_normal_migration);
    fprintf(fid1, "%s %s %lf\n" ,"beta_distribution_expected_for_normal_migration", "=", beta_distribution_expected_for_normal_migration);
    fprintf(fid1, "%s %s %lf\n" ,"migration_rate_r_mean_quia", "=", migration_rate_r_mean_quia);
    fprintf(fid1, "%s %s %lf\n" ,"migration_rate_K_mean", "=", migration_rate_K_mean);
    fprintf(fid1, "%s %s %d\n" ,"chemotaxis", "=", chemotaxis);
    fprintf(fid1, "%s %s %lf\n" ,"death_time_range_r", "=", death_time_range_r);
    fprintf(fid1, "%s %s %lf\n" ,"death_time_range_K", "=", death_time_range_K);
    fprintf(fid1, "%s %s %lf\n" ,"deltah", "=", deltah);
    fprintf(fid1, "%s %s %lf\n" ,"muhatr", "=", muhatr);
    fprintf(fid1, "%s %s %lf\n" ,"sigmahatr", "=", sigmahatr);
    fprintf(fid1, "%s %s %lf\n" ,"muhatK", "=", muhatK);
    fprintf(fid1, "%s %s %lf\n" ,"sigmahatK", "=", sigmahatK);
    fprintf(fid1, "%s %s %lf\n" ,"migration_time_range", "=", migration_time_range);
    fprintf(fid1, "%s %s %lf\n" ,"bunderD", "=", bunderD);
    fprintf(fid1, "%s %s %d\n" ,"DDM", "=", DDM);
    fprintf(fid1, "%s %s %d\n" ,"utralsmall_stage", "=", utralsmall);
    fprintf(fid1, "%s %s %d\n" ,"output_all_PNGs", "=", allpng);
    fprintf(fid1, "%s %s %d\n" ,"free_living", "=", free_living);
    fclose(fid1);
    ////////////////////////////////////////////////////////////////////migration and proliferation//////////////////////////////////////////////////////////////
    cell_array.resize(N0,28);
    cell_array=0;
    cell_array(all,all)=cell_array0(all,all);
    migrate_activation(cell_array, bunderD, sub_visual, Visual_range,migration_time_range, migration_rate_r_mean_quia,beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration, DDM);
    density_growth_rate_calculation_1(Visual_range_x, Visual_range_y, N00, N01, r_limit, K_limit, lambda_r, lambda_K, alpha, beta, carrying_capacity_r, carrying_capacity_K, Cr, CK,death_time_range_r,death_time_range_K,cell_array, sub_visual, Visual_range);
    sort_time_per_generation(cell_array, cell_array1);
    double h=0;
    int T=0;
    double migration_judgement=0;
    for (int  H=0; H<1000000000; H++)
    {
        const gsl_rng_type *T10;
        gsl_rng *r10;
        gsl_rng_env_setup();
        T10 = gsl_rng_ranlxs0;
        gsl_rng_default_seed = ((unsigned long)(time(NULL)));
        r10 = gsl_rng_alloc(T10);
        if (h>time_interval)
        {
            break;
        }
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        time_t rawtime;
        struct tm * timeinfo;
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        cout << "h = " << h << "     time : " << asctime (timeinfo) << endl;
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        death_judgement(Visual_range_x, Visual_range_y, N00, N01, r_limit, K_limit, lambda_r, lambda_K, alpha, beta, carrying_capacity_r, carrying_capacity_K, Cr, CK, death_time_range_r,death_time_range_K, deltah, h, cell_array, cell_array_temp, sub_visual, Visual_range, deathjudge);
        sort_cell_type(cell_array, cell_array1);
        stage_convert(Visual_range_x, Visual_range_y, cell_array, Visual_range, cell_label,utralsmall);
        deltah_recalculation(deltah, cell_array, MMR, DDM);
        save_data(Visual_range_x, Visual_range_y, N0, N00, N01, MMR, H, T, alpha, beta, cell_array,migration_judgement, deltah, colorspace,DDM, allpng,free_living);
        sort_time_division(cell_array,cell_array1);
        int C1=cell_array.rows();
        for (int i=C1; i>=1; i--)
        {
            if (cell_array(i,1)>=100 && cell_array(i,5) >=100 && cell_array(i,1)<=Visual_range_x+100 && cell_array(i,5)<=Visual_range_y+100)
            {
                if (cell_array(i,11)>deathjudge)
                {
                    if (cell_array(i,16)<cell_array(i,17))
                    {
                        double expected_division_time=24/cell_array(i,11);
                        double undividing_time=0.9*expected_division_time;
                        if (cell_array(i,16)<=undividing_time)
                        {
                            if (cell_array(i,25)==0)
                            {
                                if (DDM==1)
                                {
                                    double Dr=density_calculation(i, sub_visual, Visual_range, cell_array);
                                    if (Dr>=bunderD)
                                    {
                                        cell_array(i,25)=1;
                                        double inherent_migration_speed=cell_array(i,12);
                                        cell_array(i,28)=inherent_migration_speed;
                                        cell_array(i,26)=cell_array(i,17);
                                    }
                                    cell_array(i,21)=1/cell_array(i,28);
                                }
                                else
                                {
                                    cell_array(i,21)=1/cell_array(i,28);
                                }
                                if (cell_array(i,20)>=cell_array(i,21))
                                {
                                    if (chemotaxis==0)
                                    {
                                        random_migration(i, deltah, cell_array, Visual_range, cor_big, area_square, sub_area_square, cor_small, area_square_s, sub_area_square_s,migration_judgement);
                                    }
                                    else
                                    {
                                        migration(i, deltah,cell_array, Visual_range, cor_big, area_square, sub_area_square, cor_small, area_square_s, sub_area_square_s,migration_judgement);
                                    }
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
                                        if (chemotaxis==0)
                                        {
                                            random_migration(i, deltah, cell_array, Visual_range, cor_big, area_square, sub_area_square, cor_small, area_square_s, sub_area_square_s,migration_judgement);
                                        }
                                        else
                                        {
                                            migration(i, deltah,cell_array, Visual_range, cor_big, area_square, sub_area_square, cor_small, area_square_s, sub_area_square_s,migration_judgement);
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
                                        if (chemotaxis==0)
                                        {
                                            random_migration(i, deltah, cell_array, Visual_range, cor_big, area_square, sub_area_square, cor_small, area_square_s, sub_area_square_s,migration_judgement);
                                        }
                                        else
                                        {
                                            migration(i, deltah,cell_array, Visual_range, cor_big, area_square, sub_area_square, cor_small, area_square_s, sub_area_square_s,migration_judgement);
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
                            
                            
                        }
                        cell_array(i,16)=cell_array(i,16)+deltah;
                    }
                    else
                    {
                        division_low(i, max_growth_rate_r, max_growth_rate_K, cell_array, cell_array_temp, Visual_range, cor_big_1, cor_big_1_change_shape, cor_small_1, proliferation_loci, cell_temp,cell_label,deltah,utralsmall,beta_distribution_alpha_for_normal_migration, beta_distribution_beta_for_normal_migration,migration_rate_K_mean, uniup_K,unilow_K,sigmahatK,muhatK,K_label,sub_visual,beta_distribution_alpha, beta_distribution_beta, migration_rate_r_mean, migration_rate_r_mean_quia, beta_distribution_expected_for_normal_migration);
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
                            if (chemotaxis==0)
                            {
                                random_migration(i, deltah, cell_array, Visual_range, cor_big, area_square, sub_area_square, cor_small, area_square_s, sub_area_square_s,migration_judgement);
                            }
                            else
                            {
                                migration(i, deltah,cell_array, Visual_range, cor_big, area_square, sub_area_square, cor_small, area_square_s, sub_area_square_s,migration_judgement);
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
        h=h+deltah;
        gsl_rng_free(r10);
    }
}
#endif /* Low_density_initial_growth_hpp */
