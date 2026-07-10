//
//  Low_density_initial_growth.hpp
//  CCDS
//
//  Created by Taolee on 2019/4/13.
//  Copyright © 2019 Tao Lee. All rights reserved.
//

#include "Low_density_initial_growth.hpp"

#include <gsl/gsl_cdf.h>

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
#include "cell_columns.hpp"
#include "cell_store.hpp"
#include "visual_range.hpp"
#include "outer_corr.hpp"
#include "outer_cell_count.hpp"
#include "deltah_calculation.hpp"
#include "outer_initiation_array_low_density.hpp"
#include "out_initiation_visualrange.hpp"
#include "inner_count_low_density.hpp"
#include "inner_initiation_array.hpp"
#include "density_growth_rate_calculation_1.hpp"
#include "stage_convert.hpp"
#include "death_judgement.hpp"
#include "save_data.hpp"
#include "random_migration.hpp"
#include "migration.hpp"
#include "stateless_rng.hpp"
#include "division.hpp"
#include "migrate_activation.hpp"
#include "density_calculation.hpp"
#include "deltah_recalculation.hpp"
#include "int_grid.hpp"
#include <omp.h>

using std::cout;
using std::endl;


void Low_density_initial_growth(int Visual_range_x, int Visual_range_y, double R0, double R1, double mix_ratio_initial, double alpha, double beta, int DDM, int chemotaxis, double migration_rate_r_mean, double migration_rate_r_mean_quia, double migration_rate_K_mean, double deathjudge, double time_interval, int utralsmall, int allpng,double bunderD,double beta_distribution_alpha, double beta_distribution_expected, double beta_distribution_alpha_mig_time, double beta_distribution_expected_mig_time,int threads)
{
    time_t raw_initial_time;
    struct tm * initial_time;
    time ( &raw_initial_time );
    initial_time = localtime ( &raw_initial_time );
    ///////////////////////////////////////////////////////// parameters definition//////////////////////////////////////////////////////////////////////////////////////
    double r_limit=0;
    double K_limit=0;
    double carrying_capacity_r=0;
    double carrying_capacity_K=0;
    double beta_distribution_beta=(beta_distribution_alpha*(1-beta_distribution_expected))/beta_distribution_expected;/////////////*************************
    double beta_distribution_beta_mig_time=(beta_distribution_alpha_mig_time*(1-beta_distribution_expected_mig_time))/beta_distribution_expected_mig_time;
    double beta_distribution_alpha_for_normal_migration=5;
    double beta_distribution_expected_for_normal_migration=0.5;
    double beta_distribution_beta_for_normal_migration=(beta_distribution_alpha_for_normal_migration*(1-beta_distribution_expected_for_normal_migration))/beta_distribution_expected_for_normal_migration;
    double death_time_range_r=48;
    double death_time_range_K=120;
    double deltah;
    double migration_time_range=24;///*********************

    double muhatr=1.1832;
    double sigmahatr=0.2441;
    double muhatK=0.6832;
    double sigmahatK=0.3764;
    double min_growth_rate_r=1.0722619;
    double min_growth_rate_K=0.33963482;
    double max_growth_rate_r=1.3171805;
    double max_growth_rate_K=0.99505180;


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

    int Col=cell_col::kStandardColumnCount;
    int borderx=Visual_range_x+100;
    int bordery=Visual_range_y+100;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    const long rng_context = 20004;
    //////////////////////////////////////////////////////////////array definition///////////////////////////////////////////////////////////////////////
    VisualRange Visual_range(Vx,Vy);
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
    //    $23: last migration direction
    //    $24: migration lable: 1=follow  0=initial
    //    $25: migration judgement lables:  0: non_migration  1: migration
    //    $26: migration lasted time
    //    $27: passed time of migration
    //    $28: migration rate
    CellRowBuffer cell_temp(1, Col);
    IntGrid A(Visual_range_x/2, Visual_range_y/2);
    A=0;
    int NNy=Visual_range_x*Visual_range_y;
    ColorSpace colorspace(NNy, 4);
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
                colorspace(i,j)=0.05 + (0.90 * stateless_uniform(rng_context, 0, 10000 + i * 10 + j));
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
    ////////////////////////////////////////////////////////////////////////outer initiation////////////////////////////////////////////////////////////////
    outer_corr(Visual_range_x,Visual_range_y,R0,R1,A);
    outer_cell_count(Visual_range_x,Visual_range_y,N0,R0,R1);
    N0r=N0*mix_ratio_initial;
    N0K=N0-N0r;
    double migration_rate_r[N0r];
    double migration_rate_K[N0K];
    double unilow_r=gsl_cdf_gaussian_P(min_growth_rate_r-muhatr, sigmahatr );
    double uniup_r=gsl_cdf_gaussian_P(max_growth_rate_r-muhatr, sigmahatr );
    double unilow_K=gsl_cdf_gaussian_P(min_growth_rate_K-muhatK, sigmahatK );
    double uniup_K=gsl_cdf_gaussian_P(max_growth_rate_K-muhatK, sigmahatK );
    for (int x=1;x<=N0r;x++)
    {
        double mig=stateless_beta(rng_context, 0, 20000 + x, beta_distribution_alpha,beta_distribution_beta)*migration_rate_r_mean;
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
        migration_rate_K[x-1]=stateless_beta(rng_context, 0, 21000 + x, beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_K_mean;
    }
    /////////////////////////Initiation////////////////////////////////
    CellStore cells = outer_initiation_low_density_cell_store(N0, Visual_range_x, Visual_range_y, A, uniup_r, unilow_r, sigmahatr, muhatr, uniup_K, unilow_K, sigmahatK, muhatK, N0r, N0K, migration_rate_r, migration_rate_K);
    N0=cells.rows();
    N0r=0;
    N0K=0;
    for (int x=1; x<=N0; x++)
    {
        if (cells.type()[x - 1]==1)
        {
            N0r++;
        }
        else if (cells.type()[x - 1]==2)
        {
            N0K++;
        }
    }
    //////////* Outer deltah calculation*///////////////////////////////////////////
    double deltah1=deltah_calculation(N0, migration_rate_r,N0r,MMR1,DDM);
    Visual_range=outer_initiation_visualrange(cells, N0, Vx, Vy, cell_label);
    ////////////////////////////////////////////////////////////////////////* Inner cells initiation*////////////////////////////////////////////////////////////////
    N01=inner_count_low_density(Visual_range_x, Visual_range_y, Visual_range, N01, R1);
    int NN=N0+N01;
    //////////////////////////*Parameters calculation*////////////////////
    int N0r1=N01*mix_ratio_initial;
    int N0K1=N01-N0r1;
    double migration_rate_r1[N0r1];
    double migration_rate_K1[N0K1];
    double unilow_r1=gsl_cdf_gaussian_P(min_growth_rate_r-muhatr, sigmahatr );
    double uniup_r1=gsl_cdf_gaussian_P(max_growth_rate_r-muhatr, sigmahatr );
    double unilow_K1=gsl_cdf_gaussian_P(min_growth_rate_K-muhatK, sigmahatK );
    double uniup_K1=gsl_cdf_gaussian_P(max_growth_rate_K-muhatK, sigmahatK );

    /////////////////////////*Migration seepd generation*/////////////////////
    for (int x=1;x<=N0r1;x++)
    {
        double mig=stateless_beta(rng_context, 0, 22000 + x, beta_distribution_alpha,beta_distribution_beta)*migration_rate_r_mean;
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
        migration_rate_K1[x-1]=stateless_beta(rng_context, 0, 23000 + x, beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_K_mean;
    }
    //////////* Inner deltah calculation*///////////////////////////////////////////
    double deltah2=deltah_calculation(N01, migration_rate_r1, N0r1,MMR2,DDM);
    /////////////////////////*initiation*////////////////////////////////
    CellStore inner_cells = inner_initiation_cell_store(N0, N01, R1+5,Visual_range_x, Visual_range_y, Visual_range, uniup_r1, unilow_r1, sigmahatr, muhatr, uniup_K1, unilow_K1, sigmahatK, muhatK, N0r1, N0K1, migration_rate_r1, migration_rate_K1, Col);
    for (int x=1; x<=N01; x++)
    {
        int row = x - 1;
        int x1 = inner_cells.x1()[row];
        int y1 = inner_cells.y1()[row];
        int cell_array_index=inner_cells.id()[row];
        int cell_array_stage=inner_cells.stage()[row];
        Visual_range.write_site(x1, y1, cell_array_index, cell_array_stage, cell_label);
        cell_label=cell_label+1;
    }
    //////////////////////////////////////////////////////////////////////////cell_array combine////////////////////////////////////////////////////////////
    cells.append_from(inner_cells);
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


    //////////////////////////////////////////////////////////////////////////output parameters/////////////////////////////////////////////////////
    char dirname [100] = {'\0'};
    snprintf(dirname, sizeof(dirname), "mkdir ./a_%.1f_b_%.1f",alpha,beta);
    system(dirname);
    char dirname1 [100] = {'\0'};
    snprintf(dirname1, sizeof(dirname1), "mkdir ./a_%.1f_b_%.1f_pics",alpha,beta);
    system(dirname1);
    char dirname3 [100] = {'\0'};
    snprintf(dirname3, sizeof(dirname3), "mkdir ./a_%.1f_b_%.1f_clonepics",alpha,beta);
    system(dirname3);
    if (allpng==1)
    {
        char dirname2 [100] = {'\0'};
        snprintf(dirname2, sizeof(dirname2), "mkdir ./a_%.1f_b_%.1f_picsall",alpha,beta);
        system(dirname2);
        char dirname4 [100] = {'\0'};
        snprintf(dirname4, sizeof(dirname4), "mkdir ./a_%.1f_b_%.1f_clonepicsall",alpha,beta);
        system(dirname4);
    }
    char filedir [100] = {'\0'};
    snprintf(filedir, sizeof(filedir), "./Parameters.txt");
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
    fprintf(fid1, "%s %s %lf\n" ,"beta_distribution_alpha", "=", beta_distribution_alpha);
    fprintf(fid1, "%s %s %lf\n" ,"beta_distribution_expected", "=", beta_distribution_expected);
    fprintf(fid1, "%s %s %lf\n" ,"beta_distribution_alpha_mig_time", "=", beta_distribution_alpha_mig_time);
    fprintf(fid1, "%s %s %lf\n" ,"beta_distribution_expected_mig_time", "=", beta_distribution_expected_mig_time);
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
    fclose(fid1);
    ////////////////////////////////////////////////////////////////////migration and proliferation//////////////////////////////////////////////////////////////
    migrate_activation(cells, bunderD, Visual_range,migration_time_range, migration_rate_r_mean_quia,beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration, beta_distribution_alpha_mig_time, beta_distribution_beta_mig_time,DDM, 0);
    density_growth_rate_calculation_1(Visual_range_x, Visual_range_y, N00, N01, r_limit, K_limit, lambda_r, lambda_K, alpha, beta, carrying_capacity_r, carrying_capacity_K, Cr, CK,death_time_range_r,death_time_range_K,cells, Visual_range, 0);
    cells.sort_by_column(cell_col::kDivisionTime);///sort time per generation
    double h=0;
    int T=0;
    double migration_judgement=0;
    for (int  H=0; H<1000000000; H++)
    {
        if (h>time_interval)
        {
            break;
        }
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        time_t rawtime;
        struct tm * timeinfo;
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );

        double input_seconds=difftime(rawtime,raw_initial_time);
        double seconds,minutes,hours,days;
        int seconds02,minutes02, hours02,days02;
        days = input_seconds / 60 / 60 / 24;
        days02=(int)days;

        hours = input_seconds / 60 / 60;
        hours02 = hours - 24 * (double)days02;

        minutes = input_seconds / 60;
        minutes02 = minutes - (60 * (double)hours02)- (24*60* (double)days02);

        seconds = (24*60*60*days02)+(60 * 60 * hours02) + (60 * minutes02);
        seconds02 = input_seconds - seconds;

        double completeness=(h/time_interval)*100;

        cout << "Completeness: "<< completeness << "%" << "\n  ||  h = " << h <<"\n  ||  Cost time (D:H:M:S): "<< days02<<":"<< hours02 <<":"<< minutes02 <<":"<< seconds02 <<endl;
        fflush(stdout);
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///
        int nthreads = 1;
        death_judgement(Visual_range_x, Visual_range_y, N00, N01, r_limit, K_limit, lambda_r, lambda_K, alpha, beta, carrying_capacity_r, carrying_capacity_K, Cr, CK, death_time_range_r,death_time_range_K, deltah, h, cells, Visual_range, deathjudge,Col,nthreads,H);
        cells.sort_by_column(cell_col::kType);///sort cell type
        stage_convert(Visual_range_x, Visual_range_y, cells, Visual_range, cell_label,utralsmall,H);
        deltah_recalculation(deltah, cells, MMR, DDM);
        cells.sort_by_column(cell_col::kDivisionElapsed);///sort time division
        save_data(Visual_range_x, Visual_range_y, N0, N00, N01, MMR, H, T, alpha, beta, cells,migration_judgement, deltah, colorspace,DDM, allpng);
        int C1=cells.rows();
        auto &x1 = cells.x1();
        auto &y1 = cells.y1();
        auto &growth_rate = cells.growth_rate();
        auto &density_growth_rate = cells.density_growth_rate();
        auto &migration_rate_base = cells.migration_rate_base();
        auto &id = cells.id();
        auto &division_elapsed = cells.division_elapsed();
        auto &division_time = cells.division_time();
        auto &death_time = cells.death_time();
        auto &death_elapsed = cells.death_elapsed();
        auto &migration_elapsed = cells.migration_elapsed();
        auto &migration_interval = cells.migration_interval();
        auto &migration_active = cells.migration_active();
        auto &migration_duration = cells.migration_duration();
        auto &migration_passed = cells.migration_passed();
        auto &migration_rate = cells.migration_rate();
        for (int i=C1; i>=1; i--)
        {
            int row = i - 1;
            long cell_rng_id = (long)id[row];
            if (cell_rng_id == 0)
            {
                cell_rng_id = i;
            }
            long rng_event = 100;
            if (x1[row]>=100 && y1[row] >=100 && x1[row]<=borderx && y1[row]<=bordery)
            {
                if (density_growth_rate[row]>deathjudge)
                {
                    if (division_elapsed[row]<division_time[row])
                    {
                        double expected_division_time=24/density_growth_rate[row];
                        double undividing_time=0.9*expected_division_time;
                        if (division_elapsed[row]<=undividing_time)
                        {
                            if (migration_active[row]==0)
                            {
                                if (DDM==1)
                                {
                                    double Dr=density_calculation(i, Visual_range, cells);
                                    if (Dr>=bunderD)
                                    {
                                        migration_active[row]=1;
                                        double inherent_migration_speed=migration_rate_base[row];
                                        migration_rate[row]=inherent_migration_speed;
                                        migration_duration[row]=stateless_beta(cell_rng_id, H, rng_event++,beta_distribution_alpha_mig_time,beta_distribution_beta_mig_time)*(division_time[row]-division_elapsed[row]);
                                    }
                                    migration_interval[row]=1/migration_rate[row];
                                }
                                else
                                {
                                    migration_interval[row]=1/migration_rate[row];
                                }
                                if (migration_elapsed[row]>=migration_interval[row])
                                {
                                    random_migration(i, deltah, cells, Visual_range, migration_judgement, H, 1000 + rng_event++);
                                }
                                else
                                {
                                    migration_elapsed[row]=migration_elapsed[row]+deltah;
                                }
                            }
                            else
                            {
                                if (migration_passed[row]>=migration_duration[row])
                                {
                                    migration_active[row]=0;
                                    migration_duration[row]=0;
                                    migration_passed[row]=0;
                                    migration_rate[row]=stateless_beta(cell_rng_id, H, rng_event++,beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration)*migration_rate_r_mean_quia;
                                    migration_interval[row]=1/migration_rate[row];
                                    if (migration_elapsed[row]>=migration_interval[row])
                                    {
                                        if (chemotaxis==0)
                                        {
                                            random_migration(i, deltah, cells, Visual_range, migration_judgement, H, 1000 + rng_event++);
                                        }
                                        else
                                        {
                                            migration(i, deltah,cells, Visual_range, migration_judgement, H, 2000 + rng_event++);
                                        }
                                    }
                                    else
                                    {
                                        migration_elapsed[row]=migration_elapsed[row]+deltah;
                                    }
                                }
                                else
                                {
                                    if (migration_elapsed[row]>=migration_interval[row])
                                    {
                                        if (chemotaxis==0)
                                        {
                                            random_migration(i, deltah, cells, Visual_range, migration_judgement, H, 1000 + rng_event++);
                                        }
                                        else
                                        {
                                            migration(i, deltah,cells, Visual_range, migration_judgement, H, 2000 + rng_event++);
                                        }
                                    }
                                    else
                                    {
                                        migration_elapsed[row]=migration_elapsed[row]+deltah;
                                        migration_passed[row]=migration_passed[row]+deltah;
                                        migration_rate[row]=migration_rate_base[row];
                                    }
                                }
                            }


                        }
                        division_elapsed[row]=division_elapsed[row]+deltah;
                    }
                    else
                    {
                        division(i, max_growth_rate_r, max_growth_rate_K, cells, Visual_range, cell_temp,cell_label,deltah,utralsmall,Col, H);
                    }
                }
                else
                {
                    double D_time_1=1.5*(24/growth_rate[row]);
                    double D_time_2=0.9*death_time[row];
                    double D_time = 0;
                    if (D_time_1<=D_time_2)
                    {
                        D_time = D_time_1;
                    }
                    else
                    {
                        D_time = D_time_2;
                    }
                    if (death_elapsed[row]<=D_time)
                    {
                        if (migration_elapsed[row]>=migration_interval[row])
                        {
                            if (chemotaxis==0)
                            {
                                random_migration(i, deltah, cells, Visual_range, migration_judgement, H, 1000 + rng_event++);
                            }
                            else
                            {
                                if(migration_active[row]==0)
                                {
                                    random_migration(i, deltah, cells, Visual_range, migration_judgement, H, 1000 + rng_event++);
                                }
                                else
                                {
                                    migration(i, deltah,cells, Visual_range, migration_judgement, H, 2000 + rng_event++);
                                }
                            }
                        }
                        else
                        {
                            migration_elapsed[row]=migration_elapsed[row]+deltah;
                        }
                    }
                    else
                    {
                        migration_elapsed[row]=migration_elapsed[row]+deltah;
                    }
                }
            }
        }
        h=h+deltah;
    }
}
