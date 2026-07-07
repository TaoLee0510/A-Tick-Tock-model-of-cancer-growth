//
//  free_living_growth_recovery.hpp
//  ATCG
//
//  Created by Tao Lee on 5/24/24.
//  Copyright © 2024 Tao Lee. All rights reserved.
//

#ifndef free_living_growth_recovery_hpp
#define free_living_growth_recovery_hpp


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
#include "cell_columns.hpp"
#include "cell_store.hpp"
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
#include "free_living_division.hpp"
#include "free_living_division_single_thread.hpp"
#include "migrate_activation.hpp"
#include "density_calculation.hpp"
#include "deltah_recalculation.hpp"
#include <omp.h>
#include "CellMigration.hpp"
#include "CellDivisionSingleCell.hpp"
#include "CellDivision.hpp"
#include "CellMigrationDivisionSingleCell.hpp"
#include "SavePNGSingleCell.hpp"
#include "SaveCellArraySingleCell.hpp"
#include "SaveCellTraceArray.hpp"
#include "SaveAllPNG.hpp"
#include "read_files.hpp"
#include "cell_trace.hpp"
#include "recovery_parameters.hpp"
#include "int_grid.hpp"
#include "CountLines.hpp"

#include <chrono>

using std::chrono::high_resolution_clock;

void free_living_growth_recovery(int Visual_range_x, int Visual_range_y, double R0, double R1, double mix_ratio_initial, double alpha, double beta, int DDM, int chemotaxis, double migration_rate_r_mean, double migration_rate_r_mean_quia, double migration_rate_K_mean, double deathjudge, double time_interval, int utralsmall, int allpng,int Single_cell,double K_formation_rate,double bunderD, double beta_distribution_alpha, double beta_distribution_expected, double beta_distribution_alpha_mig_time, double beta_distribution_expected_mig_time,int threads,int DynamicThreads, string Cell_arry_file, string Cell_trace_arry_file, string Parameters,int Hours)
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
    double migration_time_range=24;///72*********************
    
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
    int cell_label=(Visual_range_x+200)*(Visual_range_y+200)+1;
    
    long r_label=1;
    long K_label=500000000;
    long cell_index=0;
  //int generation=0;
    
    int Col=cell_col::kFreeLivingColumnCount;
    
    int borderx=Visual_range_x+100;
    int bordery=Visual_range_y+100;
    
//    int division_interval=170; //every 51 mins; 170*deltah=0.85 hours, equ 51mins. 0.85/17=0.05. 51mins ~5% cell generation(fast)
    int division_interval=50; //15 mins
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    const long rng_context = 20002;
    //////////////////////////////////////////////////////////////array definition///////////////////////////////////////////////////////////////////////
    Range all = Range::all();
    Array<long,3> Visual_range(Vx,Vy,4,FortranArray<3>());
    Visual_range(all,all,all)=0;
    CellStore cells(Col);
    CellRowBuffer cell_temp(2, Col);
    IntGrid A(Visual_range_x/2, Visual_range_y/2);
    A=0;
    
    
    CellTraceStore cell_trace(1,150);
    cell_trace=0;
    
    CellTraceStore cell_trace_temp(1,150);
    cell_trace_temp=0;
    
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
    double unilow_K=0;
    double uniup_K=0;
    
    
    A(Visual_range_x/4,Visual_range_y/4)=1;
    N0=1;
    N0r=N0*mix_ratio_initial;
    N0K=1;
    unilow_K=gsl_cdf_gaussian_P(min_growth_rate_K-muhatK, sigmahatK );
    uniup_K=gsl_cdf_gaussian_P(max_growth_rate_K-muhatK, sigmahatK );
    ///////////////////////////////////////////////////////// read files //////////////////////////////////////////////////////////////////////////////////////
    
    RecoveryParameters parameters(39);
    parameters=0;
    
    read_file(cells,cell_trace, parameters, Cell_arry_file,Cell_trace_arry_file, Parameters);
    
    

    Vx=parameters(3)+200;
    Vy=parameters(4)+200;
    
    division_interval=parameters(39);
    
    
    int Cell_number = cells.rows();
    
    for (int i=1; i<=Cell_number; ++i)
    {
        if(cells(i,2)==0 && cells(i,6)==0)
        {
            Visual_range(int(cells(i,1)),int(cells(i,5)),1)=1;
            Visual_range(int(cells(i,1)),int(cells(i,5)),2)=cells(i,15);
            Visual_range(int(cells(i,1)),int(cells(i,5)),3)=cells(i,14);
            Visual_range(int(cells(i,1)),int(cells(i,5)),4)=cell_label;
            cell_label=cell_label+1;
        }
        else
        {
            Visual_range(int(cells(i,1)),int(cells(i,5)),1)=1;
            Visual_range(int(cells(i,2)),int(cells(i,6)),1)=1;
            Visual_range(int(cells(i,3)),int(cells(i,7)),1)=1;
            Visual_range(int(cells(i,4)),int(cells(i,8)),1)=1;
            Visual_range(int(cells(i,1)),int(cells(i,5)),2)=cells(i,15);
            Visual_range(int(cells(i,2)),int(cells(i,6)),2)=cells(i,15);
            Visual_range(int(cells(i,3)),int(cells(i,7)),2)=cells(i,15);
            Visual_range(int(cells(i,4)),int(cells(i,8)),2)=cells(i,15);
            Visual_range(int(cells(i,1)),int(cells(i,5)),3)=cells(i,14);
            Visual_range(int(cells(i,2)),int(cells(i,6)),3)=cells(i,14);
            Visual_range(int(cells(i,3)),int(cells(i,7)),3)=cells(i,14);
            Visual_range(int(cells(i,4)),int(cells(i,8)),3)=cells(i,14);
            Visual_range(int(cells(i,1)),int(cells(i,5)),4)=cell_label;
            Visual_range(int(cells(i,2)),int(cells(i,6)),4)=cell_label;
            Visual_range(int(cells(i,3)),int(cells(i,7)),4)=cell_label;
            Visual_range(int(cells(i,4)),int(cells(i,8)),4)=cell_label;
            cell_label=cell_label+1;
        }

    }

    
    /////////////////////////Initiation////////////////////////////////

    deltah=0.005;
    
    N00=N0;
    MMR=1/deltah;

    
    //////////////////////////////////////////////////////////////////////////output parameters/////////////////////////////////////////////////////
    char dirname [100] = {'\0'};
    sprintf(dirname, "mkdir ./a_%.1f_b_%.1f",alpha,beta);
    system(dirname);
    char dirname1 [100] = {'\0'};
    sprintf(dirname1, "mkdir ./a_%.1f_b_%.1f_pics",alpha,beta);
    system(dirname1);
    char dirname2 [100] = {'\0'};
    sprintf(dirname2, "mkdir ./a_%.1f_b_%.1f_CellTrace",alpha,beta);
    system(dirname2);
    if (allpng==1)
    {
        char dirname2 [100] = {'\0'};
        sprintf(dirname2, "mkdir ./a_%.1f_b_%.1f_picsall",alpha,beta);
        system(dirname2);

    }
    char filedir [100] = {'\0'};
    sprintf(filedir, "./Parameters_recovery.txt");
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
    fprintf(fid1, "%s %s %d\n" ,"Single_cell", "=", Single_cell);
    fprintf(fid1, "%s %s %lf\n" ,"K_formation_rate", "=", K_formation_rate);
    fprintf(fid1, "%s %s %d\n" ,"division_interval", "=", division_interval);
    fclose(fid1);
    
    char filedir1 [100] = {'\0'};
    sprintf(filedir1, "./Logs_recovery.txt");
    FILE * fid2;
    fid2=fopen (filedir1,"w+");
    ////////////////////////////////////////////////////////////////////migration and proliferation//////////////////////////////////////////////////////////////
    double h=Hours;
    int T=Hours;
    double migration_judgement=0;
    int MaxThread=omp_get_max_threads();
    int nthreads;
    
    cells.sort_by_column(cell_col::kCellTraceLabel);///sort cell index
    
    cell_index=cells(cells.rows(),cell_col::kCellTraceLabel);

    if(MaxThread<threads)
    {
        cout <<"Warning: custom threads larger than the system maximal threads, reset it as system maximal threads!" <<endl;
        threads=MaxThread;
        nthreads=threads;
    }
    else
    {
        nthreads = DynamicThreads == 1 ? 1 : threads;
    }
    for (int  H=0; H!=1000000000; ++H)
    {
        double start00=omp_get_wtime();
        
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
        
        if (h>time_interval)
        {
            fprintf(fid2, "%s\n" ,"Simulation Completed!");
            fprintf(fid2, "%s %d %s %d %s %d %s %d \n" ,"Cost time (D:H:M:S)  :",days02 ,":", hours02 ,":",minutes02,":", minutes02);
            fprintf(fid2, "%s\n" ,"*********************************************************");
            fclose(fid2);
            break;
        }
        
        
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        int C2=cells.rows();
        if(DynamicThreads==1)
        {
            if (nthreads<threads)
            {
                int CNumber=ceil(C2/5000);
                switch (CNumber)
                {
                    case 0:
                    {
                        nthreads=std::max(1, (int)(threads*0.5));
                        break;
                    }
                    case 1:
                    {
                        nthreads=std::max(1, (int)(threads*0.6));
                        break;
                    }
                    case 2:
                    {
                        nthreads=std::max(1, (int)(threads*0.7));
                        break;
                    }
                    case 3:
                    {
                        nthreads=std::max(1, (int)(threads*0.8));
                        break;
                    }
                    case 4:
                    {
                        nthreads=std::max(1, (int)(threads*0.9));
                        break;
                    }
                    case 5:
                    {
                        nthreads=threads;
                        break;
                    }
                    default:
                    {
                        nthreads=threads;
                        break;
                    }
                }
            }
        }
        double start12(0);
        double end12(0);
        double start14(0);
        double end14(0);
        
        int JU =H%division_interval;


        start12=omp_get_wtime();
        death_judgement(Visual_range_x, Visual_range_y, N00, N01, r_limit, K_limit, lambda_r, lambda_K, alpha, beta, carrying_capacity_r, carrying_capacity_K, Cr, CK, death_time_range_r,death_time_range_K, deltah, h, cells, Visual_range, deathjudge,Col,nthreads,H);
        end12=omp_get_wtime();
        
        start14=omp_get_wtime();
        cells.sort_by_column(cell_col::kType);///sort cell type
        end14=omp_get_wtime();
        stage_convert(Visual_range_x, Visual_range_y, cells, Visual_range, cell_label,utralsmall,H);
        
        
        int C1=cells.rows();
        
        double start04(0);
        double end04(0);
        double start05(0);
        double end05(0);
        double start06(0);
        double end06(0);
        double end07(0);
        double start08(0);
        double end08(0);
        double start09(0);
        double end09(0);
        double start10(0);
        double end10(0);
        double start11(0);
        double end11(0);
        double start13(0);
        double end13(0);
      //double start16=0;
      //double end16=0;
        double start15(0);
        double end15(0);
        double programTimes16(0);
        
        double start17(0);
        double end17(0);
        
        double start18(0);
        double end18(0);
        
        double start19(0);
        double end19(0);
        omp_set_num_threads(nthreads);
        #pragma omp parallel
        {
            #pragma omp sections
            {
                #pragma omp section
                {
                    start05=omp_get_wtime();
                    if (H%MMR==0)
                    {
                        SaveCellTraceArray(T, alpha, beta, cell_trace);
                    }
                    end05=omp_get_wtime();
                }
                #pragma omp section
                {
                    start09=omp_get_wtime();
                    if (H%MMR==0)
                    {
                        SavePNGSingleCell(Visual_range_x, Visual_range_y, T, alpha, beta, cells);
                    }
                    end09=omp_get_wtime();
                }
                #pragma omp section
                {
                    start10=omp_get_wtime();
                    if (H%MMR==0)
                    {
                        SaveCellArraySingleCell(T, alpha, beta, cells , Col);
                    }
                    end10=omp_get_wtime();
                }
                #pragma omp section
                if (allpng==1)
                {
                    start11=omp_get_wtime();
                    SaveAllPNG( Visual_range_x,  Visual_range_y, cells,  H,  T,  alpha,  beta, deltah);
                    end11=omp_get_wtime();
                }
                
            }
        }
        
        if (nthreads != 1)
        {
            if (JU !=0 )
            {
                start15=omp_get_wtime();
                cells.sort_by_column(cell_col::kRandomLabel);///sort random label
                end15=omp_get_wtime();
                
                start04=omp_get_wtime();
                omp_set_num_threads(nthreads);
                #pragma omp parallel for schedule(dynamic)
                {
                    for (int i=C1; i!=0; --i)
                    {
                        CellMigration(DDM, i, deltah,cells, Visual_range, migration_judgement,deathjudge, beta_distribution_alpha_mig_time, beta_distribution_beta_mig_time, chemotaxis, bunderD,borderx,bordery,beta_distribution_alpha_for_normal_migration, migration_rate_r_mean_quia, beta_distribution_beta_for_normal_migration, H);
                        cells(i,16)=cells(i,16)+deltah;// add detalh
                    }
                }
                end04=omp_get_wtime();
                
            }
            else
            {
                start15=omp_get_wtime();
                cells.sort_by_column(cell_col::kRandomLabel);///sort random label
                end15=omp_get_wtime();
                start04=omp_get_wtime();
                omp_set_num_threads(nthreads);
                #pragma omp parallel for schedule(dynamic)
                {
                    for (int i=C1; i!=0; --i)
                    {
                        CellMigration(DDM, i, deltah,cells, Visual_range, migration_judgement,deathjudge, beta_distribution_alpha_mig_time, beta_distribution_beta_mig_time, chemotaxis, bunderD,borderx,bordery,beta_distribution_alpha_for_normal_migration, migration_rate_r_mean_quia, beta_distribution_beta_for_normal_migration, H);
                    }
                }
                end04=omp_get_wtime();
                
                start13=omp_get_wtime();
                cells.sort_by_column(cell_col::kDivisionElapsed);///sort time division
                end13=omp_get_wtime();
                
                double programTimes15 = 0;
                start06=omp_get_wtime();

                CellTraceStore cell_trace_ndcells(1,150);
                cell_trace_ndcells=0;
                CellRowBuffer cell_array_ndcells(1,Col);
                cell_array_ndcells=0;
                int ndcells(0);
                
                for (int i=C1; i!=0; --i)
                {
                    CellDivisionSingleCell(i, max_growth_rate_r,  max_growth_rate_K, cells, Visual_range, cell_temp,cell_label, deltah, utralsmall,  beta_distribution_alpha_for_normal_migration, beta_distribution_beta_for_normal_migration, migration_rate_K_mean, uniup_K, unilow_K,sigmahatK, muhatK, K_label, beta_distribution_alpha,  beta_distribution_beta,  migration_rate_r_mean, migration_rate_r_mean_quia, beta_distribution_expected_for_normal_migration,cell_trace,cell_trace_temp, cell_index, r_label, Col, K_formation_rate, deathjudge, borderx, bordery,fid2,nthreads,programTimes15,cell_trace_ndcells,ndcells,cell_array_ndcells, H);
                    programTimes16 = programTimes16 + programTimes15;
                }
                
//                //////////////// add offsprings to the cell trace array and cell array
                if (ndcells!=0)
                {
                    cell_trace.append_from(cell_trace_ndcells);
                    
                    int current_size_ndcells=cell_array_ndcells.rows();
                    cells.reserve(C1+current_size_ndcells-1);
                    for (int new_cell=2; new_cell<=current_size_ndcells; ++new_cell)
                    {
                        cell_store_append_row_from_array(cells, cell_array_ndcells, new_cell, Col);
                    }
                }
                
                end06=omp_get_wtime();
            }
        }
        else
        {
            start08=omp_get_wtime();
            cells.sort_by_column(cell_col::kDivisionElapsed);///sort time division
            save_data_free_living(Visual_range_x, Visual_range_y, N0, N00, N01, MMR, H, T, alpha, beta, cells,migration_judgement, deltah, colorspace,DDM, allpng, Col, cell_trace);
            for (int i=C1; i!=0; --i)
            {
                CellMigrationDivisionSingleCell(DDM, i, deltah,cells, Visual_range, migration_judgement, max_growth_rate_r,  max_growth_rate_K, cell_temp, cell_label , utralsmall,  beta_distribution_alpha_for_normal_migration, beta_distribution_beta_for_normal_migration, migration_rate_K_mean, uniup_K,  unilow_K, sigmahatK, muhatK, K_label,beta_distribution_alpha,  beta_distribution_beta,  migration_rate_r_mean, migration_rate_r_mean_quia, beta_distribution_expected_for_normal_migration, cell_trace, cell_trace_temp,  cell_index,  r_label, Col, K_formation_rate, deathjudge,  beta_distribution_alpha_mig_time, beta_distribution_beta_mig_time, chemotaxis, bunderD, borderx, bordery,fid2,nthreads, H);
            }
            end08=omp_get_wtime();
        }
        
        
        if (H%MMR==0)
        {
            T++;
        }
        
        end07=omp_get_wtime();
        double programTimes04 = end04 - start04;
        double programTimes05 = end05 - start05;
        double programTimes06 = end06 - start06;
        double programTimes00 = end07 -start00;
        double programTimes07 = end08 -start08;
        double programTimes08 = end08 -start00;
        double programTimes09 = end09 -start09;
        double programTimes10 = end10 -start10;
        double programTimes11 = end11 -start11;
        double programTimes12 = end12 -start12;
        double programTimes13 = end13 -start13;
        double programTimes14 = end14 -start14;
        double programTimes15 = end15 -start15;
        switch (nthreads)
        {
            case 1:
            {
                cout << "Completeness: "<< completeness << "%" << endl;
                cout << "  =>  h = " << h << endl;
                cout << "  =>  Cost time (D:H:M:S): "<< days02 <<":"<< hours02 <<":"<< minutes02 <<":"<< seconds02 << endl;
                cout << "  =>  time per main loop   :  "<< programTimes07 << endl;
                cout << "  =>  time per delta h  :  "<< programTimes08 << endl;
                cout << "  =>  Cell number  :  "<< C1 << endl;
                cout << "  =>  threads  :  "<< nthreads << endl;
                cout << "**********************************************************"<< endl;
                break;
            }
            default:
            {
                switch (allpng)
                {
                    case 1:
                    {
                        cout << "Completeness: "<< completeness << "%" << endl;
                        cout << "  =>  h = " << h << endl;
                        cout << "  =>  Cost time (D:H:M:S): "<< days02 <<":"<< hours02 <<":"<< minutes02 <<":"<< seconds02 << endl;
                        cout << "  =>  time save CellTraceArray   :  "<< programTimes05 << endl;
                        cout << "  =>  time save PNG   :  "<< programTimes09 << endl;
                        cout << "  =>  time save CellArray   :  "<< programTimes10 << endl;
                        cout << "  =>  time save allPNG   :  "<< programTimes11 << endl;
                        cout << "  =>  time Dead Cell eliminateion   :  "<< programTimes12 << endl;
                        cout << "  =>  time Sort cell type   :  "<< programTimes14 << endl;
                        cout << "  =>  time Sort Random Number   :  "<< programTimes15 << endl;
                        cout << "  =>  time per Migration loop   :  "<< programTimes04 << endl;
                        cout << "  =>  time Sort Cell division   :  "<< programTimes13 << endl;
                        cout << "  =>  time costs in Cell Division   :  "<< programTimes16 << endl;
                        cout << "  =>  time per Division loop   :  "<< programTimes06 << endl;
                        cout << "  =>  time per delta h  :  "<< programTimes00 << endl;
                        cout << "  =>  Cell number  :  "<< C1 << endl;
                        cout << "  =>  threads  :  "<< nthreads << endl;
                        cout << "**********************************************************"<< endl;
                        break;
                    }
                    default:
                    {
                        cout << "Completeness: "<< completeness << "%" << endl;
                        cout << "  =>  h = " << h << endl;
                        cout << "  =>  Cost time (D:H:M:S): "<< days02 <<":"<< hours02 <<":"<< minutes02 <<":"<< seconds02 << endl;
                        cout << "  =>  time save CellTraceArray   :  "<< programTimes05 << endl;
                        cout << "  =>  time save PNG   :  "<< programTimes09 << endl;
                        cout << "  =>  time save CellArray   :  "<< programTimes10 << endl;
                        cout << "  =>  time Dead Cell eliminateion   :  "<< programTimes12 << endl;
                        cout << "  =>  time Sort cell type   :  "<< programTimes14 << endl;
                        cout << "  =>  time Sort Random Number   :  "<< programTimes15 << endl;
                        cout << "  =>  time per Migration loop   :  "<< programTimes04 << endl;
                        cout << "  =>  time Sort Cell division   :  "<< programTimes13 << endl;
                        cout << "  =>  time costs in Cell Division   :  "<< programTimes16 << endl;
                        cout << "  =>  time per Division loop   :  "<< programTimes06 << endl;
                        cout << "  =>  time per delta h  :  "<< programTimes00 << endl;
                        cout << "  =>  Cell number  :  "<< C1 << endl;
                        cout << "  =>  threads  :  "<< nthreads << endl;
                        cout << "**********************************************************"<< endl;
                        break;
                    }
                }
                break;
            }
        }
        h=h+deltah;
    }
}

#endif /* free_living_growth_recovery_hpp */
