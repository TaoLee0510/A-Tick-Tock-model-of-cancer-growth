//
//  inner_initiation_array.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef inner_initiation_array_hpp
#define inner_initiation_array_hpp

#include <stdio.h>
#include <random>
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
#include "stateless_rng.hpp"
#include "cell_store.hpp"
using namespace blitz;

template <typename CellArray>
inline void fill_inner_initiation_cells(int N0,int N01,int R0,int Visual_range_x, int Visual_range_y, CellArray &cell_array_inner, const Array<long,3> &Visual_range, double uniup_r1, double unilow_r1, double sigmahatr,double muhatr, double uniup_K1, double unilow_K1, double sigmahatK,double muhatK, int N0r1,int N0K1, double *migration_rate_r1, double *migration_rate_K1, int Col)
{
    const long rng_context = 10003;
    double initial_r_growth_rate[N0r1];
    double initial_K_growth_rate[N0K1];
    Array<double,2> cor(2,N01+1,FortranArray<2>());
    cor=0;
    int number_cor=1;
    for (int x=1; x<=Visual_range_x; x++)
    {
        for (int y=1; y<=Visual_range_y; y++)
        {
            if(pow((x-Visual_range_x/2),2)+pow((y-Visual_range_y/2),2)<=pow(R0-5,2))
            {
                if (Visual_range(x,y,1)==0)
                {
                    cor(1,number_cor)=x;
                    cor(2,number_cor)=y;
                    number_cor=number_cor+1;
                }
            }
        }
    }
    double random_cor[N01];
    for (int x=0; x<N01; x++)
    {
        random_cor[x]=x+1;
    }
    stateless_shuffle(random_cor, random_cor+N01, rng_context, N01, 1);
    Array<double,2> cell_array_cor(2,N01,FortranArray<2>());
    for (int x=1; x<=N01; x++)
    {
        int seed = random_cor[x-1];
        cell_array_cor(1,x)=cor(1,seed);
        cell_array_cor(2,x)=cor(2,seed);
    }
    double rangr2 = uniup_r1 - unilow_r1;
    for (int x=1; x<=N0r1; x++)
    {
        double rand1 = (stateless_uniform(rng_context, 0, 1000 + x)*rangr2)+unilow_r1;
        initial_r_growth_rate[x-1]=gsl_cdf_gaussian_Pinv(rand1, sigmahatr) + muhatr;
    }
    
    double rangK2 = uniup_K1 - unilow_K1;
    for (int x=1; x<=N0K1; x++)
    {
        double rand2 = (stateless_uniform(rng_context, 0, 2000 + x)*rangK2)+unilow_K1;
        initial_K_growth_rate[x-1]=gsl_cdf_gaussian_Pinv(rand2, sigmahatK) + muhatK;
    }
    int a=0;
    for (int i=1;i<=N01;i++)
    {
        if (i<=N0r1)
        {
            int x1=cell_array_cor(1,i);
            int y1=cell_array_cor(2,i);
            cell_array_inner(i,1)=x1;
            cell_array_inner(i,5)=y1;
            cell_array_inner(i,9)=1;
            cell_array_inner(i,10)=initial_r_growth_rate[i-1];
            cell_array_inner(i,11)=initial_r_growth_rate[i-1];
            cell_array_inner(i,12)=migration_rate_r1[i-1];
            cell_array_inner(i,14)=1;
            cell_array_inner(i,15)=i+N0;
            cell_array_inner(i,22)=1;
            if (Col>=cell_col::kCellTraceLabel)
            {
                cell_array_inner(i,cell_col::kCellTraceLabel)=i+N0;
            }
        }
        else
        {
            int x1=cell_array_cor(1,i);
            int y1=cell_array_cor(2,i);
            cell_array_inner(i,1)=x1;
            cell_array_inner(i,5)=y1;
            cell_array_inner(i,9)=2;
            cell_array_inner(i,10)=initial_K_growth_rate[a];
            cell_array_inner(i,11)=initial_K_growth_rate[a];
            cell_array_inner(i,12)=migration_rate_K1[a];
            cell_array_inner(i,14)=1;
            cell_array_inner(i,15)=i+N0;
            cell_array_inner(i,22)=1;
            if (Col>=cell_col::kCellTraceLabel)
            {
                cell_array_inner(i,cell_col::kCellTraceLabel)=i+N0;
            }
            a=a+1;
        }
    }
}

inline CellStore inner_initiation_cell_store(int N0,int N01,int R0,int Visual_range_x, int Visual_range_y, const Array<long,3> &Visual_range, double uniup_r1, double unilow_r1, double sigmahatr,double muhatr, double uniup_K1, double unilow_K1, double sigmahatK,double muhatK, int N0r1,int N0K1, double *migration_rate_r1, double *migration_rate_K1, int Col)
{
    CellStore cells(Col);
    cells.resize(N01);
    fill_inner_initiation_cells(N0, N01, R0, Visual_range_x, Visual_range_y, cells, Visual_range, uniup_r1, unilow_r1, sigmahatr, muhatr, uniup_K1, unilow_K1, sigmahatK, muhatK, N0r1, N0K1, migration_rate_r1, migration_rate_K1, Col);
    return cells;
}

#endif /* inner_initiation_array_hpp */
