//
//  outer_initiation_array_low_density.hpp
//  CCDS
//
//  Created by Taolee on 2019/4/13.
//  Copyright © 2019 Tao Lee. All rights reserved.
//

#ifndef outer_initiation_array_low_density_hpp
#define outer_initiation_array_low_density_hpp

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

inline CellStore outer_initiation_low_density_cell_store(int N0, int Visual_range_x, int Visual_range_y, const Array<int,2> &A, double uniup_r, double unilow_r, double sigmahatr,double muhatr, double uniup_K, double unilow_K, double sigmahatK,double muhatK, int N0r,int N0K, double *migration_rate_r, double *migration_rate_K)
{
    const long rng_context = 10002;
    double initial_r_growth_rate[N0r];
    double initial_K_growth_rate[N0K];
    Array<double,2> cor(2,N0+1,FortranArray<2>());
    cor=0;

    int number_cor=1;
    for (int x=1; x<=Visual_range_x/2; x++)
    {
        for (int y=1; y<=Visual_range_y/2; y++)
        {
            if (A(x,y)==1)
            {
                cor(1,number_cor)=x;
                cor(2,number_cor)=y;
                number_cor=number_cor+1;
            }
        }
    }
    double random_cor[N0];
    for (int x=0; x<N0; x++)
    {
        random_cor[x]=x+1;
    }
    stateless_shuffle(random_cor, random_cor+N0, rng_context, N0, 1);
    Array<double,2> cell_array_cor(2,N0,FortranArray<2>());
    for (int x=1; x<=N0; x++)
    {
        int seed = random_cor[x-1];
        cell_array_cor(1,x)=cor(1,seed);
        cell_array_cor(2,x)=cor(2,seed);
    }
    double rangr2 = uniup_r - unilow_r;
    for (int x=1; x<=N0r; x++)
    {
        double rand1 = (stateless_uniform(rng_context, 0, 1000 + x)*rangr2)+unilow_r;
        initial_r_growth_rate[x-1]=gsl_cdf_gaussian_Pinv(rand1, sigmahatr) + muhatr;
    }

    double rangK2 = uniup_K - unilow_K;
    for (int x=1; x<=N0K; x++)
    {
        double rand2 = (stateless_uniform(rng_context, 0, 2000 + x)*rangK2)+unilow_K;
        initial_K_growth_rate[x-1]=gsl_cdf_gaussian_Pinv(rand2, sigmahatK) + muhatK;
    }

    CellStore cells(cell_col::kStandardColumnCount);
    cells.reserve(N0);
    int a=0;
    int Nr=N0r*0.6;
    int Nk=N0*0.8;
    for (int x=1; x<=N0; x++)
    {
        int type=0;
        double growth_rate=0;
        double migration_rate=0;
        if (x<=Nr)
        {
            type=1;
            growth_rate=initial_r_growth_rate[x-1];
            migration_rate=migration_rate_r[x-1];
        }
        else if (x>N0r && x<=Nk)
        {
            type=2;
            growth_rate=initial_K_growth_rate[a];
            migration_rate=migration_rate_K[a];
            a++;
        }
        else
        {
            continue;
        }

        int x1 = (2*cell_array_cor(1,x))-1;
        int x2 = x1;
        int x3 = x1+1;
        int x4 = x1+1;
        int y1 = (2*cell_array_cor(2,x))-1;
        int y2 = y1+1;
        int y3 = y1+1;
        int y4 = y1;
        if (x1<1 || y1<1)
        {
            continue;
        }

        cells.push_empty();
        int target_row=cells.rows();
        cells(target_row,1)=x1;
        cells(target_row,2)=x2;
        cells(target_row,3)=x3;
        cells(target_row,4)=x4;
        cells(target_row,5)=y1;
        cells(target_row,6)=y2;
        cells(target_row,7)=y3;
        cells(target_row,8)=y4;
        cells(target_row,9)=type;
        cells(target_row,10)=growth_rate;
        cells(target_row,11)=growth_rate;
        cells(target_row,12)=migration_rate;
        cells(target_row,14)=0;
        cells(target_row,15)=target_row;
        cells(target_row,22)=1;
    }
    return cells;
}

#endif /* outer_initiation_array_low_density_hpp */
