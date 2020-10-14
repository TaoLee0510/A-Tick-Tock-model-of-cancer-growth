//
//  outer_initiation_array.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef outer_initiation_array_hpp
#define outer_initiation_array_hpp

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
#include "random_uniform.hpp"
using namespace blitz;

Array<float,2> outer_initiation_array(int N0, int Visual_range_x, int Visual_range_y, Array<int,2> A, double uniup_r, double unilow_r, double sigmahatr,double muhatr, double uniup_K, double unilow_K, double sigmahatK,double muhatK, int N0r,int N0K, double *migration_rate_r, double *migration_rate_K)
{
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 RNG(seed);
    
    double initial_r_growth_rate[N0r];
    double initial_K_growth_rate[N0K];
    Array<float,2> cor(2,N0+1,FortranArray<2>());
    cor=0;
    Array<float,2> cell_array_out_1(N0,28,FortranArray<2>());
    cell_array_out_1=0;
    
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
    float random_cor[N0];
    for (int x=0; x<N0; x++)
    {
        random_cor[x]=x+1;
    }
    shuffle(random_cor, random_cor+N0,RNG);
    Array<float,2> cell_array_cor(2,N0,FortranArray<2>());
    for (int x=1; x<=N0; x++)
    {
        int seed = random_cor[x-1];
        cell_array_cor(1,x)=cor(1,seed);
        cell_array_cor(2,x)=cor(2,seed);
    }
    Array<float,2> radom_number(1,N0,FortranArray<2>());
    radom_number=random_uniform(N0);
    double rangr2 = uniup_r - unilow_r;
    for (int x=1; x<=N0r; x++)
    {
        double rand1 = (radom_number(1,x)*rangr2)+unilow_r;
        initial_r_growth_rate[x-1]=gsl_cdf_gaussian_Pinv(rand1, sigmahatr) + muhatr;
    }
    radom_number=random_uniform(N0);
    
    double rangK2 = uniup_K - unilow_K;
    for (int x=1; x<=N0K; x++)
    {
        double rand2 = (radom_number(1,x)*rangK2)+unilow_K;
        initial_K_growth_rate[x-1]=gsl_cdf_gaussian_Pinv(rand2, sigmahatK) + muhatK;
    }
    int a=0;
    for (int x=1; x<=N0; x++)
    {
        if (x<=N0r)
        {
            int x1 = (2*cell_array_cor(1,x))-1;
            int x2 = x1;
            int x3 = x1+1;
            int x4 = x1+1;
            int y1 = (2*cell_array_cor(2,x))-1;
            int y2 = y1+1;
            int y3 = y1+1;
            int y4 = y1;
            cell_array_out_1(x,1)=x1;
            cell_array_out_1(x,2)=x2;
            cell_array_out_1(x,3)=x3;
            cell_array_out_1(x,4)=x4;
            cell_array_out_1(x,5)=y1;
            cell_array_out_1(x,6)=y2;
            cell_array_out_1(x,7)=y3;
            cell_array_out_1(x,8)=y4;
            cell_array_out_1(x,9)=1;
            cell_array_out_1(x,10)=initial_r_growth_rate[x-1];
            cell_array_out_1(x,11)=initial_r_growth_rate[x-1];
            cell_array_out_1(x,12)=migration_rate_r[x-1];
            cell_array_out_1(x,14)=0;
            cell_array_out_1(x,15)=x;
            cell_array_out_1(x,22)=1;
        }
        else
        {
            int x1 = (2*cell_array_cor(1,x))-1;
            int x2 = x1;
            int x3 = x1+1;
            int x4 = x1+1;
            int y1 = (2*cell_array_cor(2,x))-1;
            int y2 = y1+1;
            int y3 = y1+1;
            int y4 = y1;
            cell_array_out_1(x,1)=x1;
            cell_array_out_1(x,2)=x2;
            cell_array_out_1(x,3)=x3;
            cell_array_out_1(x,4)=x4;
            cell_array_out_1(x,5)=y1;
            cell_array_out_1(x,6)=y2;
            cell_array_out_1(x,7)=y3;
            cell_array_out_1(x,8)=y4;
            cell_array_out_1(x,9)=2;
            cell_array_out_1(x,10)=initial_K_growth_rate[a];
            cell_array_out_1(x,11)=initial_K_growth_rate[a];
            cell_array_out_1(x,12)=migration_rate_K[a];
            cell_array_out_1(x,14)=0;
            cell_array_out_1(x,15)=x;
            cell_array_out_1(x,22)=1;
            a++;
        }
    }
    return cell_array_out_1;
}
#endif /* outer_initiation_array_hpp */

