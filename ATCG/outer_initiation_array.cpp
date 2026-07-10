//
//  outer_initiation_array.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#include "outer_initiation_array.hpp"

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
#include "int_grid.hpp"
using namespace blitz;

void fill_outer_initiation_cells(CellStore &cell_array_out_1, int N0, int Visual_range_x, int Visual_range_y, const IntGrid &A, double uniup_r, double unilow_r, double sigmahatr,double muhatr, double uniup_K, double unilow_K, double sigmahatK,double muhatK, int N0r,int N0K, double *migration_rate_r, double *migration_rate_K,int Col)
{
    const long rng_context = 10001;
    double initial_r_growth_rate[N0r];
    double initial_K_growth_rate[N0K];
    std::vector<int> cor_x(N0 + 1, 0);
    std::vector<int> cor_y(N0 + 1, 0);
    
    int number_cor=1;
    for (int x=1; x<=Visual_range_x/2; x++)
    {
        for (int y=1; y<=Visual_range_y/2; y++)
        {
            if (A(x,y)==1)
            {
                cor_x[number_cor]=x;
                cor_y[number_cor]=y;
                number_cor=number_cor+1;
            }
        }
    }
    std::vector<int> random_cor(N0);
    for (int x=0; x<N0; x++)
    {
        random_cor[x]=x+1;
    }
    stateless_shuffle(random_cor.data(), random_cor.data() + random_cor.size(), rng_context, N0, 1);
    std::vector<int> cell_cor_x(N0 + 1, 0);
    std::vector<int> cell_cor_y(N0 + 1, 0);
    for (int x=1; x<=N0; x++)
    {
        int seed = random_cor[x-1];
        cell_cor_x[x]=cor_x[seed];
        cell_cor_y[x]=cor_y[seed];
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
    int a=0;
    for (int x=1; x<=N0; x++)
    {
        if (x<=N0r)
        {
            int x1 = (2*cell_cor_x[x])-1;
            int x2 = x1;
            int x3 = x1+1;
            int x4 = x1+1;
            int y1 = (2*cell_cor_y[x])-1;
            int y2 = y1+1;
            int y3 = y1+1;
            int y4 = y1;
            int row = x - 1;
            cell_array_out_1.x1()[row]=x1;
            cell_array_out_1.x2()[row]=x2;
            cell_array_out_1.x3()[row]=x3;
            cell_array_out_1.x4()[row]=x4;
            cell_array_out_1.y1()[row]=y1;
            cell_array_out_1.y2()[row]=y2;
            cell_array_out_1.y3()[row]=y3;
            cell_array_out_1.y4()[row]=y4;
            cell_array_out_1.type()[row]=1;
            cell_array_out_1.growth_rate()[row]=initial_r_growth_rate[x-1];
            cell_array_out_1.density_growth_rate()[row]=initial_r_growth_rate[x-1];
            cell_array_out_1.migration_rate_base()[row]=migration_rate_r[x-1];
            cell_array_out_1.random_label()[row]=stateless_uniform(rng_context, 0, 3000 + x);
            cell_array_out_1.stage()[row]=0;
            cell_array_out_1.id()[row]=x;
            cell_array_out_1.viability()[row]=1;
            if (Col>=cell_col::kCellTraceLabel)
            {
                cell_array_out_1.cell_trace_label()[row]=x;
            }
        }
        else
        {
            int x1 = (2*cell_cor_x[x])-1;
            int x2 = x1;
            int x3 = x1+1;
            int x4 = x1+1;
            int y1 = (2*cell_cor_y[x])-1;
            int y2 = y1+1;
            int y3 = y1+1;
            int y4 = y1;
            int row = x - 1;
            cell_array_out_1.x1()[row]=x1;
            cell_array_out_1.x2()[row]=x2;
            cell_array_out_1.x3()[row]=x3;
            cell_array_out_1.x4()[row]=x4;
            cell_array_out_1.y1()[row]=y1;
            cell_array_out_1.y2()[row]=y2;
            cell_array_out_1.y3()[row]=y3;
            cell_array_out_1.y4()[row]=y4;
            cell_array_out_1.type()[row]=2;
            cell_array_out_1.growth_rate()[row]=initial_K_growth_rate[a];
            cell_array_out_1.density_growth_rate()[row]=initial_K_growth_rate[a];
            cell_array_out_1.migration_rate_base()[row]=migration_rate_K[a];
            cell_array_out_1.random_label()[row]=stateless_uniform(rng_context, 0, 3000 + x);
            cell_array_out_1.stage()[row]=0;
            cell_array_out_1.id()[row]=x;
            cell_array_out_1.viability()[row]=1;
            if (Col>=cell_col::kCellTraceLabel)
            {
                cell_array_out_1.cell_trace_label()[row]=x;
            }
            a++;
        }
    }
}

CellStore outer_initiation_cell_store(int N0, int Visual_range_x, int Visual_range_y, const IntGrid &A, double uniup_r, double unilow_r, double sigmahatr,double muhatr, double uniup_K, double unilow_K, double sigmahatK,double muhatK, int N0r,int N0K, double *migration_rate_r, double *migration_rate_K,int Col)
{
    CellStore cells(Col);
    cells.resize(N0);
    fill_outer_initiation_cells(cells, N0, Visual_range_x, Visual_range_y, A, uniup_r, unilow_r, sigmahatr, muhatr, uniup_K, unilow_K, sigmahatK, muhatK, N0r, N0K, migration_rate_r, migration_rate_K, Col);
    return cells;
}
