//
//  density_calculation.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef density_calculation_hpp
#define density_calculation_hpp

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
#include "cell_columns.hpp"
#include "cell_store.hpp"
#include "visual_range.hpp"
#include "deltah_calculation.hpp"

inline double density_calculation_from_position(int x1, int y1, int cell_stage, const VisualRange &Visual_range)
{
    int ar=70;
//    int ar=30;
    int xar=(ar/2)-1;
    int yar=ar/2;
    int cell_small=ar*ar;
    int cell_number_limit=ar*ar;
    int cell_big=cell_small*0.25;
    long cell_count[cell_number_limit];////********
    int cc=0;
    for (int cx=0; cx<ar; cx++)
    {
        for(int cy=0; cy<ar; cy++)
        {
            cell_count[cc]=Visual_range.cell_label(x1 - xar + cx, y1 - xar + cy);
            cc++;
        }
    }
    vector<int> mycellcount (cell_count, cell_count+cell_number_limit);
    sort(mycellcount.begin(),mycellcount.end());
    mycellcount.erase(unique(mycellcount.begin(), mycellcount.end()), mycellcount.end());
    long cells_number=0;
    cells_number = mycellcount.size();
    if (mycellcount[0]==0)
    {
        cells_number=cells_number-1;
    }
    long cell_count_small=0;
    for (int cx=0; cx<ar; cx++)
    {
        for(int cy=0; cy<ar; cy++)
        {
            if (Visual_range.stage(x1 - xar + cx, y1 - xar + cy)==2)
            {
            cell_count_small=cell_count_small+1;
            }
        }
    }
    long cell_number_final = cells_number + cell_count_small;
    
    double density;
//    double density=cell_number_final/(double)cell_number_limit;
    
    if (cell_stage==0)
    {
        density=cell_number_final/(double)cell_big;
    }
    else
    {
        density=cell_number_final/(double)cell_small;
    }
    return density;
}

inline double density_calculation(int i, const VisualRange &Visual_range, const CellStore &cells)
{
    int row = i - 1;
    int x1 = (int)cells.x1()[row];
    int y1 = (int)cells.y1()[row];
    int cell_stage = (int)cells.stage()[row];
    return density_calculation_from_position(x1, y1, cell_stage, Visual_range);
}
#endif /* density_calculation_hpp */
