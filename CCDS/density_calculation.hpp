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
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "deltah_calculation.hpp"

double density_calculation(int i, Array<int, 3> sub_visual, Array<int,3> Visual_range,Array<float, 2> cell_array)
{
    Range all = Range::all();
    int ar=70;
    int xar=(ar/2)-1;
    int yar=ar/2;
    int cell_small=ar*ar;
    int cell_big=cell_small*0.25;
    sub_visual.resize(ar,ar,4);
    sub_visual=0;
    sub_visual(all,all,all)=Visual_range(Range(cell_array(i,1)-xar,cell_array(i,1)+yar),Range(cell_array(i,5)-xar,cell_array(i,5)+yar),all);
    int cell_count[6000];
    int cc=0;
    for (int cx=0; cx<ar; cx++)
    {
        for(int cy=0; cy<ar; cy++)
        {
            cell_count[cc]=sub_visual(cx+1,cy+1,4);
            cc++;
        }
    }
    vector<int> mycellcount (cell_count, cell_count+cell_small);
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
            if (sub_visual(cx+1,cy+1,3)==2)
            {
            cell_count_small=cell_count_small+1;
            }
        }
    }
    long cell_number_final = cells_number + cell_count_small;
    double density;
    if (cell_array(i,14)==0)
    {
        density=cell_number_final/(double)cell_big;
    }
    else
    {
        density=cell_number_final/(double)cell_small;
    }
    return density;
}
#endif /* density_calculation_hpp */
