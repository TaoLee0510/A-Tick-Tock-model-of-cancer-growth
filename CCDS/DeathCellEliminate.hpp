//
//  DeathCellEliminate.hpp
//  CCDS
//
//  Created by Tao Lee on 12/11/22.
//  Copyright © 2022 Tao Lee. All rights reserved.
//

#ifndef DeathCellEliminate_hpp
#define DeathCellEliminate_hpp

#include <stdio.h>
#include <omp.h>
#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#define BZ_THREADSAFE
#define BZ_THREADSAFE_USE_OPENMP
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "sortRow.hpp"
#include "sortRowReverse.hpp"


using namespace blitz;

void DeathCellEliminate(Array<int,3> &Visual_range, Array<float, 2> &cell_array, Array<float,2> &cell_array_temp,Array<float, 2> cell_array1,int alive_cell_number, int Col, int nthreads)
{
    Range all=Range::all();
    int current_size=cell_array.rows();
    
//    cell_array_temp.resize(alive_cell_number,Col);
//    cell_array_temp=0;
//
    sortRowReverse(cell_array, cell_array1,Col,22,nthreads);///sort cell type
    ///

    for (int site=current_size; site>alive_cell_number; site--)
    {
        if(cell_array(site,14)==0)
        {
            Visual_range(cell_array(site,1),cell_array(site,5),all)=0;
            Visual_range(cell_array(site,2),cell_array(site,6),all)=0;
            Visual_range(cell_array(site,3),cell_array(site,7),all)=0;
            Visual_range(cell_array(site,4),cell_array(site,8),all)=0;
        }
        else if(cell_array(site,14)==1)
        {
            Visual_range(cell_array(site,1),cell_array(site,5),all)=0;
        }
    }
    
    cell_array.resizeAndPreserve(alive_cell_number,Col);
    
}



#endif /* DeathCellEliminate_hpp */
