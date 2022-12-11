//
//  AliveCellCount.hpp
//  CCDS
//
//  Created by Tao Lee on 12/11/22.
//  Copyright © 2022 Tao Lee. All rights reserved.
//

#ifndef AliveCellCount_hpp
#define AliveCellCount_hpp

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
using namespace blitz;

int AliveCellCount(Array<float, 2> cell_array,int nthreads)
{
    int current_size=cell_array.rows();
    int sum =0;
    omp_set_num_threads(nthreads);
    #pragma omp parallel for schedule(dynamic) reduction(+:sum)
    {
        for (int CN=1; CN<=current_size; CN++)
        {
            if (cell_array(CN,22)==1)
            {
                sum++;
            }
        }
    }
    return sum;
}
#endif /* AliveCellCount_hpp */
