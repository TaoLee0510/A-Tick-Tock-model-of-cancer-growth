//
//  random_uniform.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef random_uniform_hpp
#define random_uniform_hpp

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "stateless_rng.hpp"
using namespace blitz;
Array<double,2> random_uniform (int N0, long rng_context = 1, long rng_time_step = 0, long rng_event_base = 0)
{
    Array<double,2> A(1,N0,FortranArray<2>());
    for(int x=1;x<=N0;x++)
    {
        A(1,x) = stateless_uniform(rng_context, rng_time_step, rng_event_base + x);
    }
    return A;
}
#endif /* random_uniform_hpp */
