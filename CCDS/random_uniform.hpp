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
#include <ctime>
#include <blitz/blitz.h>
#include <blitz/array.h>
using namespace blitz;
Array<double,2> random_uniform (int N0)
{
    const gsl_rng_type *T1;
    gsl_rng *r1;
    gsl_rng_env_setup();
    T1 = gsl_rng_ranlxs0;
    gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    r1 = gsl_rng_alloc(T1);
    Array<double,2> A(1,N0,FortranArray<2>());
    for(int x=1;x<=N0;x++)
    {
        A(1,x) = gsl_rng_uniform(r1);
    }
    return A;
    gsl_rng_free(r1);
}
#endif /* random_uniform_hpp */
