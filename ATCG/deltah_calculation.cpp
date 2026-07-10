//
//  deltah_calculation.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#include "deltah_calculation.hpp"

#include <cmath>
#include <cstdlib>

namespace
{
int cmp( const void *a , const void *b )
{
    return *(double *)a > *(double *)b ? 1 : -1;
}
}

double deltah_calculation(int N0,double *migration_rate_r,int N0r, int &MMR, int DDM)
{
    qsort(migration_rate_r,N0r,sizeof(migration_rate_r[0]),cmp);
    double max_mig_r=migration_rate_r[N0r-1];
    double aT=1;
    MMR=(int)ceil(max_mig_r);
    double h;
    h = aT/(double)MMR; 
    if (DDM==0)
    {
        if(h>0.2)
        {
            MMR=5;
            h = aT/(double)MMR;
        }
    }
    return h;
}
