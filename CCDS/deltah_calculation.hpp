//
//  deltah_calculation.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef deltah_calculation_hpp
#define deltah_calculation_hpp

#include <stdio.h>
struct element
{
    double data;
    size_t index;
};
int compare(const void *a, const void *b)
{
    return (*(const element*)a).data - (*(const element*)b).data;
}

int cmp( const void *a , const void *b )
{
    return *(double *)a > *(double *)b ? 1 : -1;
}

template<typename T>
int count(T& x)
{
    int s1 = sizeof(x);
    int s2 = sizeof(x[0]);
    int result = s1 / s2;
    return result;
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
#endif /* deltah_calculation_hpp */
