//
//  DeletPreviousCellTraceArray.hpp
//  CCDS
//
//  Created by Tao Lee on 12/11/22.
//  Copyright © 2022 Tao Lee. All rights reserved.
//

#ifndef DeletPreviousCellTraceArray_hpp
#define DeletPreviousCellTraceArray_hpp

#include <stdio.h>

#include <blitz/blitz.h>
#include <blitz/array.h>
#include <pngwriter.h>
#include <cmath>

using namespace blitz;
void DeletPreviousCellTraceArray(int &T, float alpha, float beta)
{
    char dirname2 [100] = {'\0'};
    int tt=T-1;
    sprintf(dirname2, "rm ./a_%.1f_b_%.1f_CellTrace/Cell_Trace_%.1d.txt",alpha,beta,tt);
    system(dirname2);
}

#endif /* DeletPreviousCellTraceArray_hpp */
