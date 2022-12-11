//
//  SaveCellTraceArray.hpp
//  CCDS
//
//  Created by Tao Lee on 12/11/22.
//  Copyright © 2022 Tao Lee. All rights reserved.
//

#ifndef SaveCellTraceArray_hpp
#define SaveCellTraceArray_hpp

#include <stdio.h>

#include <blitz/blitz.h>
#include <blitz/array.h>
#include <pngwriter.h>
#include <cmath>

using namespace blitz;
void SaveCellTraceArray(int &T, float alpha, float beta, Array<float,2> cell_trace)
{
    char filedir1 [100] = {'\0'};
    sprintf(filedir1, "./a_%.1f_b_%.1f_CellTrace/Cell_Trace_%.1d.txt",alpha,beta,T);
    FILE * fid8;
    fid8=fopen (filedir1,"w+");
    int C01 = cell_trace.rows();
    for (int i=1;i<=C01;i++)
    {
        for(int co=1;co<=1000;co++)
        {
            if(co<1000)
            {
                fprintf(fid8,"%g\t",cell_trace(i,co));
            }
            else
            {
                fprintf(fid8,"%g\n",cell_trace(i,co));
            }
        }
    }
    fclose(fid8);

    
    char dirname2 [100] = {'\0'};
    int tt=T-1;
    sprintf(dirname2, "rm ./a_%.1f_b_%.1f_CellTrace/Cell_Trace_%.1d.txt",alpha,beta,tt);
    system(dirname2);
}




#endif /* SaveCellTraceArray_hpp */
