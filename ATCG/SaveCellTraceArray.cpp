//
//  SaveCellTraceArray.hpp
//  CCDS
//
//  Created by Tao Lee on 12/11/22.
//  Copyright © 2022 Tao Lee. All rights reserved.
//

#include "SaveCellTraceArray.hpp"

#include <stdio.h>

#include <pngwriter.h>
#include <cmath>
#include "cell_trace.hpp"

void SaveCellTraceArray(int T, double alpha, double beta, const CellTraceStore &cell_trace)
{
    char filedir1 [100] = {'\0'};
    snprintf(filedir1, sizeof(filedir1), "./a_%.1f_b_%.1f_CellTrace/Cell_Trace_%.1d.txt",alpha,beta,T);
    FILE * fid8;
    fid8=fopen (filedir1,"w+");
    int C01 = cell_trace.rows();
    for (int i=1;i<=C01;i++)
    {
        for(int co=1;co<=150;co++)
        {
            if(co<150)
            {
                fprintf(fid8,"%ld\t",cell_trace(i,co));
            }
            else
            {
                fprintf(fid8,"%ld\n",cell_trace(i,co));
            }
        }
    }
    fclose(fid8);


    if (T > 0)
    {
        char filedir_prev [100] = {'\0'};
        int tt=T-1;
        snprintf(filedir_prev, sizeof(filedir_prev), "./a_%.1f_b_%.1f_CellTrace/Cell_Trace_%.1d.txt",alpha,beta,tt);
        remove(filedir_prev);
    }
}
