//
//  SaveCellArray.hpp
//  CCDS
//
//  Created by Tao Lee on 1/1/23.
//  Copyright © 2023 Tao Lee. All rights reserved.
//

#include "SaveCellArray.hpp"

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <pngwriter.h>
#include <cmath>

using namespace blitz;
void SaveCellArray(int &T, double alpha, double beta, const CellStore &cell_array ,int Col)
{
    if(Col>28)
    {
        char filedir3 [100] = {'\0'};
        sprintf(filedir3, "./a_%.1f_b_%.1f/Cell_array_a_%.1f_b_%.1f_h_%.1d.txt",alpha,beta,alpha,beta,T);
        FILE * fid3;
        fid3=fopen (filedir3,"w+");
        int C0 = cell_array.rows();
        for (int i=1;i<=C0;i++)
        {
            for(int co=1;co<=Col;co++)
            {
                if(co<29)
                {
                    fprintf(fid3,"%g\t",cell_array.column(co)[i - 1]);
                }
                else if (co>=29 & co<Col)
                {
                    fprintf(fid3,"%ld\t",(long)cell_array.column(co)[i - 1]);
                }
                else
                {
                    fprintf(fid3,"%ld\n",(long)cell_array.column(co)[i - 1]);
                }
            }
        }
        fclose(fid3);
    }
    else
    {
        char filedir3 [100] = {'\0'};
        sprintf(filedir3, "./a_%.1f_b_%.1f/Cell_array_a_%.1f_b_%.1f_h_%.1d.txt",alpha,beta,alpha,beta,T);
        FILE * fid3;
        fid3=fopen (filedir3,"w+");
        int C0 = cell_array.rows();
        for (int i=1;i<=C0;i++)
        {
            for(int co=1;co<=Col;co++)
            {
                if(co<Col)
                {
                    fprintf(fid3,"%g\t",cell_array.column(co)[i - 1]);
                }
                else
                {
                    fprintf(fid3,"%g\n",cell_array.column(co)[i - 1]);
                }
            }
        }
        fclose(fid3);
    }
}
