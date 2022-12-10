//
//  SaveCellArraySingleCell.hpp
//  CCDS
//
//  Created by Tao Lee on 12/11/22.
//  Copyright © 2022 Tao Lee. All rights reserved.
//

#ifndef SaveCellArraySingleCell_hpp
#define SaveCellArraySingleCell_hpp

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <pngwriter.h>
#include <cmath>

using namespace blitz;
void SaveCellArraySingleCell(int &T, float alpha, float beta, Array<float,2> cell_array ,int Col)
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
                fprintf(fid3,"%g\t",cell_array(i,co));
            }
            else
            {
                fprintf(fid3,"%g\n",cell_array(i,co));
            }
        }
    }
    fclose(fid3);
}
#endif /* SaveCellArraySingleCell_hpp */
