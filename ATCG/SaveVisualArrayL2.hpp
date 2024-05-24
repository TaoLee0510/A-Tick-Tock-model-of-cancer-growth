//
//  SaveVisualArrayL2.hpp
//  ATCG
//
//  Created by Tao Lee on 5/23/24.
//  Copyright © 2024 Tao Lee. All rights reserved.
//

#ifndef SaveVisualArrayL2_hpp
#define SaveVisualArrayL2_hpp

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <pngwriter.h>
#include <cmath>
using namespace blitz;
void SaveVisualArrayL2(int T, double alpha, double beta, Array<long,3> Visual_range, int Vx, int Vy)
{
    char filedir1 [100] = {'\0'};
    sprintf(filedir1, "./a_%.1f_b_%.1f_Visual_range/Visual_range_layer_2_%.1d.txt",alpha,beta,T);
    FILE * fid8;
    fid8=fopen (filedir1,"w+");
    for (int i=1;i<=Vx;i++)
    {
        for(int co=1;co<=Vy;co++)
        {
            if(co<Vy)
            {
                fprintf(fid8,"%ld\t",Visual_range(i,co,2));
            }
            else
            {
                fprintf(fid8,"%ld\n",Visual_range(i,co,2));
            }
        }
    }
    fclose(fid8);
}

#endif /* SaveVisualArrayL2_hpp */
