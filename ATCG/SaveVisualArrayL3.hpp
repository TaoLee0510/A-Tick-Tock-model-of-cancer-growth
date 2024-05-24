//
//  SaveVisualArrayL3.hpp
//  ATCG
//
//  Created by Tao Lee on 5/23/24.
//  Copyright © 2024 Tao Lee. All rights reserved.
//

#ifndef SaveVisualArrayL3_hpp
#define SaveVisualArrayL3_hpp

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <pngwriter.h>
#include <cmath>
using namespace blitz;
void SaveVisualArrayL3(int T, double alpha, double beta, Array<long,3> Visual_range, int Vx, int Vy)
{
    char filedir1 [100] = {'\0'};
    sprintf(filedir1, "./a_%.1f_b_%.1f_Visual_range/Visual_range_layer_3_%.1d.txt",alpha,beta,T);
    FILE * fid8;
    fid8=fopen (filedir1,"w+");
    for (int i=1;i<=Vx;i++)
    {
        for(int co=1;co<=Vy;co++)
        {
            if(co<Vy)
            {
                fprintf(fid8,"%ld\t",Visual_range(i,co,3));
            }
            else
            {
                fprintf(fid8,"%ld\n",Visual_range(i,co,3));
            }
        }
    }
    fclose(fid8);
}

#endif /* SaveVisualArrayL3_hpp */
