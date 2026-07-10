//
//  SaveVisualArrayL3.hpp
//  ATCG
//
//  Created by Tao Lee on 5/23/24.
//  Copyright © 2024 Tao Lee. All rights reserved.
//

#include "SaveVisualArrayL3.hpp"

#include <stdio.h>
#include "visual_range.hpp"
#include <pngwriter.h>
#include <cmath>
void SaveVisualArrayL3(int T, double alpha, double beta, const VisualRange &Visual_range, int Vx, int Vy)
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
                fprintf(fid8,"%ld\t",Visual_range.stage(i,co));
            }
            else
            {
                fprintf(fid8,"%ld\n",Visual_range.stage(i,co));
            }
        }
    }
    fclose(fid8);
}
