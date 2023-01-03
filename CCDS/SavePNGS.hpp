//
//  SavePNGS.hpp
//  CCDS
//
//  Created by Tao Lee on 1/1/23.
//  Copyright © 2023 Tao Lee. All rights reserved.
//

#ifndef SavePNGS_hpp
#define SavePNGS_hpp


#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <pngwriter.h>
#include <cmath>

using namespace blitz;
void SavePNGS(int Visual_range_x, int Visual_range_y, int &T, double alpha, double beta, Array<double,2> cell_array)
{
    /////////////////////////////////////////////////PNG//////////////////////////////////////////////////////////////
    char filedir4 [100] = {'\0'};
    sprintf(filedir4, "./a_%.1f_b_%.1f_pics/%.1d.png",alpha,beta,T);
    char filedir5 [100] = {'\0'};
    sprintf(filedir5, "%.04d h",T);
    char filedir6 [100] = {'\0'};
    sprintf(filedir6, "/Users/taolee/Library/Fonts/Calisto MT.ttf");
    FILE * fid4;
    fid4=fopen (filedir4,"wb");
    pngwriter image(Visual_range_x, Visual_range_y, 0, filedir4);
    int C0 = cell_array.rows();
    for (int i=1;i<=C0;i++)
    {
        int x= cell_array(i,1);
        int y= cell_array(i,5);
        int cell_type=cell_array(i,9);
        int cell_stage=cell_array(i,14);
        
        if(cell_stage==0)
        {
            if(cell_type==1)
            {
                image.plot(x, y, 0.0, 1.0, 0.0);
                image.plot(x+1, y+1, 0.0, 1.0, 0.0);
                image.plot(x+1, y, 0.0, 1.0, 0.0);
                image.plot(x, y+1, 0.0, 1.0, 0.0);
            }
            else
            {
                image.plot(x, y, 1.0, 0.0, 0.0);
                image.plot(x+1, y+1, 1.0, 0.0, 0.0);
                image.plot(x+1, y, 1.0, 0.0, 0.0);
                image.plot(x, y+1, 1.0, 0.0, 0.0);
            }
        }
        else
        {
            if(cell_type==1)
            {
                image.plot(x, y, 0.0, 1.0, 0.0);
            }
            else
            {
                image.plot(x, y, 1.0, 0.0, 0.0);
            }
        }
    }
    int posx=Visual_range_x-(Visual_range_x * 0.15);
    int posy=Visual_range_y-(Visual_range_y * 0.1);
    image.plot_text(filedir6, 30, posx, posy, 0.0, filedir5, 1.0, 1.0, 1.0);
    image.close();
    fclose(fid4);
}

#endif /* SavePNGS_hpp */
