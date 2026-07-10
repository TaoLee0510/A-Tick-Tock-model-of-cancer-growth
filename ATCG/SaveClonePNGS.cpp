//
//  SaveClonePNGS.hpp
//  CCDS
//
//  Created by Tao Lee on 1/1/23.
//  Copyright © 2023 Tao Lee. All rights reserved.
//

#include "SaveClonePNGS.hpp"

#include <stdio.h>
#include <pngwriter.h>
#include <cmath>
#include "color_space.hpp"

void SaveClonePNGS(int Visual_range_x, int Visual_range_y, int &T, double alpha, double beta, const CellStore &cell_array, const ColorSpace &colorspace)
{
    char filedir10 [100] = {'\0'};
    sprintf(filedir10, "./a_%.1f_b_%.1f_clonepics/%.1d.png",alpha,beta,T);
    char filedir11 [100] = {'\0'};
    sprintf(filedir11, "%.04d h",T);
    char filedir12 [100] = {'\0'};
    sprintf(filedir12, "/Users/taolee/Library/Fonts/Calisto MT.ttf");
    FILE * fid6;
    fid6=fopen (filedir10,"wb");
    pngwriter image1(Visual_range_x, Visual_range_y, 0, filedir10);
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int C0 = cell_array.rows();
    for (int i=1;i<=C0;i++)
    {
        int x= cell_array.x1()[i - 1];
        int y= cell_array.y1()[i - 1];
//        int cell_type=cell_array.type()[i - 1];
        int cell_stage=cell_array.stage()[i - 1];
        int cell_index=cell_array.id()[i - 1];

        if(cell_stage==0)
        {
            image1.plot(x, y, colorspace(cell_index,2), colorspace(cell_index,3), colorspace(cell_index,4));
            image1.plot(x+1, y+1, colorspace(cell_index,2), colorspace(cell_index,3), colorspace(cell_index,4));
            image1.plot(x+1, y, colorspace(cell_index,2), colorspace(cell_index,3), colorspace(cell_index,4));
            image1.plot(x, y+1, colorspace(cell_index,2), colorspace(cell_index,3), colorspace(cell_index,4));
        }
        else
        {
            image1.plot(x, y, colorspace(cell_index,2), colorspace(cell_index,3), colorspace(cell_index,4));
        }
    }
    int posx=Visual_range_x-(Visual_range_x * 0.15);
    int posy=Visual_range_y-(Visual_range_y * 0.1);
    image1.plot_text(filedir12, 30, posx, posy, 0.0, filedir11, 1.0, 1.0, 1.0);
    image1.close();
    fclose(fid6);
}
