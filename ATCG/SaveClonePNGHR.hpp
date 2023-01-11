//
//  SaveClonePNGHR.hpp
//  CCDS
//
//  Created by Tao Lee on 1/1/23.
//  Copyright © 2023 Tao Lee. All rights reserved.
//

#ifndef SaveClonePNGHR_hpp
#define SaveClonePNGHR_hpp


#include <blitz/blitz.h>
#include <blitz/array.h>
#include <pngwriter.h>
#include <cmath>

using namespace blitz;
void SaveClonePNGHR(int Visual_range_x, int Visual_range_y, Array<double,2> cell_array, int H, int &T, double alpha, double beta,double deltah, Array<double,2> colorspace)
{
    double TT=deltah*(double)3600*(double)H;
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    char filedir13 [100] = {'\0'};
    sprintf(filedir13, "./a_%.1f_b_%.1f_clonepicsall/%.1d.png",alpha,beta,H);
    char filedir14 [100] = {'\0'};
    sprintf(filedir14, "%.08d s",(int)TT);;
    char filedir15 [100] = {'\0'};
    sprintf(filedir15, "/Users/taolee/Library/Fonts/Calisto MT.ttf");
    FILE * fid7;
    fid7=fopen (filedir13,"wb");
    pngwriter image3(Visual_range_x, Visual_range_y, 0, filedir13);
    
    int C0 = cell_array.rows();
    for (int i=1;i<=C0;i++)
    {
        int x= cell_array(i,1);
        int y= cell_array(i,5);
        int cell_stage=cell_array(i,14);
        int cell_index=cell_array(i,15);
        
        if(cell_stage==0)
        {
            image3.plot(x, y, colorspace(cell_index,2), colorspace(cell_index,3), colorspace(cell_index,4));
            image3.plot(x+1, y+1, colorspace(cell_index,2), colorspace(cell_index,3), colorspace(cell_index,4));
            image3.plot(x+1, y, colorspace(cell_index,2), colorspace(cell_index,3), colorspace(cell_index,4));
            image3.plot(x, y+1, colorspace(cell_index,2), colorspace(cell_index,3), colorspace(cell_index,4));
        }
        else
        {
            image3.plot(x, y, colorspace(cell_index,2), colorspace(cell_index,3), colorspace(cell_index,4));
        }
        
    }
    int posx=Visual_range_x-(Visual_range_x * 0.15);
    int posy=Visual_range_y-(Visual_range_y * 0.1);
    image3.plot_text(filedir15, 30, posx, posy, 0.0, filedir14, 1.0, 1.0, 1.0);
    image3.close();
    fclose(fid7);
}


#endif /* SaveClonePNGHR_hpp */
