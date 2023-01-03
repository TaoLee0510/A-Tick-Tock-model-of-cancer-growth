//
//  SaveAllPNG.hpp
//  CCDS
//
//  Created by Tao Lee on 12/11/22.
//  Copyright © 2022 Tao Lee. All rights reserved.
//

#ifndef SaveAllPNG_hpp
#define SaveAllPNG_hpp

#include <blitz/blitz.h>
#include <blitz/array.h>
#include <pngwriter.h>
#include <cmath>

using namespace blitz;
void SaveAllPNG(int Visual_range_x, int Visual_range_y, Array<double,2> cell_array, int H, int &T, double alpha, double beta,double deltah)
{
    int DELTA=10;
    if (H%DELTA==0)
    {
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////PNG//////////////////////////////////////////////////////////////
        double TT=deltah*(double)3600*(double)H;
        int HH=H/DELTA;
        char filedir7 [100] = {'\0'};
        sprintf(filedir7, "./a_%.1f_b_%.1f_picsall/%.1d.png",alpha,beta,HH);
        char filedir8 [100] = {'\0'};
        sprintf(filedir8, "%.08d s",(int)TT);
        char filedir9 [100] = {'\0'};
        sprintf(filedir9, "/Users/taolee/Library/Fonts/Calisto MT.ttf");
        FILE * fid5;
        fid5=fopen (filedir7,"wb");
        pngwriter image2(Visual_range_x, Visual_range_y, 0, filedir7);
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
                    image2.plot(x, y, 0.0, 1.0, 0.0);
                    image2.plot(x+1, y+1, 0.0, 1.0, 0.0);
                    image2.plot(x+1, y, 0.0, 1.0, 0.0);
                    image2.plot(x, y+1, 0.0, 1.0, 0.0);
                }
                else
                {
                    image2.plot(x, y, 1.0, 0.0, 0.0);
                    image2.plot(x+1, y+1, 1.0, 0.0, 0.0);
                    image2.plot(x+1, y, 1.0, 0.0, 0.0);
                    image2.plot(x, y+1, 1.0, 0.0, 0.0);
                }
            }
            else
            {
                if(cell_type==1)
                {
                    image2.plot(x, y, 0.0, 1.0, 0.0);
                }
                else
                {
                    image2.plot(x, y, 1.0, 0.0, 0.0);
                }
            }
            
        }
        int posx=Visual_range_x-(Visual_range_x * 0.15);
        int posy=Visual_range_y-(Visual_range_y * 0.1);
        image2.plot_text(filedir9, 30, posx, posy, 0.0, filedir8, 1.0, 1.0, 1.0);
        image2.close();
        fclose(fid5);
    }

}

#endif /* SaveAllPNG_hpp */
