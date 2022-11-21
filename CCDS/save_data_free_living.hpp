//
//  save_data_free_living.hpp
//  CCDS
//
//  Created by Tao Lee on 11/10/22.
//  Copyright © 2022 Tao Lee. All rights reserved.
//

#ifndef save_data_free_living_hpp
#define save_data_free_living_hpp

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <pngwriter.h>
#include <cmath>

using namespace blitz;
void save_data_free_living(int Visual_range_x, int Visual_range_y, int N0, int N00, int N01, int MMR, int H, int &T, float alpha, float beta, Array<float,2> cell_array, int migration_judgement,double deltah, Array<double,2> colorspace, int DDM, int allpng,int Col)
{
    if (H%MMR==0)
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
        /////////////////////////////////////////////////PNG//////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        char filedir3 [100] = {'\0'};
        sprintf(filedir3, "./a_%.1f_b_%.1f/Cell_array_a_%.1f_b_%.1f_h_%.1d.txt",alpha,beta,alpha,beta,T);
        FILE * fid3;
        fid3=fopen (filedir3,"w+");
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
        int posx=Visual_range_x-(Visual_range_x * 0.15);
        int posy=Visual_range_y-(Visual_range_y * 0.1);
        image.plot_text(filedir6, 30, posx, posy, 0.0, filedir5, 1.0, 1.0, 1.0);
        image.close();
        fclose(fid4);
        fclose(fid3);
        T++;
    }
    if (allpng==1)
    {
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////PNG//////////////////////////////////////////////////////////////
        double TT=deltah*(double)3600*(double)H;
        char filedir7 [100] = {'\0'};
        sprintf(filedir7, "./a_%.1f_b_%.1f_picsall/%.1d.png",alpha,beta,H);
        char filedir8 [100] = {'\0'};
        sprintf(filedir8, "%.08d s",(int)TT);
        char filedir9 [100] = {'\0'};
        sprintf(filedir9, "/Users/taolee/Library/Fonts/Calisto MT.ttf");
        FILE * fid5;
        fid5=fopen (filedir7,"wb");
        pngwriter image2(Visual_range_x, Visual_range_y, 0, filedir7);
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        char filedir16 [100] = {'\0'};
        sprintf(filedir16, "./a_%.1f_b_%.1f_all/Cell_array_%.1d.txt",alpha,beta,H);
        FILE * fid8;
        fid8=fopen (filedir16,"w+");
        
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
            for(int co=1;co<=Col;co++)
            {
                if(co<Col)
                {
                    fprintf(fid8,"%g\t",cell_array(i,co));
                }
                else
                {
                    fprintf(fid8,"%g\n",cell_array(i,co));
                }
            }
        }
        int posx=Visual_range_x-(Visual_range_x * 0.15);
        int posy=Visual_range_y-(Visual_range_y * 0.1);
        image2.plot_text(filedir9, 30, posx, posy, 0.0, filedir8, 1.0, 1.0, 1.0);
        image2.close();
        fclose(fid5);
        fclose(fid8);
    }
}

#endif /* save_data_free_living_hpp */
