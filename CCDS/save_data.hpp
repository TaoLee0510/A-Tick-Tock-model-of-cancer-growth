//
//  save_data.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef save_data_hpp
#define save_data_hpp

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <pngwriter.h>
#include <cmath>

using namespace blitz;
void save_data(int Visual_range_x, int Visual_range_y, int N0, int N00, int N01, int MMR, int H, int &T, float alpha, float beta, Array<float,2> cell_array, int migration_judgement,double deltah, Array<double,2> colorspace, int DDM, int allpng,int free_living)
{
    if (H%MMR==0)
    {
        if (DDM==1 && T>=1  && migration_judgement==0 && free_living==0)
        {
            cout << "error: No cells migrated during the first hour, simulation aborted" <<endl;
            exit(0);
        }
        //int TT=int(deltah*(double)3600*(double)H);
        //if(TT%3600==0)
        //{
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
                int cell_index=cell_array(i,15);
                
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
                    image1.plot(x, y, colorspace(cell_index,2), colorspace(cell_index,3), colorspace(cell_index,4));
                    image1.plot(x+1, y+1, colorspace(cell_index,2), colorspace(cell_index,3), colorspace(cell_index,4));
                    image1.plot(x+1, y, colorspace(cell_index,2), colorspace(cell_index,3), colorspace(cell_index,4));
                    image1.plot(x, y+1, colorspace(cell_index,2), colorspace(cell_index,3), colorspace(cell_index,4));
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
                    image1.plot(x, y, colorspace(cell_index,2), colorspace(cell_index,3), colorspace(cell_index,4));
                }
                for(int co=1;co<=28;co++)
                {
                    if(co<28)
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
            image1.plot_text(filedir12, 30, posx, posy, 0.0, filedir11, 1.0, 1.0, 1.0);
            image1.close();
            fclose(fid6);
            fclose(fid3);
            T++;
        //}
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
        char filedir13 [100] = {'\0'};
        sprintf(filedir13, "./a_%.1f_b_%.1f_clonepicsall/%.1d.png",alpha,beta,H);
        char filedir14 [100] = {'\0'};
        sprintf(filedir14, "%.08d s",(int)TT);;
        char filedir15 [100] = {'\0'};
        sprintf(filedir15, "/Users/taolee/Library/Fonts/Calisto MT.ttf");
        FILE * fid7;
        fid7=fopen (filedir13,"wb");
        pngwriter image3(Visual_range_x, Visual_range_y, 0, filedir13);
        
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
            int cell_index=cell_array(i,15);
            
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
                image3.plot(x, y, colorspace(cell_index,2), colorspace(cell_index,3), colorspace(cell_index,4));
                image3.plot(x+1, y+1, colorspace(cell_index,2), colorspace(cell_index,3), colorspace(cell_index,4));
                image3.plot(x+1, y, colorspace(cell_index,2), colorspace(cell_index,3), colorspace(cell_index,4));
                image3.plot(x, y+1, colorspace(cell_index,2), colorspace(cell_index,3), colorspace(cell_index,4));
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
                image3.plot(x, y, colorspace(cell_index,2), colorspace(cell_index,3), colorspace(cell_index,4));
            }
            for(int co=1;co<=28;co++)
            {
                if(co<28)
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
        image3.plot_text(filedir15, 30, posx, posy, 0.0, filedir14, 1.0, 1.0, 1.0);
        image3.close();
        fclose(fid7);
        fclose(fid8);
    }
}
#endif /* save_data_hpp */
