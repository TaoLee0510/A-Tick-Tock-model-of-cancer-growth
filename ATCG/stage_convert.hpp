//
//  stage_convert.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef stage_convert_hpp
#define stage_convert_hpp

#include <stdio.h>
#include <random>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_matrix.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
using namespace blitz;
void stage_convert(int Visual_range_x, int Visual_range_y, Array<double,2> &cell_array, Array<long,3> &Visual_range, int &cell_label,int utralsmall)
{
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 RNG(seed);
    Range all = Range::all();
    int C00= cell_array.rows();
    for (int x=1; x<=C00; ++x)
    {
        int stage_cor[4]={0};
        int direction[8]={0};
        if (cell_array(x,1)>=100 && cell_array(x,5) >=100  && cell_array(x,1)<=Visual_range_x+100 && cell_array(x,5)<=Visual_range_y+100)
        {
            if(cell_array(x,14)==1)
            {
                int x1 = cell_array(x,1);
                int y1 = cell_array(x,5);
                long loci_cor[9]={0};
                int loci_direction=0;
                for (int xx=x1-1;xx<=x1+1;xx++)
                {
                    for (int yy=y1-1;yy<=y1+1;yy++)
                    {
                        loci_cor[loci_direction]=Visual_range(xx,yy,1);
                        loci_direction=loci_direction+1;
                    }
                }
                if (loci_cor[0]==0 && loci_cor[1]==0 && loci_cor[3]==0)
                {
                    stage_cor[0]=1;
                }
                if (loci_cor[1]==0 && loci_cor[2]==0 && loci_cor[5]==0)
                {
                    stage_cor[1]=2;
                }
                if (loci_cor[5]==0 && loci_cor[7]==0 && loci_cor[8]==0)
                {
                    stage_cor[2]=3;
                }
                if (loci_cor[3]==0 && loci_cor[6]==0 && loci_cor[7]==0)
                {
                    stage_cor[3]=4;
                }
                int size=0;
                for(int a=0;a<4;a++)
                {
                    if(stage_cor[a]!=0)
                    {
                        size=size+1;
                    }
                }
                if (size>0)
                {
                    int *stage_cor_1=new int[size];
                    int num=0;
                    for(int a=0;a<4;a++)
                    {
                        if(stage_cor[a]!=0)
                        {
                            stage_cor_1[num]=stage_cor[a];
                            num=num+1;
                        }
                    }
                    shuffle(stage_cor_1,stage_cor_1+num,RNG);
                    int scor=stage_cor_1[0];
                    

                    if (scor==1)
                    {
                        cell_array(x,14)=0;
                        cell_array(x,1)=x1-1;
                        cell_array(x,2)=x1-1;
                        cell_array(x,3)=x1;
                        cell_array(x,4)=x1;
                        cell_array(x,5)=y1-1;
                        cell_array(x,6)=y1;
                        cell_array(x,7)=y1;
                        cell_array(x,8)=y1-1;
                        Visual_range((int)cell_array(x,1),(int)cell_array(x,5),1)=1;
                        Visual_range((int)cell_array(x,2),(int)cell_array(x,6),1)=1;
                        Visual_range((int)cell_array(x,3),(int)cell_array(x,7),1)=1;
                        Visual_range((int)cell_array(x,4),(int)cell_array(x,8),1)=1;
                        Visual_range((int)cell_array(x,1),(int)cell_array(x,5),2)=cell_array(x,15);
                        Visual_range((int)cell_array(x,2),(int)cell_array(x,6),2)=cell_array(x,15);
                        Visual_range((int)cell_array(x,3),(int)cell_array(x,7),2)=cell_array(x,15);
                        Visual_range((int)cell_array(x,4),(int)cell_array(x,8),2)=cell_array(x,15);
                        Visual_range((int)cell_array(x,1),(int)cell_array(x,5),3)=0;
                        Visual_range((int)cell_array(x,2),(int)cell_array(x,6),3)=0;
                        Visual_range((int)cell_array(x,3),(int)cell_array(x,7),3)=0;
                        Visual_range((int)cell_array(x,4),(int)cell_array(x,8),3)=0;
                        Visual_range((int)cell_array(x,1),(int)cell_array(x,5),4)=Visual_range(x1,y1,4);
                        Visual_range((int)cell_array(x,2),(int)cell_array(x,6),4)=Visual_range(x1,y1,4);
                        Visual_range((int)cell_array(x,3),(int)cell_array(x,7),4)=Visual_range(x1,y1,4);
                        Visual_range((int)cell_array(x,4),(int)cell_array(x,8),4)=Visual_range(x1,y1,4);
                        
 
                    }
                    else if (scor==2)
                    {
                        cell_array(x,14)=0;
                        cell_array(x,1)=x1-1;
                        cell_array(x,2)=x1-1;
                        cell_array(x,3)=x1;
                        cell_array(x,4)=x1;
                        cell_array(x,5)=y1;
                        cell_array(x,6)=y1+1;
                        cell_array(x,7)=y1+1;
                        cell_array(x,8)=y1;
                        Visual_range((int)cell_array(x,1),(int)cell_array(x,5),1)=1;
                        Visual_range((int)cell_array(x,2),(int)cell_array(x,6),1)=1;
                        Visual_range((int)cell_array(x,3),(int)cell_array(x,7),1)=1;
                        Visual_range((int)cell_array(x,4),(int)cell_array(x,8),1)=1;
                        Visual_range((int)cell_array(x,1),(int)cell_array(x,5),2)=cell_array(x,15);
                        Visual_range((int)cell_array(x,2),(int)cell_array(x,6),2)=cell_array(x,15);
                        Visual_range((int)cell_array(x,3),(int)cell_array(x,7),2)=cell_array(x,15);
                        Visual_range((int)cell_array(x,4),(int)cell_array(x,8),2)=cell_array(x,15);
                        Visual_range((int)cell_array(x,1),(int)cell_array(x,5),3)=0;
                        Visual_range((int)cell_array(x,2),(int)cell_array(x,6),3)=0;
                        Visual_range((int)cell_array(x,3),(int)cell_array(x,7),3)=0;
                        Visual_range((int)cell_array(x,4),(int)cell_array(x,8),3)=0;
                        Visual_range((int)cell_array(x,1),(int)cell_array(x,5),4)=Visual_range(x1,y1,4);
                        Visual_range((int)cell_array(x,2),(int)cell_array(x,6),4)=Visual_range(x1,y1,4);
                        Visual_range((int)cell_array(x,3),(int)cell_array(x,7),4)=Visual_range(x1,y1,4);
                        Visual_range((int)cell_array(x,4),(int)cell_array(x,8),4)=Visual_range(x1,y1,4);
                             
                        
                    }
                    else if (scor==3)
                    {
                        cell_array(x,14)=0;
                        cell_array(x,1)=x1;
                        cell_array(x,2)=x1;
                        cell_array(x,3)=x1+1;
                        cell_array(x,4)=x1+1;
                        cell_array(x,5)=y1;
                        cell_array(x,6)=y1+1;
                        cell_array(x,7)=y1+1;
                        cell_array(x,8)=y1;
                        Visual_range((int)cell_array(x,1),(int)cell_array(x,5),1)=1;
                        Visual_range((int)cell_array(x,2),(int)cell_array(x,6),1)=1;
                        Visual_range((int)cell_array(x,3),(int)cell_array(x,7),1)=1;
                        Visual_range((int)cell_array(x,4),(int)cell_array(x,8),1)=1;
                        Visual_range((int)cell_array(x,1),(int)cell_array(x,5),2)=cell_array(x,15);
                        Visual_range((int)cell_array(x,2),(int)cell_array(x,6),2)=cell_array(x,15);
                        Visual_range((int)cell_array(x,3),(int)cell_array(x,7),2)=cell_array(x,15);
                        Visual_range((int)cell_array(x,4),(int)cell_array(x,8),2)=cell_array(x,15);
                        Visual_range((int)cell_array(x,1),(int)cell_array(x,5),3)=0;
                        Visual_range((int)cell_array(x,2),(int)cell_array(x,6),3)=0;
                        Visual_range((int)cell_array(x,3),(int)cell_array(x,7),3)=0;
                        Visual_range((int)cell_array(x,4),(int)cell_array(x,8),3)=0;
                        Visual_range((int)cell_array(x,1),(int)cell_array(x,5),4)=Visual_range(x1,y1,4);
                        Visual_range((int)cell_array(x,2),(int)cell_array(x,6),4)=Visual_range(x1,y1,4);
                        Visual_range((int)cell_array(x,3),(int)cell_array(x,7),4)=Visual_range(x1,y1,4);
                        Visual_range((int)cell_array(x,4),(int)cell_array(x,8),4)=Visual_range(x1,y1,4);
   
                    }
                    else if (scor==4)
                    {
                        cell_array(x,14)=0;
                        cell_array(x,1)=x1;
                        cell_array(x,2)=x1;
                        cell_array(x,3)=x1+1;
                        cell_array(x,4)=x1+1;
                        cell_array(x,5)=y1-1;
                        cell_array(x,6)=y1;
                        cell_array(x,7)=y1;
                        cell_array(x,8)=y1-1;
                        Visual_range((int)cell_array(x,1),(int)cell_array(x,5),1)=1;
                        Visual_range((int)cell_array(x,2),(int)cell_array(x,6),1)=1;
                        Visual_range((int)cell_array(x,3),(int)cell_array(x,7),1)=1;
                        Visual_range((int)cell_array(x,4),(int)cell_array(x,8),1)=1;
                        Visual_range((int)cell_array(x,1),(int)cell_array(x,5),2)=cell_array(x,15);
                        Visual_range((int)cell_array(x,2),(int)cell_array(x,6),2)=cell_array(x,15);
                        Visual_range((int)cell_array(x,3),(int)cell_array(x,7),2)=cell_array(x,15);
                        Visual_range((int)cell_array(x,4),(int)cell_array(x,8),2)=cell_array(x,15);
                        Visual_range((int)cell_array(x,1),(int)cell_array(x,5),3)=0;
                        Visual_range((int)cell_array(x,2),(int)cell_array(x,6),3)=0;
                        Visual_range((int)cell_array(x,3),(int)cell_array(x,7),3)=0;
                        Visual_range((int)cell_array(x,4),(int)cell_array(x,8),3)=0;
                        Visual_range((int)cell_array(x,1),(int)cell_array(x,5),4)=Visual_range(x1,y1,4);
                        Visual_range((int)cell_array(x,2),(int)cell_array(x,6),4)=Visual_range(x1,y1,4);
                        Visual_range((int)cell_array(x,3),(int)cell_array(x,7),4)=Visual_range(x1,y1,4);
                        Visual_range((int)cell_array(x,4),(int)cell_array(x,8),4)=Visual_range(x1,y1,4);
                           
                        
                    }
                    delete[] stage_cor_1;
                    stage_cor_1 = NULL;
                }
            }
            else if (cell_array(x,14)==2)
            {
                int x1 = cell_array(x,1);
                int y1 = cell_array(x,5);
                Array<double,2> cor_small(3,3,FortranArray<2>());
                cor_small=0;
                cor_small(all,all)=Visual_range(Range(x1-1,x1+1),Range(y1-1,y1+1),1);
                if (cor_small(1,1)==0)
                {
                    direction[0]=1;
                }
                if (cor_small(1,2)==0)
                {
                    direction[1]=2;
                }
                if (cor_small(1,3)==0)
                {
                    direction[2]=3;
                }
                if (cor_small(2,3)==0)
                {
                    direction[3]=4;
                }
                if (cor_small(3,3)==0)
                {
                    direction[4]=5;
                }
                if (cor_small(3,2)==0)
                {
                    direction[5]=6;
                }
                if (cor_small(3,1)==0)
                {
                    direction[6]=7;
                }
                if (cor_small(2,1)==0)
                {
                    direction[7]=8;
                }
                int mloci=0;
                for (int mlo=0; mlo<8; mlo++)
                {
                    if (direction[mlo]!=0)
                    {
                        mloci++;
                    }
                }
                if (mloci>0)
                {
                    int *direction1=new int[mloci];
                    int new_loci=0;
                    for (int loci=0; loci<8; loci++)
                    {
                        if (direction[loci]!=0)
                        {
                            direction1[new_loci]=direction[loci];
                            new_loci++;
                        }
                    }
                    int order=0;
                    long length_dir=new_loci;
                    if (length_dir==0)
                    {
                        order=direction1[0];
                    }
                    else if (length_dir>0)
                    {
                        shuffle(direction1, direction1+new_loci,RNG);
                        order=direction1[0];
                    }
                    if (order==1)
                    {
                        int cell_label_1=cell_label+1;
                        Visual_range(x1-1,y1-1,1)=1;
                        Visual_range(x1-1,y1-1,2)=cell_array(x,15);
                        Visual_range(x1-1,y1-1,3)=1;
                        Visual_range(x1-1,y1-1,4)=cell_label_1;
                        cell_array(x,1)=x1-1;
                        cell_array(x,5)=y1-1;
                        cell_array(x,14)=1;
                        Visual_range(x1,y1,3)=1;
                        int size1=cell_array.rows();
                        for (int cus=1;cus<=size1;cus++)
                        {
                            if(cell_array(cus,1)==x1 && cell_array(cus,5)==y1)
                            {
                                cell_array(cus,14)=1;
                            }
                        }
                    }
                    else if (order==2)
                    {
                        int cell_label_1=cell_label+1;
                        Visual_range(x1-1,y1,1)=1;
                        Visual_range(x1-1,y1,2)=cell_array(x,15);
                        Visual_range(x1-1,y1,3)=1;
                        Visual_range(x1-1,y1,4)=cell_label_1;
                        cell_array(x,1)=x1-1;
                        cell_array(x,14)=1;
                        Visual_range(x1,y1,3)=1;
                        int size1=cell_array.rows();
                        for (int cus=1;cus<=size1;cus++)
                        {
                            if(cell_array(cus,1)==x1 && cell_array(cus,5)==y1)
                            {
                                cell_array(cus,14)=1;
                            }
                        }
                    }
                    else if (order==3)
                    {
                        int cell_label_1=cell_label+1;
                        Visual_range(x1-1,y1+1,1)=1;
                        Visual_range(x1-1,y1+1,2)=cell_array(x,15);
                        Visual_range(x1-1,y1+1,3)=1;
                        Visual_range(x1-1,y1+1,4)=cell_label_1;
                        cell_array(x,1)=x1-1;
                        cell_array(x,5)=y1+1;
                        cell_array(x,14)=1;
                        Visual_range(x1,y1,3)=1;
                        int size1=cell_array.rows();
                        for (int cus=1;cus<=size1;cus++)
                        {
                            if(cell_array(cus,1)==x1 && cell_array(cus,5)==y1)
                            {
                                cell_array(cus,14)=1;
                            }
                        }
                    }
                    else if (order==4)
                    {
                        int cell_label_1=cell_label+1;
                        Visual_range(x1,y1+1,1)=1;
                        Visual_range(x1,y1+1,2)=cell_array(x,15);
                        Visual_range(x1,y1+1,3)=1;
                        Visual_range(x1,y1+1,4)=cell_label_1;
                        cell_array(x,5)=y1+1;
                        cell_array(x,14)=1;
                        Visual_range(x1,y1,3)=1;
                        int size1=cell_array.rows();
                        for (int cus=1;cus<=size1;cus++)
                        {
                            if(cell_array(cus,1)==x1 && cell_array(cus,5)==y1)
                            {
                                cell_array(cus,14)=1;
                            }
                        }
                    }
                    else if (order==5)
                    {
                        int cell_label_1=cell_label+1;
                        Visual_range(x1+1,y1+1,1)=1;
                        Visual_range(x1+1,y1+1,2)=cell_array(x,15);
                        Visual_range(x1+1,y1+1,3)=1;
                        Visual_range(x1+1,y1+1,4)=cell_label_1;
                        cell_array(x,1)=x1+1;
                        cell_array(x,5)=y1+1;
                        cell_array(x,14)=1;
                        Visual_range(x1,y1,3)=1;
                        int size1=cell_array.rows();
                        for (int cus=1;cus<=size1;cus++)
                        {
                            if(cell_array(cus,1)==x1 && cell_array(cus,5)==y1)
                            {
                                cell_array(cus,14)=1;
                            }
                        }
                    }
                    else if (order==6)
                    {
                        int cell_label_1=cell_label+1;
                        Visual_range(x1+1,y1,1)=1;
                        Visual_range(x1+1,y1,2)=cell_array(x,15);
                        Visual_range(x1+1,y1,3)=1;
                        Visual_range(x1+1,y1,4)=cell_label_1;
                        cell_array(x,1)=x1+1;
                        cell_array(x,14)=1;
                        Visual_range(x1,y1,3)=1;
                        int size1=cell_array.rows();
                        for (int cus=1;cus<=size1;cus++)
                        {
                            if(cell_array(cus,1)==x1 && cell_array(cus,5)==y1)
                            {
                                cell_array(cus,14)=1;
                            }
                        }
                    }
                    else if (order==7)
                    {
                        int cell_label_1=cell_label+1;
                        Visual_range(x1+1,y1-1,1)=1;
                        Visual_range(x1+1,y1-1,2)=cell_array(x,15);
                        Visual_range(x1+1,y1-1,3)=1;
                        Visual_range(x1+1,y1-1,4)=cell_label_1;
                        cell_array(x,1)=x1+1;
                        cell_array(x,5)=y1-1;
                        cell_array(x,14)=1;
                        Visual_range(x1,y1,3)=1;
                        int size1=cell_array.rows();
                        for (int cus=1;cus<=size1;cus++)
                        {
                            if(cell_array(cus,1)==x1 && cell_array(cus,5)==y1)
                            {
                                cell_array(cus,14)=1;
                            }
                        }
                    }
                    else if (order==8)
                    {
                        int cell_label_1=cell_label+1;
                        Visual_range(x1,y1-1,1)=1;
                        Visual_range(x1,y1-1,2)=cell_array(x,15);
                        Visual_range(x1,y1-1,3)=1;
                        Visual_range(x1,y1-1,4)=cell_label_1;
                        cell_array(x,5)=y1-1;
                        cell_array(x,14)=1;
                        Visual_range(x1,y1,3)=1;
                        int size1=cell_array.rows();
                        for (int cus=1;cus<=size1;cus++)
                        {
                            if(cell_array(cus,1)==x1 && cell_array(cus,5)==y1)
                            {
                                cell_array(cus,14)=1;
                            }
                        }
                    }
                    delete[] direction1;
                    direction1 = NULL;
                }
            }
        }
    }
    if (utralsmall==1)
    {
        for (int x=1; x<=C00; x++)
        {
            int stage_cor[4]={0};
            if (cell_array(x,1)>=100 && cell_array(x,5) >=100  && cell_array(x,1)<=Visual_range_x+100 && cell_array(x,5)<=Visual_range_y+100)
            {
                if(cell_array(x,14)==1)
                {
                    int x1 = cell_array(x,1);
                    int y1 = cell_array(x,5);
                    long loci_cor[9]={0};
                    int loci_direction=0;
                    for (int xx=x1-1;xx<=x1+1;xx++)
                    {
                        for (int yy=y1-1;yy<=y1+1;yy++)
                        {
                            loci_cor[loci_direction]=Visual_range(xx,yy,1);
                            loci_direction=loci_direction+1;
                        }
                    }
                    if (loci_cor[0]==0 && loci_cor[1]==0 && loci_cor[3]==0)
                    {
                        stage_cor[0]=1;
                    }
                    if (loci_cor[1]==0 && loci_cor[2]==0 && loci_cor[5]==0)
                    {
                        stage_cor[1]=2;
                    }
                    if (loci_cor[5]==0 && loci_cor[7]==0 && loci_cor[8]==0)
                    {
                        stage_cor[2]=3;
                    }
                    if (loci_cor[3]==0 && loci_cor[6]==0 && loci_cor[7]==0)
                    {
                        stage_cor[3]=4;
                    }
                    int size=0;
                    for(int a=0;a<4;a++)
                    {
                        if(stage_cor[a]!=0)
                        {
                            size=size+1;
                        }
                    }
                    if (size>0)
                    {
                        int *stage_cor_1=new int[size];
                        int num=0;
                        for(int a=0;a<4;a++)
                        {
                            if(stage_cor[a]!=0)
                            {
                                stage_cor_1[num]=stage_cor[a];
                                num=num+1;
                            }
                        }
                        shuffle(stage_cor_1,stage_cor_1+num,RNG);
                        int scor=stage_cor_1[0];
                        if (scor==1)
                        {
                            cell_array(x,14)=0;
                            cell_array(x,1)=x1-1;
                            cell_array(x,2)=x1-1;
                            cell_array(x,3)=x1;
                            cell_array(x,4)=x1;
                            cell_array(x,5)=y1-1;
                            cell_array(x,6)=y1;
                            cell_array(x,7)=y1;
                            cell_array(x,8)=y1-1;
                            Visual_range((int)cell_array(x,1),(int)cell_array(x,5),1)=1;
                            Visual_range((int)cell_array(x,2),(int)cell_array(x,6),1)=1;
                            Visual_range((int)cell_array(x,3),(int)cell_array(x,7),1)=1;
                            Visual_range((int)cell_array(x,4),(int)cell_array(x,8),1)=1;
                            Visual_range((int)cell_array(x,1),(int)cell_array(x,5),2)=cell_array(x,15);
                            Visual_range((int)cell_array(x,2),(int)cell_array(x,6),2)=cell_array(x,15);
                            Visual_range((int)cell_array(x,3),(int)cell_array(x,7),2)=cell_array(x,15);
                            Visual_range((int)cell_array(x,4),(int)cell_array(x,8),2)=cell_array(x,15);
                            Visual_range((int)cell_array(x,1),(int)cell_array(x,5),3)=0;
                            Visual_range((int)cell_array(x,2),(int)cell_array(x,6),3)=0;
                            Visual_range((int)cell_array(x,3),(int)cell_array(x,7),3)=0;
                            Visual_range((int)cell_array(x,4),(int)cell_array(x,8),3)=0;
                            Visual_range((int)cell_array(x,1),(int)cell_array(x,5),4)=Visual_range(x1,y1,4);
                            Visual_range((int)cell_array(x,2),(int)cell_array(x,6),4)=Visual_range(x1,y1,4);
                            Visual_range((int)cell_array(x,3),(int)cell_array(x,7),4)=Visual_range(x1,y1,4);
                            Visual_range((int)cell_array(x,4),(int)cell_array(x,8),4)=Visual_range(x1,y1,4);
                        }
                        else if (scor==2)
                        {
                            cell_array(x,14)=0;
                            cell_array(x,1)=x1-1;
                            cell_array(x,2)=x1-1;
                            cell_array(x,3)=x1;
                            cell_array(x,4)=x1;
                            cell_array(x,5)=y1;
                            cell_array(x,6)=y1+1;
                            cell_array(x,7)=y1+1;
                            cell_array(x,8)=y1;
                            Visual_range((int)cell_array(x,1),(int)cell_array(x,5),1)=1;
                            Visual_range((int)cell_array(x,2),(int)cell_array(x,6),1)=1;
                            Visual_range((int)cell_array(x,3),(int)cell_array(x,7),1)=1;
                            Visual_range((int)cell_array(x,4),(int)cell_array(x,8),1)=1;
                            Visual_range((int)cell_array(x,1),(int)cell_array(x,5),2)=cell_array(x,15);
                            Visual_range((int)cell_array(x,2),(int)cell_array(x,6),2)=cell_array(x,15);
                            Visual_range((int)cell_array(x,3),(int)cell_array(x,7),2)=cell_array(x,15);
                            Visual_range((int)cell_array(x,4),(int)cell_array(x,8),2)=cell_array(x,15);
                            Visual_range((int)cell_array(x,1),(int)cell_array(x,5),3)=0;
                            Visual_range((int)cell_array(x,2),(int)cell_array(x,6),3)=0;
                            Visual_range((int)cell_array(x,3),(int)cell_array(x,7),3)=0;
                            Visual_range((int)cell_array(x,4),(int)cell_array(x,8),3)=0;
                            Visual_range((int)cell_array(x,1),(int)cell_array(x,5),4)=Visual_range(x1,y1,4);
                            Visual_range((int)cell_array(x,2),(int)cell_array(x,6),4)=Visual_range(x1,y1,4);
                            Visual_range((int)cell_array(x,3),(int)cell_array(x,7),4)=Visual_range(x1,y1,4);
                            Visual_range((int)cell_array(x,4),(int)cell_array(x,8),4)=Visual_range(x1,y1,4);
                        }
                        else if (scor==3)
                        {
                            cell_array(x,14)=0;
                            cell_array(x,1)=x1;
                            cell_array(x,2)=x1;
                            cell_array(x,3)=x1+1;
                            cell_array(x,4)=x1+1;
                            cell_array(x,5)=y1;
                            cell_array(x,6)=y1+1;
                            cell_array(x,7)=y1+1;
                            cell_array(x,8)=y1;
                            Visual_range((int)cell_array(x,1),(int)cell_array(x,5),1)=1;
                            Visual_range((int)cell_array(x,2),(int)cell_array(x,6),1)=1;
                            Visual_range((int)cell_array(x,3),(int)cell_array(x,7),1)=1;
                            Visual_range((int)cell_array(x,4),(int)cell_array(x,8),1)=1;
                            Visual_range((int)cell_array(x,1),(int)cell_array(x,5),2)=cell_array(x,15);
                            Visual_range((int)cell_array(x,2),(int)cell_array(x,6),2)=cell_array(x,15);
                            Visual_range((int)cell_array(x,3),(int)cell_array(x,7),2)=cell_array(x,15);
                            Visual_range((int)cell_array(x,4),(int)cell_array(x,8),2)=cell_array(x,15);
                            Visual_range((int)cell_array(x,1),(int)cell_array(x,5),3)=0;
                            Visual_range((int)cell_array(x,2),(int)cell_array(x,6),3)=0;
                            Visual_range((int)cell_array(x,3),(int)cell_array(x,7),3)=0;
                            Visual_range((int)cell_array(x,4),(int)cell_array(x,8),3)=0;
                            Visual_range((int)cell_array(x,1),(int)cell_array(x,5),4)=Visual_range(x1,y1,4);
                            Visual_range((int)cell_array(x,2),(int)cell_array(x,6),4)=Visual_range(x1,y1,4);
                            Visual_range((int)cell_array(x,3),(int)cell_array(x,7),4)=Visual_range(x1,y1,4);
                            Visual_range((int)cell_array(x,4),(int)cell_array(x,8),4)=Visual_range(x1,y1,4);
                        }
                        else if (scor==4)
                        {
                            cell_array(x,14)=0;
                            cell_array(x,1)=x1;
                            cell_array(x,2)=x1;
                            cell_array(x,3)=x1+1;
                            cell_array(x,4)=x1+1;
                            cell_array(x,5)=y1-1;
                            cell_array(x,6)=y1;
                            cell_array(x,7)=y1;
                            cell_array(x,8)=y1-1;
                            Visual_range((int)cell_array(x,1),(int)cell_array(x,5),1)=1;
                            Visual_range((int)cell_array(x,2),(int)cell_array(x,6),1)=1;
                            Visual_range((int)cell_array(x,3),(int)cell_array(x,7),1)=1;
                            Visual_range((int)cell_array(x,4),(int)cell_array(x,8),1)=1;
                            Visual_range((int)cell_array(x,1),(int)cell_array(x,5),2)=cell_array(x,15);
                            Visual_range((int)cell_array(x,2),(int)cell_array(x,6),2)=cell_array(x,15);
                            Visual_range((int)cell_array(x,3),(int)cell_array(x,7),2)=cell_array(x,15);
                            Visual_range((int)cell_array(x,4),(int)cell_array(x,8),2)=cell_array(x,15);
                            Visual_range((int)cell_array(x,1),(int)cell_array(x,5),3)=0;
                            Visual_range((int)cell_array(x,2),(int)cell_array(x,6),3)=0;
                            Visual_range((int)cell_array(x,3),(int)cell_array(x,7),3)=0;
                            Visual_range((int)cell_array(x,4),(int)cell_array(x,8),3)=0;
                            Visual_range((int)cell_array(x,1),(int)cell_array(x,5),4)=Visual_range(x1,y1,4);
                            Visual_range((int)cell_array(x,2),(int)cell_array(x,6),4)=Visual_range(x1,y1,4);
                            Visual_range((int)cell_array(x,3),(int)cell_array(x,7),4)=Visual_range(x1,y1,4);
                            Visual_range((int)cell_array(x,4),(int)cell_array(x,8),4)=Visual_range(x1,y1,4);
                        }
                        delete[] stage_cor_1;
                        stage_cor_1 = NULL;
                    }
                }
            }
        }
    }
}
#endif /* stage_convert_hpp */
