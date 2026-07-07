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
#include "stateless_rng.hpp"
#include "cell_store.hpp"
using namespace blitz;
void stage_convert(int Visual_range_x, int Visual_range_y, Array<double,2> &cell_array, Array<long,3> &Visual_range, int &cell_label,int utralsmall,long rng_time_step)
{
    Range all = Range::all();
    int C00= cell_array.rows();
    for (int x=1; x<=C00; ++x)
    {
        int stage_cor[4]={0};
        int direction[8]={0};
        long cell_rng_id = (long)cell_array(x,15);
        if (cell_rng_id == 0)
        {
            cell_rng_id = x;
        }
        long rng_event = 100;
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
                    stateless_shuffle(stage_cor_1,stage_cor_1+num,cell_rng_id,rng_time_step,rng_event++);
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
                        stateless_shuffle(direction1, direction1+new_loci,cell_rng_id,rng_time_step,rng_event++);
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
            long cell_rng_id = (long)cell_array(x,15);
            if (cell_rng_id == 0)
            {
                cell_rng_id = x;
            }
            long rng_event = 1000;
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
                        stateless_shuffle(stage_cor_1,stage_cor_1+num,cell_rng_id,rng_time_step,rng_event++);
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

namespace stage_convert_detail
{
inline bool within_visual_range(int x, int y, int Visual_range_x, int Visual_range_y)
{
    return x >= 100 && y >= 100 && x <= Visual_range_x + 100 && y <= Visual_range_y + 100;
}

inline int choose_large_stage_square(int x1, int y1, const Array<long,3> &Visual_range, long cell_rng_id, long rng_time_step, long &rng_event)
{
    int candidates[4] = {0};
    int count = 0;
    long loci_cor[9] = {0};
    int loci_direction = 0;
    for (int xx = x1 - 1; xx <= x1 + 1; ++xx)
    {
        for (int yy = y1 - 1; yy <= y1 + 1; ++yy)
        {
            loci_cor[loci_direction] = Visual_range(xx, yy, 1);
            ++loci_direction;
        }
    }

    if (loci_cor[0] == 0 && loci_cor[1] == 0 && loci_cor[3] == 0)
    {
        candidates[count++] = 1;
    }
    if (loci_cor[1] == 0 && loci_cor[2] == 0 && loci_cor[5] == 0)
    {
        candidates[count++] = 2;
    }
    if (loci_cor[5] == 0 && loci_cor[7] == 0 && loci_cor[8] == 0)
    {
        candidates[count++] = 3;
    }
    if (loci_cor[3] == 0 && loci_cor[6] == 0 && loci_cor[7] == 0)
    {
        candidates[count++] = 4;
    }

    if (count == 0)
    {
        return 0;
    }
    stateless_shuffle(candidates, candidates + count, cell_rng_id, rng_time_step, rng_event++);
    return candidates[0];
}

inline void square_loci(int stage_square, int x1, int y1, int xs[4], int ys[4])
{
    if (stage_square == 1)
    {
        xs[0] = x1 - 1; xs[1] = x1 - 1; xs[2] = x1;     xs[3] = x1;
        ys[0] = y1 - 1; ys[1] = y1;     ys[2] = y1;     ys[3] = y1 - 1;
    }
    else if (stage_square == 2)
    {
        xs[0] = x1 - 1; xs[1] = x1 - 1; xs[2] = x1;     xs[3] = x1;
        ys[0] = y1;     ys[1] = y1 + 1; ys[2] = y1 + 1; ys[3] = y1;
    }
    else if (stage_square == 3)
    {
        xs[0] = x1;     xs[1] = x1;     xs[2] = x1 + 1; xs[3] = x1 + 1;
        ys[0] = y1;     ys[1] = y1 + 1; ys[2] = y1 + 1; ys[3] = y1;
    }
    else
    {
        xs[0] = x1;     xs[1] = x1;     xs[2] = x1 + 1; xs[3] = x1 + 1;
        ys[0] = y1 - 1; ys[1] = y1;     ys[2] = y1;     ys[3] = y1 - 1;
    }
}

inline void apply_large_stage_square(int row, CellStore &cells, Array<long,3> &Visual_range, int stage_square, int x1, int y1)
{
    int xs[4] = {0};
    int ys[4] = {0};
    square_loci(stage_square, x1, y1, xs, ys);

    cells.stage()[row - 1] = 0;
    cells.x1()[row - 1] = xs[0];
    cells.x2()[row - 1] = xs[1];
    cells.x3()[row - 1] = xs[2];
    cells.x4()[row - 1] = xs[3];
    cells.y1()[row - 1] = ys[0];
    cells.y2()[row - 1] = ys[1];
    cells.y3()[row - 1] = ys[2];
    cells.y4()[row - 1] = ys[3];

    const long cell_id = (long)cells.id()[row - 1];
    const long cell_label = Visual_range(x1, y1, 4);
    for (int idx = 0; idx < 4; ++idx)
    {
        Visual_range(xs[idx], ys[idx], 1) = 1;
        Visual_range(xs[idx], ys[idx], 2) = cell_id;
        Visual_range(xs[idx], ys[idx], 3) = 0;
        Visual_range(xs[idx], ys[idx], 4) = cell_label;
    }
}

inline int choose_small_stage_direction(int x1, int y1, const Array<long,3> &Visual_range, long cell_rng_id, long rng_time_step, long &rng_event)
{
    int candidates[8] = {0};
    int count = 0;
    if (Visual_range(x1 - 1, y1 - 1, 1) == 0) { candidates[count++] = 1; }
    if (Visual_range(x1 - 1, y1,     1) == 0) { candidates[count++] = 2; }
    if (Visual_range(x1 - 1, y1 + 1, 1) == 0) { candidates[count++] = 3; }
    if (Visual_range(x1,     y1 + 1, 1) == 0) { candidates[count++] = 4; }
    if (Visual_range(x1 + 1, y1 + 1, 1) == 0) { candidates[count++] = 5; }
    if (Visual_range(x1 + 1, y1,     1) == 0) { candidates[count++] = 6; }
    if (Visual_range(x1 + 1, y1 - 1, 1) == 0) { candidates[count++] = 7; }
    if (Visual_range(x1,     y1 - 1, 1) == 0) { candidates[count++] = 8; }

    if (count == 0)
    {
        return 0;
    }
    stateless_shuffle(candidates, candidates + count, cell_rng_id, rng_time_step, rng_event++);
    return candidates[0];
}

inline void direction_locus(int direction, int x1, int y1, int &x2, int &y2)
{
    static const int dx[9] = {0, -1, -1, -1, 0, 1, 1, 1, 0};
    static const int dy[9] = {0, -1, 0, 1, 1, 1, 0, -1, -1};
    x2 = x1 + dx[direction];
    y2 = y1 + dy[direction];
}

inline void apply_small_stage_direction(int row, CellStore &cells, Array<long,3> &Visual_range, int &cell_label, int direction, int x1, int y1)
{
    int x2 = x1;
    int y2 = y1;
    direction_locus(direction, x1, y1, x2, y2);

    const int cell_label_1 = cell_label + 1;
    Visual_range(x2, y2, 1) = 1;
    Visual_range(x2, y2, 2) = (long)cells.id()[row - 1];
    Visual_range(x2, y2, 3) = 1;
    Visual_range(x2, y2, 4) = cell_label_1;

    cells.x1()[row - 1] = x2;
    cells.y1()[row - 1] = y2;
    cells.stage()[row - 1] = 1;
    Visual_range(x1, y1, 3) = 1;

    const int row_count = cells.rows();
    CellStore::Column &xs = cells.x1();
    CellStore::Column &ys = cells.y1();
    CellStore::Column &stages = cells.stage();
    for (int other = 1; other <= row_count; ++other)
    {
        if ((int)xs[other - 1] == x1 && (int)ys[other - 1] == y1)
        {
            stages[other - 1] = 1;
        }
    }
}
}

inline void stage_convert(int Visual_range_x, int Visual_range_y, CellStore &cells, Array<long,3> &Visual_range, int &cell_label,int utralsmall,long rng_time_step)
{
    const int row_count = cells.rows();
    for (int row = 1; row <= row_count; ++row)
    {
        long cell_rng_id = (long)cells.id()[row - 1];
        if (cell_rng_id == 0)
        {
            cell_rng_id = row;
        }
        long rng_event = 100;
        const int x1 = (int)cells.x1()[row - 1];
        const int y1 = (int)cells.y1()[row - 1];
        if (!stage_convert_detail::within_visual_range(x1, y1, Visual_range_x, Visual_range_y))
        {
            continue;
        }

        if ((int)cells.stage()[row - 1] == 1)
        {
            int stage_square = stage_convert_detail::choose_large_stage_square(x1, y1, Visual_range, cell_rng_id, rng_time_step, rng_event);
            if (stage_square != 0)
            {
                stage_convert_detail::apply_large_stage_square(row, cells, Visual_range, stage_square, x1, y1);
            }
        }
        else if ((int)cells.stage()[row - 1] == 2)
        {
            int direction = stage_convert_detail::choose_small_stage_direction(x1, y1, Visual_range, cell_rng_id, rng_time_step, rng_event);
            if (direction != 0)
            {
                stage_convert_detail::apply_small_stage_direction(row, cells, Visual_range, cell_label, direction, x1, y1);
            }
        }
    }

    if (utralsmall == 1)
    {
        for (int row = 1; row <= row_count; ++row)
        {
            long cell_rng_id = (long)cells.id()[row - 1];
            if (cell_rng_id == 0)
            {
                cell_rng_id = row;
            }
            long rng_event = 1000;
            const int x1 = (int)cells.x1()[row - 1];
            const int y1 = (int)cells.y1()[row - 1];
            if (!stage_convert_detail::within_visual_range(x1, y1, Visual_range_x, Visual_range_y))
            {
                continue;
            }
            if ((int)cells.stage()[row - 1] == 1)
            {
                int stage_square = stage_convert_detail::choose_large_stage_square(x1, y1, Visual_range, cell_rng_id, rng_time_step, rng_event);
                if (stage_square != 0)
                {
                    stage_convert_detail::apply_large_stage_square(row, cells, Visual_range, stage_square, x1, y1);
                }
            }
        }
    }
}
#endif /* stage_convert_hpp */
