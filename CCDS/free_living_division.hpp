//
//  free_living_division.hpp
//  CCDS
//
//  Created by Tao Lee on 11/10/22.
//  Copyright © 2022 Tao Lee. All rights reserved.
//

#ifndef free_living_division_hpp
#define free_living_division_hpp


#include <stdio.h>
#include <random>
#include <memory>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <ctime>
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
#define BZ_THREADSAFE
#define BZ_THREADSAFE_USE_OPENMP
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "deltah_calculation.hpp"
#include "cell_type_transform.hpp"
#include <chrono>


using std::chrono::high_resolution_clock;
using namespace std;
using namespace blitz;
void free_living_division(int i, double max_growth_rate_r, double max_growth_rate_K, Array<float, 2> &cell_array, Array<float,2> cell_array_temp, Array<int, 3> &Visual_range, Array<int,2> cor_big_1, Array<int, 2> cor_big_1_change_shape, Array<int, 2> cor_small_1, Array<int, 2> proliferation_loci, Array<float, 2> cell_temp,int &cell_label, double &deltah,int utralsmall, double beta_distribution_alpha_for_normal_migration,double beta_distribution_beta_for_normal_migration,double migration_rate_K_mean,double uniup_K, double unilow_K,double sigmahatK,double muhatK,int &K_label,Array<int, 3> sub_visual,double beta_distribution_alpha, double beta_distribution_beta, double migration_rate_r_mean,double migration_rate_r_mean_quia,double beta_distribution_expected_for_normal_migration,Array<int,2> &cell_trace,Array<int,2> cell_trace_temp, int &cell_index,int &generation,int &r_label,int Col,double K_formation_rate)
{
    auto start = std::chrono::high_resolution_clock::now();
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 RNG(seed);
    Range all = Range::all();
 
    int pro_loci[8]={0};
    int pro_loci1[4]={0};
    int pro_loci2[4]={0};
    int pro_loci1_new[4]={0};
    int pro_loci2_new[4]={0};
    cor_big_1.resize(2, 16);
    cor_big_1=0;
    cor_big_1_change_shape.resize(2, 16);
    cor_big_1_change_shape=0;
    cor_small_1.resize(2, 8);
    cor_small_1=0;
    proliferation_loci.resize(2, 4);
    proliferation_loci=0;
    cell_temp.resize(1, Col);
    cell_temp=0;
    int x=(int)cell_array(i,1);
    int y=(int)cell_array(i,5);
    int A=x-2;
    int B=x-1;
    int C=x;
    int D=x+1;
    int E=x+2;
    int a=y-2;
    int b=y-1;
    int c=y;
    int d=y+1;
    int e=y+2;
    cor_big_1(1,1)=A;
    cor_big_1(2,1)=a;
    cor_big_1(1,2)=A;
    cor_big_1(2,2)=b;
    cor_big_1(1,3)=A;
    cor_big_1(2,3)=c;
    cor_big_1(1,4)=A;
    cor_big_1(2,4)=d;
    cor_big_1(1,5)=A;
    cor_big_1(2,5)=e;
    cor_big_1(1,6)=B;
    cor_big_1(2,6)=e;
    cor_big_1(1,7)=C;
    cor_big_1(2,7)=e;
    cor_big_1(1,8)=D;
    cor_big_1(2,8)=e;
    cor_big_1(1,9)=E;
    cor_big_1(2,9)=e;
    cor_big_1(1,10)=E;
    cor_big_1(2,10)=d;
    cor_big_1(1,11)=E;
    cor_big_1(2,11)=c;
    cor_big_1(1,12)=E;
    cor_big_1(2,12)=b;
    cor_big_1(1,13)=E;
    cor_big_1(2,13)=a;
    cor_big_1(1,14)=D;
    cor_big_1(2,14)=a;
    cor_big_1(1,15)=C;
    cor_big_1(2,15)=a;
    cor_big_1(1,16)=B;
    cor_big_1(2,16)=a;
    cor_big_1_change_shape(all,1)=B,b;
    cor_big_1_change_shape(all,2)=B,c;
    cor_big_1_change_shape(all,3)=B,d;
    cor_big_1_change_shape(all,4)=B,e;
    cor_big_1_change_shape(all,5)=C,e;
    cor_big_1_change_shape(all,6)=D,e;
    cor_big_1_change_shape(all,7)=E,e;
    cor_big_1_change_shape(all,8)=E,d;
    cor_big_1_change_shape(all,9)=E,c;
    cor_big_1_change_shape(all,10)=E,b;
    cor_big_1_change_shape(all,11)=D,b;
    cor_big_1_change_shape(all,12)=C,b;
    cor_big_1_change_shape(all,13)=C,c;
    cor_big_1_change_shape(all,14)=C,d;
    cor_big_1_change_shape(all,15)=D,d;
    cor_big_1_change_shape(all,16)=D,c;
    cor_small_1(all,1)=B,b;
    cor_small_1(all,2)=B,c;
    cor_small_1(all,3)=B,d;
    cor_small_1(all,4)=C,d;
    cor_small_1(all,5)=D,d;
    cor_small_1(all,6)=D,c;
    cor_small_1(all,7)=D,b;
    cor_small_1(all,8)=C,b;
    int cor_temp1[16]={0};
    
    const gsl_rng_type *T7;
    gsl_rng *r7;
    gsl_rng_env_setup();
    T7 = gsl_rng_ranlxs0;
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start);
    gsl_rng_default_seed = (duration.count());
    r7 = gsl_rng_alloc(T7);
    
    int cell_type=(int)cell_array(i,14);
    switch (cell_type)
    {
        case 0:
        {
            K_label=K_label+1;
            r_label=r_label+1;
            double growth_rate_inherent=cell_array(i,10);
            float X1=growth_rate_inherent*(1-0.05);
            float X2=growth_rate_inherent*(1+0.05);
            cell_array(i,10)=(X2-X1)*gsl_rng_uniform(r7)+X1;
            cell_array(i,13)=cell_array(i,13)+1;
            
            if(cell_array(i,9)==1)
            {
                cell_array(i,15)=r_label;
            }
            else if (cell_array(i,9)==1)
            {
                cell_array(i,15)=K_label;
            }
            
            
            cell_type_transform(cell_temp, beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration,migration_rate_K_mean,uniup_K,unilow_K, sigmahatK, muhatK,K_label,i,sub_visual,Visual_range,cell_array, beta_distribution_alpha, beta_distribution_beta, migration_rate_r_mean, migration_rate_r_mean_quia, beta_distribution_expected_for_normal_migration,r_label,K_formation_rate);
            
            cell_index=cell_index+1;
            generation=generation+1;
            
            cell_array(i,30)=cell_array(i,29);
            cell_array(i,29)=cell_index;
            
            cell_index=cell_index+1;
            cell_temp(1,29)=cell_index;
            cell_temp(1,30)=cell_array(i,30);
            
            cell_array(i,Col)=cell_array(i,13);
            cell_temp(1,Col)=cell_array(i,13);
            
            
            int cor_temp_length=0;
            for (int s=1; s<=16; s++)
            {
                int x1=cor_big_1(1,s);
                int y1=cor_big_1(2,s);
                if (Visual_range(x1,y1,1)==0 && Visual_range(x1,y1+1,1)==0 && Visual_range(x1+1,y1,1)==0 && Visual_range(x1+1,y1+1,1)==0)
                {
                    cor_temp1[cor_temp_length]=s;
                    cor_temp_length=cor_temp_length+1;
                }
            }
            int size_big=0;
            int cor_big_temp_1[cor_temp_length];
            for (int iiii=0;iiii<16;iiii++)
            {
                if (cor_temp1[iiii]!=0)
                {
                    cor_big_temp_1[size_big]=cor_temp1[iiii];
                    size_big=size_big+1;
                }
            }
            int loci_number=size_big;
            if (loci_number>0)
            {
                cell_label=cell_label+1;
                shuffle(cor_big_temp_1,cor_big_temp_1+size_big,RNG);
                //            gsl_ran_shuffle(r7, cor_big_temp_1, size_big, sizeof (int));
                int loci=cor_big_temp_1[0];
                cell_temp(1,1)=cor_big_1(1,loci);
                cell_temp(1,2)=cor_big_1(1,loci);
                cell_temp(1,3)=cor_big_1(1,loci)+1;
                cell_temp(1,4)=cor_big_1(1,loci)+1;
                cell_temp(1,5)=cor_big_1(2,loci);
                cell_temp(1,6)=cor_big_1(2,loci)+1;
                cell_temp(1,7)=cor_big_1(2,loci)+1;
                cell_temp(1,8)=cor_big_1(2,loci);
                cell_temp(1,14)=cell_array(i,14);
                cell_temp(1,22)=cell_array(i,22);
                cell_temp(1,24)=cell_array(i,24);
                cell_array(i,20)=0;
                cell_temp(1,20)=0;
                cell_temp(1,16)=0;
                cell_array(i,16)=0;
                int cell_index= (int)cell_temp(1,15);
                Visual_range((int)cell_temp(1,1),(int)cell_temp(1,5),1)=1;
                Visual_range((int)cell_temp(1,2),(int)cell_temp(1,6),1)=1;
                Visual_range((int)cell_temp(1,3),(int)cell_temp(1,7),1)=1;
                Visual_range((int)cell_temp(1,4),(int)cell_temp(1,8),1)=1;
                Visual_range((int)cell_temp(1,1),(int)cell_temp(1,5),2)=cell_index;
                Visual_range((int)cell_temp(1,2),(int)cell_temp(1,6),2)=cell_index;
                Visual_range((int)cell_temp(1,3),(int)cell_temp(1,7),2)=cell_index;
                Visual_range((int)cell_temp(1,4),(int)cell_temp(1,8),2)=cell_index;
                Visual_range((int)cell_temp(1,1),(int)cell_temp(1,5),3)=(int)cell_array(i,14);
                Visual_range((int)cell_temp(1,2),(int)cell_temp(1,6),3)=(int)cell_array(i,14);
                Visual_range((int)cell_temp(1,3),(int)cell_temp(1,7),3)=(int)cell_array(i,14);
                Visual_range((int)cell_temp(1,4),(int)cell_temp(1,8),3)=(int)cell_array(i,14);
                Visual_range((int)cell_temp(1,1),(int)cell_temp(1,5),4)=cell_label;
                Visual_range((int)cell_temp(1,2),(int)cell_temp(1,6),4)=cell_label;
                Visual_range((int)cell_temp(1,3),(int)cell_temp(1,7),4)=cell_label;
                Visual_range((int)cell_temp(1,4),(int)cell_temp(1,8),4)=cell_label;
            }
            else
            {
                if  (Visual_range(cor_big_1_change_shape(1,1),cor_big_1_change_shape(2,1),1)==0 && Visual_range(cor_big_1_change_shape(1,2),cor_big_1_change_shape(2,2),1)==0 && Visual_range(cor_big_1_change_shape(1,12),cor_big_1_change_shape(2,12),1)==0)
                {
                    pro_loci[0]=1;
                    pro_loci1[0]=1;
                }
                if (Visual_range(cor_big_1_change_shape(1,2),cor_big_1_change_shape(2,2),1)==0 && Visual_range(cor_big_1_change_shape(1,3),cor_big_1_change_shape(2,3),1)==0)
                {
                    pro_loci[1]=2;
                    pro_loci2[0]=2;
                }
                if (Visual_range(cor_big_1_change_shape(1,3),cor_big_1_change_shape(2,3),1)==0 && Visual_range(cor_big_1_change_shape(1,4),cor_big_1_change_shape(2,4),1)==0 && Visual_range(cor_big_1_change_shape(1,5),cor_big_1_change_shape(2,5),1)==0)
                {
                    pro_loci[2]=3;
                    pro_loci1[1]=3;
                }
                if (Visual_range(cor_big_1_change_shape(1,5),cor_big_1_change_shape(2,5),1)==0 && Visual_range(cor_big_1_change_shape(1,6),cor_big_1_change_shape(2,6),1)==0)
                {
                    pro_loci[3]=4;
                    pro_loci2[1]=4;
                }
                if (Visual_range(cor_big_1_change_shape(1,6),cor_big_1_change_shape(2,6),1)==0 && Visual_range(cor_big_1_change_shape(1,7),cor_big_1_change_shape(2,7),1)==0 && Visual_range(cor_big_1_change_shape(1,8),cor_big_1_change_shape(2,8),1)==0)
                {
                    pro_loci[4]=5;
                    pro_loci1[2]=5;
                }
                if (Visual_range(cor_big_1_change_shape(1,8),cor_big_1_change_shape(2,8),1)==0 && Visual_range(cor_big_1_change_shape(1,9),cor_big_1_change_shape(2,9),1)==0)
                {
                    pro_loci[5]=6;
                    pro_loci2[2]=6;
                }
                if (Visual_range(cor_big_1_change_shape(1,9),cor_big_1_change_shape(2,9),1)==0 && Visual_range(cor_big_1_change_shape(1,10),cor_big_1_change_shape(2,10),1)==0 && Visual_range(cor_big_1_change_shape(1,11),cor_big_1_change_shape(2,11),1)==0)
                {
                    pro_loci[6]=7;
                    pro_loci1[3]=7;
                }
                if (Visual_range(cor_big_1_change_shape(1,11),cor_big_1_change_shape(2,11),1)==0 && Visual_range(cor_big_1_change_shape(1,12),cor_big_1_change_shape(2,12),1)==0)
                {
                    pro_loci[7]=8;
                    pro_loci2[3]=8;
                }
                int siz=0;
                for(int aaa=0;aaa<8;aaa++)
                {
                    if(pro_loci[aaa]!=0)
                    {
                        siz=siz+1;
                    }
                }
                int *pro_loci_new=new int[siz];
                int num=0;
                for(int aaa=0;aaa<8;aaa++)
                {
                    if(pro_loci[aaa]!=0)
                    {
                        pro_loci_new[num]=pro_loci[aaa];
                        num=num+1;
                    }
                }
                int siz1=0;
                for(int aaa=0;aaa<4;aaa++)
                {
                    if(pro_loci1[aaa]!=0)
                    {
                        siz1=siz1+1;
                    }
                }
                int num1=0;
                for(int aaa=0;aaa<4;aaa++)
                {
                    if(pro_loci1[aaa]!=0)
                    {
                        pro_loci1_new[num1]=pro_loci1[aaa];
                        num1=num1+1;
                    }
                }
                int siz2=0;
                for(int aaa=0;aaa<4;aaa++)
                {
                    if(pro_loci2[aaa]!=0)
                    {
                        siz2=siz2+1;
                    }
                }
                int num2=0;
                for(int aaa=0;aaa<4;aaa++)
                {
                    if(pro_loci[aaa]!=0)
                    {
                        pro_loci2_new[num2]=pro_loci2[aaa];
                        num2=num2+1;
                    }
                }
                int pro_loci_number=siz;
                int pro_loci_number1=siz1;
                int pro_loci_number2=siz2;
                if (pro_loci_number>1)
                {
                    cell_array(i,20)=0;
                    cell_temp(1,20)=0;
                    cell_temp(1,16)=0;
                    cell_array(i,16)=0;
                    int cell_label_1=Visual_range(x,y,4);
                    int cellstage=Visual_range(x,y,3);
                    cell_label=cell_label+1;
                    if (pro_loci_number1 > 1 && pro_loci_number2 >1)
                    {
                        int random_pro_loci[2]={1,2};
                        shuffle(random_pro_loci,random_pro_loci+2,RNG);
                        //                    gsl_ran_shuffle(r7, random_pro_loci, 2, sizeof (int));
                        
                        int LociNumber=random_pro_loci[0];
                        switch (LociNumber)
                        {
                            case 1:
                            {
                                int sizep1=num1;
                                int random_pro_loci1[sizep1];
                                for(int aa=0;aa<sizep1;aa++)
                                {
                                    random_pro_loci1[aa]=pro_loci1_new[aa];
                                }
                                shuffle(random_pro_loci1,random_pro_loci1+sizep1,RNG);
                                //                            gsl_ran_shuffle(r7, random_pro_loci1, sizep1, sizeof (int));
                                int proloci1=random_pro_loci1[0];
                                int proloci2=random_pro_loci1[1];
                                switch (proloci1)
                                {
                                    case 1:
                                    {
                                        Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                        for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                        {
                                            int cor_cell_y=cor_cell+4;
                                            cell_array(i,cor_cell)=cell_array(i,cor_cell)-1;
                                            cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)-1;
                                        }
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                        break;
                                    }
                                    case 3:
                                    {
                                        Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                        for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                        {
                                            int cor_cell_y=cor_cell+4;
                                            cell_array(i,cor_cell)=cell_array(i,cor_cell)-1;
                                            cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)+1;
                                        }
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                        break;
                                    }
                                    case 5:
                                    {
                                        Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                        for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                        {
                                            int cor_cell_y=cor_cell+4;
                                            cell_array(i,cor_cell)=cell_array(i,cor_cell)+1;
                                            cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)+1;
                                        }
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                        break;
                                    }
                                    case 7:
                                    {
                                        Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                        for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                        {
                                            int cor_cell_y=cor_cell+4;
                                            cell_array(i,cor_cell)=cell_array(i,cor_cell)+1;
                                            cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)-1;
                                        }
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                        break;
                                    }
                                }
                                switch (proloci2)
                                {
                                    case 1:
                                    {
                                        Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                        cell_temp(1,1)=x-1;
                                        cell_temp(1,5)=y-1;
                                        cell_temp(1,2)=cell_temp(1,1);
                                        cell_temp(1,3)=cell_temp(1,1)+1;
                                        cell_temp(1,4)=cell_temp(1,1)+1;
                                        cell_temp(1,6)=cell_temp(1,5)+1;
                                        cell_temp(1,7)=cell_temp(1,5)+1;
                                        cell_temp(1,8)=cell_temp(1,5);
                                        cell_temp(1,14)=cell_array(i,14);
                                        cell_temp(1,22)=cell_array(i,22);
                                        cell_temp(1,24)=cell_array(i,24);
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                        break;
                                    }
                                    case 3:
                                    {
                                        Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                        cell_temp(1,1)=x-1;
                                        cell_temp(1,5)=y+1;
                                        cell_temp(1,2)=cell_temp(1,1);
                                        cell_temp(1,3)=cell_temp(1,1)+1;
                                        cell_temp(1,4)=cell_temp(1,1)+1;
                                        cell_temp(1,6)=cell_temp(1,5)+1;
                                        cell_temp(1,7)=cell_temp(1,5)+1;
                                        cell_temp(1,8)=cell_temp(1,5);
                                        cell_temp(1,14)=cell_array(i,14);
                                        cell_temp(1,22)=cell_array(i,22);
                                        cell_temp(1,24)=cell_array(i,24);
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                        break;
                                    }
                                    case 5:
                                    {
                                        Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                        cell_temp(1,1)=x+1;
                                        cell_temp(1,5)=y+1;
                                        cell_temp(1,2)=cell_temp(1,1);
                                        cell_temp(1,3)=cell_temp(1,1)+1;
                                        cell_temp(1,4)=cell_temp(1,1)+1;
                                        cell_temp(1,6)=cell_temp(1,5)+1;
                                        cell_temp(1,7)=cell_temp(1,5)+1;
                                        cell_temp(1,8)=cell_temp(1,5);
                                        cell_temp(1,14)=cell_array(i,14);
                                        cell_temp(1,22)=cell_array(i,22);
                                        cell_temp(1,24)=cell_array(i,24);
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                        break;
                                    }
                                    case 7:
                                    {
                                        Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                        cell_temp(1,1)=x+1;
                                        cell_temp(1,5)=y-1;
                                        cell_temp(1,2)=cell_temp(1,1);
                                        cell_temp(1,3)=cell_temp(1,1)+1;
                                        cell_temp(1,4)=cell_temp(1,1)+1;
                                        cell_temp(1,6)=cell_temp(1,5)+1;
                                        cell_temp(1,7)=cell_temp(1,5)+1;
                                        cell_temp(1,8)=cell_temp(1,5);
                                        cell_temp(1,14)=cell_array(i,14);
                                        cell_temp(1,22)=cell_array(i,22);
                                        cell_temp(1,24)=cell_array(i,24);
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                        break;
                                    }
                                }
                                break;
                            }
                            case 2:
                            {
                                int sizep2=num2;
                                int random_pro_loci2[sizep2];
                                for(int aa=0;aa<sizep2;aa++)
                                {
                                    random_pro_loci2[aa]=pro_loci2_new[aa];
                                }
                                shuffle(random_pro_loci2,random_pro_loci2+sizep2,RNG);
                                //                            gsl_ran_shuffle(r7, random_pro_loci2, sizep2, sizeof (int));
                                int proloci1=random_pro_loci2[0];
                                int proloci2=random_pro_loci2[1];
                                switch (proloci1)
                                {
                                    case 2:
                                    {
                                        Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                        for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                        {
                                            cell_array(i,cor_cell)=cell_array(i,cor_cell)-1;
                                        }
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                        break;
                                    }
                                    case 4:
                                    {
                                        Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                        for (int cor_cell_y=5; cor_cell_y<=8;cor_cell_y++)
                                        {
                                            cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)+1;
                                        }
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                        break;
                                    }
                                    case 6:
                                    {
                                        Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                        for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                        {
                                            cell_array(i,cor_cell)=cell_array(i,cor_cell)+1;
                                        }
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                        break;
                                    }
                                    case 8:
                                    {
                                        Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                        for (int cor_cell_y=5; cor_cell_y<=8;cor_cell_y++)
                                        {
                                            cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)-1;
                                        }
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                        Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                        break;
                                    }
                                }
                                switch (proloci2)
                                {
                                    case 2:
                                    {
                                        Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                        cell_temp(1,1)=x-1;
                                        cell_temp(1,5)=y;
                                        cell_temp(1,2)=cell_temp(1,1);
                                        cell_temp(1,3)=cell_temp(1,1)+1;
                                        cell_temp(1,4)=cell_temp(1,1)+1;
                                        cell_temp(1,6)=cell_temp(1,5)+1;
                                        cell_temp(1,7)=cell_temp(1,5)+1;
                                        cell_temp(1,8)=cell_temp(1,5);
                                        cell_temp(1,14)=cell_array(i,14);
                                        cell_temp(1,22)=cell_array(i,22);
                                        cell_temp(1,24)=cell_array(i,24);
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                        break;
                                    }
                                    case 4:
                                    {
                                        Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                        cell_temp(1,1)=x;
                                        cell_temp(1,5)=y+1;
                                        cell_temp(1,2)=cell_temp(1,1);
                                        cell_temp(1,3)=cell_temp(1,1)+1;
                                        cell_temp(1,4)=cell_temp(1,1)+1;
                                        cell_temp(1,6)=cell_temp(1,5)+1;
                                        cell_temp(1,7)=cell_temp(1,5)+1;
                                        cell_temp(1,8)=cell_temp(1,5);
                                        cell_temp(1,14)=cell_array(i,14);
                                        cell_temp(1,22)=cell_array(i,22);
                                        cell_temp(1,24)=cell_array(i,24);
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                        break;
                                    }
                                    case 6:
                                    {
                                        Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                        cell_temp(1,1)=x+1;
                                        cell_temp(1,5)=y;
                                        cell_temp(1,2)=cell_temp(1,1);
                                        cell_temp(1,3)=cell_temp(1,1)+1;
                                        cell_temp(1,4)=cell_temp(1,1)+1;
                                        cell_temp(1,6)=cell_temp(1,5)+1;
                                        cell_temp(1,7)=cell_temp(1,5)+1;
                                        cell_temp(1,8)=cell_temp(1,5);
                                        cell_temp(1,14)=cell_array(i,14);
                                        cell_temp(1,22)=cell_array(i,22);
                                        cell_temp(1,24)=cell_array(i,24);
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                        break;
                                    }
                                    case 8:
                                    {
                                        Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                        cell_temp(1,1)=x;
                                        cell_temp(1,5)=y-1;
                                        cell_temp(1,2)=cell_temp(1,1);
                                        cell_temp(1,3)=cell_temp(1,1)+1;
                                        cell_temp(1,4)=cell_temp(1,1)+1;
                                        cell_temp(1,6)=cell_temp(1,5)+1;
                                        cell_temp(1,7)=cell_temp(1,5)+1;
                                        cell_temp(1,8)=cell_temp(1,5);
                                        cell_temp(1,14)=cell_array(i,14);
                                        cell_temp(1,22)=cell_array(i,22);
                                        cell_temp(1,24)=cell_array(i,24);
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                        Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                        break;
                                    }
                                }
                                break;
                            }
                        }
                    }
                    else if (pro_loci_number1 > 1 && pro_loci_number2 <2)
                    {
                        int sizep1=num1;;
                        int random_pro_loci1[sizep1];
                        for(int aa=0;aa<sizep1;aa++)
                        {
                            random_pro_loci1[aa]=pro_loci1_new[aa];
                        }
                        shuffle(random_pro_loci1,random_pro_loci1+sizep1,RNG);
                        //                    gsl_ran_shuffle(r7, random_pro_loci1, sizep1, sizeof (int));
                        int proloci1=random_pro_loci1[0];
                        int proloci2=random_pro_loci1[1];
                        
                        switch (proloci1)
                        {
                                
                            case 1:
                            {
                                Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                {
                                    int cor_cell_y=cor_cell+4;
                                    cell_array(i,cor_cell)=cell_array(i,cor_cell)-1;
                                    cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)-1;
                                }
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                break;
                            }
                            case 3:
                            {
                                Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                {
                                    int cor_cell_y=cor_cell+4;
                                    cell_array(i,cor_cell)=cell_array(i,cor_cell)-1;
                                    cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)+1;
                                }
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                break;
                            }
                            case 5:
                            {
                                Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                {
                                    int cor_cell_y=cor_cell+4;
                                    cell_array(i,cor_cell)=cell_array(i,cor_cell)+1;
                                    cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)+1;
                                }
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                break;
                            }
                            case 7:
                            {
                                Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                {
                                    int cor_cell_y=cor_cell+4;
                                    cell_array(i,cor_cell)=cell_array(i,cor_cell)+1;
                                    cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)-1;
                                }
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                break;
                            }
                        }
                        switch (proloci2)
                        {
                            case 1:
                            {
                                Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                cell_temp(1,1)=x-1;
                                cell_temp(1,5)=y-1;
                                cell_temp(1,2)=cell_temp(1,1);
                                cell_temp(1,3)=cell_temp(1,1)+1;
                                cell_temp(1,4)=cell_temp(1,1)+1;
                                cell_temp(1,6)=cell_temp(1,5)+1;
                                cell_temp(1,7)=cell_temp(1,5)+1;
                                cell_temp(1,8)=cell_temp(1,5);
                                cell_temp(1,14)=cell_array(i,14);
                                cell_temp(1,22)=cell_array(i,22);
                                cell_temp(1,24)=cell_array(i,24);
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                break;
                            }
                            case 3:
                            {
                                Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                cell_temp(1,1)=x-1;
                                cell_temp(1,5)=y+1;
                                cell_temp(1,2)=cell_temp(1,1);
                                cell_temp(1,3)=cell_temp(1,1)+1;
                                cell_temp(1,4)=cell_temp(1,1)+1;
                                cell_temp(1,6)=cell_temp(1,5)+1;
                                cell_temp(1,7)=cell_temp(1,5)+1;
                                cell_temp(1,8)=cell_temp(1,5);
                                cell_temp(1,14)=cell_array(i,14);
                                cell_temp(1,22)=cell_array(i,22);
                                cell_temp(1,24)=cell_array(i,24);
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                break;
                            }
                            case 5:
                            {
                                Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                cell_temp(1,1)=x+1;
                                cell_temp(1,5)=y+1;
                                cell_temp(1,2)=cell_temp(1,1);
                                cell_temp(1,3)=cell_temp(1,1)+1;
                                cell_temp(1,4)=cell_temp(1,1)+1;
                                cell_temp(1,6)=cell_temp(1,5)+1;
                                cell_temp(1,7)=cell_temp(1,5)+1;
                                cell_temp(1,8)=cell_temp(1,5);
                                cell_temp(1,14)=cell_array(i,14);
                                cell_temp(1,22)=cell_array(i,22);
                                cell_temp(1,24)=cell_array(i,24);
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                break;
                            }
                            case 7:
                            {
                                Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                cell_temp(1,1)=x+1;
                                cell_temp(1,5)=y-1;
                                cell_temp(1,2)=cell_temp(1,1);
                                cell_temp(1,3)=cell_temp(1,1)+1;
                                cell_temp(1,4)=cell_temp(1,1)+1;
                                cell_temp(1,6)=cell_temp(1,5)+1;
                                cell_temp(1,7)=cell_temp(1,5)+1;
                                cell_temp(1,8)=cell_temp(1,5);
                                cell_temp(1,14)=cell_array(i,14);
                                cell_temp(1,22)=cell_array(i,22);
                                cell_temp(1,24)=cell_array(i,24);
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                break;
                            }
                        }
                    }
                    else if (pro_loci_number1 < 2 && pro_loci_number2 >1)
                    {
                        int sizep2=num2;
                        int random_pro_loci2[sizep2];
                        for(int aa=0;aa<sizep2;aa++)
                        {
                            random_pro_loci2[aa]=pro_loci2_new[aa];
                        }
                        shuffle(random_pro_loci2,random_pro_loci2+sizep2,RNG);
                        //                    gsl_ran_shuffle(r7, random_pro_loci2, sizep2, sizeof (int));
                        int proloci1=random_pro_loci2[0];
                        int proloci2=random_pro_loci2[1];
                        
                        switch (proloci1)
                        {
                            case 2:
                            {
                                Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                {
                                    cell_array(i,cor_cell)=cell_array(i,cor_cell)-1;
                                }
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                break;
                            }
                            case 4:
                            {
                                Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                for (int cor_cell_y=5; cor_cell_y<=8;cor_cell_y++)
                                {
                                    cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)+1;
                                }
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                break;
                            }
                            case 6:
                            {
                                Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                {
                                    cell_array(i,cor_cell)=cell_array(i,cor_cell)+1;
                                }
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                break;
                            }
                            case 8:
                            {
                                Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                for (int cor_cell_y=5; cor_cell_y<=8;cor_cell_y++)
                                {
                                    cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)-1;
                                }
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                break;
                            }
                        }
                        switch (proloci2)
                        {
                            case 2:
                            {
                                Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                cell_temp(1,1)=x-1;
                                cell_temp(1,5)=y;
                                cell_temp(1,2)=cell_temp(1,1);
                                cell_temp(1,3)=cell_temp(1,1)+1;
                                cell_temp(1,4)=cell_temp(1,1)+1;
                                cell_temp(1,6)=cell_temp(1,5)+1;
                                cell_temp(1,7)=cell_temp(1,5)+1;
                                cell_temp(1,8)=cell_temp(1,5);
                                cell_temp(1,14)=cell_array(i,14);
                                cell_temp(1,22)=cell_array(i,22);
                                cell_temp(1,24)=cell_array(i,24);
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                break;
                                
                            }
                            case 4:
                            {
                                Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                cell_temp(1,1)=x;
                                cell_temp(1,5)=y+1;
                                cell_temp(1,2)=cell_temp(1,1);
                                cell_temp(1,3)=cell_temp(1,1)+1;
                                cell_temp(1,4)=cell_temp(1,1)+1;
                                cell_temp(1,6)=cell_temp(1,5)+1;
                                cell_temp(1,7)=cell_temp(1,5)+1;
                                cell_temp(1,8)=cell_temp(1,5);
                                cell_temp(1,14)=cell_array(i,14);
                                cell_temp(1,22)=cell_array(i,22);
                                cell_temp(1,24)=cell_array(i,24);
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                break;
                            }
                            case 6:
                            {
                                Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                cell_temp(1,1)=x+1;
                                cell_temp(1,5)=y;
                                cell_temp(1,2)=cell_temp(1,1);
                                cell_temp(1,3)=cell_temp(1,1)+1;
                                cell_temp(1,4)=cell_temp(1,1)+1;
                                cell_temp(1,6)=cell_temp(1,5)+1;
                                cell_temp(1,7)=cell_temp(1,5)+1;
                                cell_temp(1,8)=cell_temp(1,5);
                                cell_temp(1,14)=cell_array(i,14);
                                cell_temp(1,22)=cell_array(i,22);
                                cell_temp(1,24)=cell_array(i,24);
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                break;
                            }
                            case 8:
                            {
                                Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                cell_temp(1,1)=x;
                                cell_temp(1,5)=y-1;
                                cell_temp(1,2)=cell_temp(1,1);
                                cell_temp(1,3)=cell_temp(1,1)+1;
                                cell_temp(1,4)=cell_temp(1,1)+1;
                                cell_temp(1,6)=cell_temp(1,5)+1;
                                cell_temp(1,7)=cell_temp(1,5)+1;
                                cell_temp(1,8)=cell_temp(1,5);
                                cell_temp(1,14)=cell_array(i,14);
                                cell_temp(1,22)=cell_array(i,22);
                                cell_temp(1,24)=cell_array(i,24);
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                break;
                            }
                        }
                    }
                    else if (pro_loci_number1 ==1 && pro_loci_number2 ==1)
                    {
                        if (abs(pro_loci_number2-pro_loci_number1)>=2 && pro_loci_number2-pro_loci_number1!=7)
                        {
                            int random_pro_loci[2]={1,2};
                            shuffle(random_pro_loci,random_pro_loci+2,RNG);
                            //                        gsl_ran_shuffle(r7, random_pro_loci, 2, sizeof (int));
                            int LociNumber=random_pro_loci[0];
                            switch (LociNumber)
                            {
                                    
                                case 1:
                                {
                                    int proloci1=pro_loci1_new[0];
                                    int proloci2=pro_loci2_new[0];
                                    switch (proloci1)
                                    {
                                        case 1:
                                        {
                                            Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                            for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                            {
                                                int cor_cell_y=cor_cell+4;
                                                cell_array(i,cor_cell)=cell_array(i,cor_cell)-1;
                                                cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)-1;
                                            }
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                            break;
                                        }
                                        case 3:
                                        {
                                            Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                            for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                            {
                                                int cor_cell_y=cor_cell+4;
                                                cell_array(i,cor_cell)=cell_array(i,cor_cell)-1;
                                                cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)+1;
                                            }
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                            break;
                                        }
                                        case 5:
                                        {
                                            Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                            for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                            {
                                                int cor_cell_y=cor_cell+4;
                                                cell_array(i,cor_cell)=cell_array(i,cor_cell)+1;
                                                cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)+1;
                                            }
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                            break;
                                        }
                                        case 7:
                                        {
                                            Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                            for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                            {
                                                int cor_cell_y=cor_cell+4;
                                                cell_array(i,cor_cell)=cell_array(i,cor_cell)+1;
                                                cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)-1;
                                            }
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                            break;
                                        }
                                    }
                                    switch (proloci2)
                                    {
                                        case 2:
                                        {
                                            Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                            cell_temp(1,1)=x-1;
                                            cell_temp(1,5)=y;
                                            cell_temp(1,2)=cell_temp(1,1);
                                            cell_temp(1,3)=cell_temp(1,1)+1;
                                            cell_temp(1,4)=cell_temp(1,1)+1;
                                            cell_temp(1,6)=cell_temp(1,5)+1;
                                            cell_temp(1,7)=cell_temp(1,5)+1;
                                            cell_temp(1,8)=cell_temp(1,5);
                                            cell_temp(1,14)=cell_array(i,14);
                                            cell_temp(1,22)=cell_array(i,22);
                                            cell_temp(1,24)=cell_array(i,24);
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                            break;
                                            
                                        }
                                        case 4:
                                        {
                                            Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                            cell_temp(1,1)=x;
                                            cell_temp(1,5)=y+1;
                                            cell_temp(1,2)=cell_temp(1,1);
                                            cell_temp(1,3)=cell_temp(1,1)+1;
                                            cell_temp(1,4)=cell_temp(1,1)+1;
                                            cell_temp(1,6)=cell_temp(1,5)+1;
                                            cell_temp(1,7)=cell_temp(1,5)+1;
                                            cell_temp(1,8)=cell_temp(1,5);
                                            cell_temp(1,14)=cell_array(i,14);
                                            cell_temp(1,22)=cell_array(i,22);
                                            cell_temp(1,24)=cell_array(i,24);
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                            break;
                                        }
                                        case 6:
                                        {
                                            Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                            cell_temp(1,1)=x+1;
                                            cell_temp(1,5)=y;
                                            cell_temp(1,2)=cell_temp(1,1);
                                            cell_temp(1,3)=cell_temp(1,1)+1;
                                            cell_temp(1,4)=cell_temp(1,1)+1;
                                            cell_temp(1,6)=cell_temp(1,5)+1;
                                            cell_temp(1,7)=cell_temp(1,5)+1;
                                            cell_temp(1,8)=cell_temp(1,5);
                                            cell_temp(1,14)=cell_array(i,14);
                                            cell_temp(1,22)=cell_array(i,22);
                                            cell_temp(1,24)=cell_array(i,24);
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                            break;
                                        }
                                        case 8:
                                        {
                                            Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                            cell_temp(1,1)=x;
                                            cell_temp(1,5)=y-1;
                                            cell_temp(1,2)=cell_temp(1,1);
                                            cell_temp(1,3)=cell_temp(1,1)+1;
                                            cell_temp(1,4)=cell_temp(1,1)+1;
                                            cell_temp(1,6)=cell_temp(1,5)+1;
                                            cell_temp(1,7)=cell_temp(1,5)+1;
                                            cell_temp(1,8)=cell_temp(1,5);
                                            cell_temp(1,14)=cell_array(i,14);
                                            cell_temp(1,22)=cell_array(i,22);
                                            cell_temp(1,24)=cell_array(i,24);
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                            break;
                                        }
                                    }
                                    break;
                                }
                                case 2:
                                {
                                    int proloci1=pro_loci2_new[0];
                                    int proloci2=pro_loci1_new[0];
                                    switch (proloci1)
                                    {
                                        case 2:
                                        {
                                            Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                            for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                            {
                                                cell_array(i,cor_cell)=cell_array(i,cor_cell)-1;
                                            }
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                            break;
                                        }
                                        case 4:
                                        {
                                            Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                            for (int cor_cell_y=5; cor_cell_y<=8;cor_cell_y++)
                                            {
                                                cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)+1;
                                            }
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                            break;
                                        }
                                        case 6:
                                        {
                                            Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                            for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                            {
                                                cell_array(i,cor_cell)=cell_array(i,cor_cell)+1;
                                            }
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                            break;
                                        }
                                        case 8:
                                        {
                                            Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                            for (int cor_cell_y=5; cor_cell_y<=8;cor_cell_y++)
                                            {
                                                cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)-1;
                                            }
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=(int)cell_array(i,15);
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                                            Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                                            break;
                                        }
                                    }
                                    switch (proloci2)
                                    {
                                        case 1:
                                        {
                                            Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                            cell_temp(1,1)=x-1;
                                            cell_temp(1,5)=y-1;
                                            cell_temp(1,2)=cell_temp(1,1);
                                            cell_temp(1,3)=cell_temp(1,1)+1;
                                            cell_temp(1,4)=cell_temp(1,1)+1;
                                            cell_temp(1,6)=cell_temp(1,5)+1;
                                            cell_temp(1,7)=cell_temp(1,5)+1;
                                            cell_temp(1,8)=cell_temp(1,5);
                                            cell_temp(1,14)=cell_array(i,14);
                                            cell_temp(1,22)=cell_array(i,22);
                                            cell_temp(1,24)=cell_array(i,24);
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                            break;
                                        }
                                        case 3:
                                        {
                                            Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                            cell_temp(1,1)=x-1;
                                            cell_temp(1,5)=y+1;
                                            cell_temp(1,2)=cell_temp(1,1);
                                            cell_temp(1,3)=cell_temp(1,1)+1;
                                            cell_temp(1,4)=cell_temp(1,1)+1;
                                            cell_temp(1,6)=cell_temp(1,5)+1;
                                            cell_temp(1,7)=cell_temp(1,5)+1;
                                            cell_temp(1,8)=cell_temp(1,5);
                                            cell_temp(1,14)=cell_array(i,14);
                                            cell_temp(1,22)=cell_array(i,22);
                                            cell_temp(1,24)=cell_array(i,24);
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                            break;
                                        }
                                        case 5:
                                        {
                                            Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                            cell_temp(1,1)=x+1;
                                            cell_temp(1,5)=y+1;
                                            cell_temp(1,2)=cell_temp(1,1);
                                            cell_temp(1,3)=cell_temp(1,1)+1;
                                            cell_temp(1,4)=cell_temp(1,1)+1;
                                            cell_temp(1,6)=cell_temp(1,5)+1;
                                            cell_temp(1,7)=cell_temp(1,5)+1;
                                            cell_temp(1,8)=cell_temp(1,5);
                                            cell_temp(1,14)=cell_array(i,14);
                                            cell_temp(1,22)=cell_array(i,22);
                                            cell_temp(1,24)=cell_array(i,24);
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                            break;
                                        }
                                        case 7:
                                        {
                                            Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                                            cell_temp(1,1)=x+1;
                                            cell_temp(1,5)=y-1;
                                            cell_temp(1,2)=cell_temp(1,1);
                                            cell_temp(1,3)=cell_temp(1,1)+1;
                                            cell_temp(1,4)=cell_temp(1,1)+1;
                                            cell_temp(1,6)=cell_temp(1,5)+1;
                                            cell_temp(1,7)=cell_temp(1,5)+1;
                                            cell_temp(1,8)=cell_temp(1,5);
                                            cell_temp(1,14)=cell_array(i,14);
                                            cell_temp(1,22)=cell_array(i,22);
                                            cell_temp(1,24)=cell_array(i,24);
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),1)=1;
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),2)=(int)cell_temp(1,15);
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),3)=cellstage;
                                            Visual_range(Range(cell_temp(1,1),cell_temp(1,1)+1),Range(cell_temp(1,5),cell_temp(1,5)+1),4)=cell_label;
                                            break;
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                    }
                    else
                    {
                        Array<int, 2> cor_pro_1_1(2,16,FortranArray<2>());
                        cor_pro_1_1=0;
                        int ssss=1;
                        for (int sss=1;sss<=16;sss++)
                        {
                            if (Visual_range(cor_big_1_change_shape(1,sss),cor_big_1_change_shape(2,sss),1)==0)
                            {
                                cor_pro_1_1(1,ssss)=cor_big_1_change_shape(1,sss);
                                cor_pro_1_1(2,ssss)=cor_big_1_change_shape(2,sss);
                                ssss=ssss+1;
                            }
                        }
                        int cor_pro_1_1_nozero_length=0;
                        for (int cor_pro_1_1_locus=1;cor_pro_1_1_locus<=16;cor_pro_1_1_locus++)
                        {
                            if (cor_pro_1_1(1,cor_pro_1_1_locus)!=0)
                            {
                                cor_pro_1_1_nozero_length=cor_pro_1_1_nozero_length+1;
                            }
                        }
                        Array<int, 2> cor_pro_1(2,cor_pro_1_1_nozero_length+4,FortranArray<2>());
                        cor_pro_1=0;
                        int nzl=1;
                        for (int aa=1;aa<=16;aa++ )
                        {
                            if (cor_pro_1_1(1,aa)!=0)
                            {
                                cor_pro_1(all,nzl)=cor_pro_1_1(all,aa);
                                nzl=nzl+1;
                            }
                        }
                        cor_pro_1(1,cor_pro_1_1_nozero_length+1)=x;
                        cor_pro_1(2,cor_pro_1_1_nozero_length+1)=y;
                        cor_pro_1(1,cor_pro_1_1_nozero_length+2)=x;
                        cor_pro_1(2,cor_pro_1_1_nozero_length+2)=y+1;
                        cor_pro_1(1,cor_pro_1_1_nozero_length+3)=x+1;
                        cor_pro_1(2,cor_pro_1_1_nozero_length+3)=y+1;
                        cor_pro_1(1,cor_pro_1_1_nozero_length+4)=x+1;
                        cor_pro_1(2,cor_pro_1_1_nozero_length+4)=y;
                        int cor_pro_1_length=cor_pro_1.columns();
                        int *random_pro_loci_1=new int[cor_pro_1_length];
                        for (int cor_pro_1_1_nozero_locus=0;cor_pro_1_1_nozero_locus<cor_pro_1_length;cor_pro_1_1_nozero_locus++)
                        {
                            random_pro_loci_1[cor_pro_1_1_nozero_locus]=cor_pro_1_1_nozero_locus+1;
                        }
                        shuffle(random_pro_loci_1,random_pro_loci_1+cor_pro_1_length,RNG);
                        //                    gsl_ran_shuffle(r7, random_pro_loci_1, cor_pro_1_length, sizeof (int));
                        Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                        cell_temp(1,1)=(float)cor_pro_1(1,random_pro_loci_1[0]);
                        cell_temp(1,5)=(float)cor_pro_1(2,random_pro_loci_1[0]);
                        cell_temp(1,14)=1;
                        cell_temp(1,15)=cell_temp(1,15);
                        cell_array(i,1)=(float)cor_pro_1(1,random_pro_loci_1[1]);
                        cell_array(i,5)=(float)cor_pro_1(2,random_pro_loci_1[1]);
                        if (cell_temp(1,1)!=0 && cell_temp(1,5)!=0)
                        {
                            Visual_range((int)cell_temp(1,1),(int)cell_temp(1,5),1)=1;
                            Visual_range((int)cell_temp(1,1),(int)cell_temp(1,5),2)=(int)cell_temp(1,15);
                            Visual_range((int)cell_temp(1,1),(int)cell_temp(1,5),3)=cellstage;
                            Visual_range((int)cell_temp(1,1),(int)cell_temp(1,5),4)=cell_label;
                            Visual_range((int)cell_array(i,1),(int)cell_array(i,5),1)=1;
                            Visual_range((int)cell_array(i,1),(int)cell_array(i,5),2)=(int)cell_array(i,15);
                            Visual_range((int)cell_array(i,1),(int)cell_array(i,5),3)=cellstage;
                            Visual_range((int)cell_array(i,1),(int)cell_array(i,5),4)=cell_label_1;
                        }
                        cell_array(i,14)=1;
                        cell_array(i,2)=0;
                        cell_array(i,3)=0;
                        cell_array(i,4)=0;
                        cell_array(i,6)=0;
                        cell_array(i,7)=0;
                        cell_array(i,8)=0;
                        delete[] random_pro_loci_1;
                        random_pro_loci_1 = NULL;
                    }
                }
                else
                {
                    cell_array(i,20)=0;
                    cell_temp(1,20)=0;
                    cell_temp(1,16)=0;
                    cell_array(i,16)=0;
                    int cell_label_1=Visual_range(x,y,4);
                    int cellstage=Visual_range(x,y,3);
                    cell_label=cell_label+1;
                    Array<int, 2> pro_loci_small_1(2,16,FortranArray<2>());
                    pro_loci_small_1=0;
                    int ddddd=1;
                    for (int loci22=1;loci22<=16;loci22++)
                    {
                        if (Visual_range(cor_big_1_change_shape(1,loci22),cor_big_1_change_shape(2,loci22),1)==0)
                        {
                            pro_loci_small_1(1,ddddd)=cor_big_1_change_shape(1,loci22);
                            pro_loci_small_1(2,ddddd)=cor_big_1_change_shape(2,loci22);
                            ddddd=ddddd+1;
                        }
                    }
                    int pro_loci_small_1_nozero_length=0;
                    for (int pro_loci_small_1_locus=1;pro_loci_small_1_locus<=16;pro_loci_small_1_locus++)
                    {
                        if (pro_loci_small_1(1,pro_loci_small_1_locus)!=0)
                        {
                            pro_loci_small_1_nozero_length=pro_loci_small_1_nozero_length+1;
                        }
                    }
                    Array<int, 2> pro_loci_small(2,pro_loci_small_1_nozero_length+4,FortranArray<2>());
                    pro_loci_small=0;
                    int nzl=1;
                    for (int aa=1;aa<=16;aa++ )
                    {
                        if (pro_loci_small_1(1,aa)!=0)
                        {
                            pro_loci_small(all,nzl)=pro_loci_small_1(all,aa);
                            nzl=nzl+1;
                        }
                    }
                    pro_loci_small(1,pro_loci_small_1_nozero_length+1)=x;
                    pro_loci_small(2,pro_loci_small_1_nozero_length+1)=y;
                    pro_loci_small(1,pro_loci_small_1_nozero_length+2)=x;
                    pro_loci_small(2,pro_loci_small_1_nozero_length+2)=y+1;
                    pro_loci_small(1,pro_loci_small_1_nozero_length+3)=x+1;
                    pro_loci_small(2,pro_loci_small_1_nozero_length+3)=y+1;
                    pro_loci_small(1,pro_loci_small_1_nozero_length+4)=x+1;
                    pro_loci_small(2,pro_loci_small_1_nozero_length+4)=y;
                    int pro_loci_small_length=pro_loci_small.columns();
                    int *random_pro_loci_1=new int[pro_loci_small_length];
                    for (int pro_loci_small_nozero_locus=0;pro_loci_small_nozero_locus<pro_loci_small_length;pro_loci_small_nozero_locus++)
                    {
                        random_pro_loci_1[pro_loci_small_nozero_locus]=pro_loci_small_nozero_locus+1;
                    }
                    shuffle(random_pro_loci_1,random_pro_loci_1+pro_loci_small_length,RNG);
                    //                gsl_ran_shuffle(r7, random_pro_loci_1, pro_loci_small_length, sizeof (int));
                    Visual_range(Range(x,x+1),Range(y,y+1),all)=0;
                    cell_temp(1,1)=(float)pro_loci_small(1,random_pro_loci_1[0]);
                    cell_temp(1,5)=(float)pro_loci_small(2,random_pro_loci_1[0]);
                    cell_temp(1,14)=1;
                    cell_temp(1,15)=cell_temp(1,15);
                    cell_array(i,1)=(float)pro_loci_small(1,random_pro_loci_1[1]);
                    cell_array(i,5)=(float)pro_loci_small(2,random_pro_loci_1[1]);
                    if (cell_temp(1,1)!=0 && cell_temp(1,5)!=0)
                    {
                        Visual_range((int)cell_temp(1,1),(int)cell_temp(1,5),1)=1;
                        Visual_range((int)cell_temp(1,1),(int)cell_temp(1,5),2)=(int)cell_temp(1,15);
                        Visual_range((int)cell_temp(1,1),(int)cell_temp(1,5),3)=cellstage;
                        Visual_range((int)cell_temp(1,1),(int)cell_temp(1,5),4)=cell_label;
                        Visual_range((int)cell_array(i,1),(int)cell_array(i,5),1)=1;
                        Visual_range((int)cell_array(i,1),(int)cell_array(i,5),2)=(int)cell_array(i,15);
                        Visual_range((int)cell_array(i,1),(int)cell_array(i,5),3)=cellstage;
                        Visual_range((int)cell_array(i,1),(int)cell_array(i,5),4)=cell_label_1;
                    }
                    cell_array(i,14)=1;
                    cell_array(i,2)=0;
                    cell_array(i,3)=0;
                    cell_array(i,4)=0;
                    cell_array(i,6)=0;
                    cell_array(i,7)=0;
                    cell_array(i,8)=0;
                    delete[] random_pro_loci_1;
                    random_pro_loci_1 = NULL;
                }
                delete[] pro_loci_new;
                pro_loci_new = NULL;
            }
            break;
        }
        case 1:
        {
            float growth_rate_inherent=cell_array(i,10);
            float X1=growth_rate_inherent*(1-0.05);
            float X2=growth_rate_inherent*(1+0.05);
            cell_array(i,10)=(X2-X1)*gsl_rng_uniform(r7)+X1;
            
            r_label=r_label+1;
            K_label=K_label+1;
            if(cell_array(i,9)==1)
            {
                cell_array(i,15)=r_label;
            }
            else if (cell_array(i,9)==1)
            {
                cell_array(i,15)=K_label;
            }
            cell_type_transform(cell_temp, beta_distribution_alpha_for_normal_migration,beta_distribution_beta_for_normal_migration,migration_rate_K_mean,uniup_K,unilow_K, sigmahatK, muhatK,K_label,i,sub_visual,Visual_range,cell_array, beta_distribution_alpha, beta_distribution_beta, migration_rate_r_mean, migration_rate_r_mean_quia, beta_distribution_expected_for_normal_migration,r_label, K_formation_rate);
            
            cell_index=cell_index+1;
            generation=generation+1;
            cell_array(i,30)=cell_array(i,29);
            cell_array(i,29)=cell_index;
            
            cell_index=cell_index+1;
            cell_temp(1,29)=cell_index;
            cell_temp(1,30)=cell_array(i,30);
            
            
            Array<int, 2> cor_temp_2(1,8,FortranArray<2>());
            cor_temp_2=0;
            int cor_temp_length=1;
            for (int s=1;s<=8;s++)
            {
                int x1=cor_small_1(1,s);
                int y1=cor_small_1(2,s);
                if (Visual_range(x1,y1,1)==0)
                {
                    cor_temp_2(1,cor_temp_length)=s;
                    cor_temp_length=cor_temp_length+1;
                }
            }
            int *cor_temp_3=new int[cor_temp_length-1];
            int ln1=0;
            for (int length_cor=0;length_cor<8;length_cor++)
            {
                if(cor_temp_2(1,length_cor)!=0)
                {
                    cor_temp_3[ln1]=cor_temp_2(1,length_cor);
                    ln1=ln1+1;
                }
            }
            int loci_number=ln1;
            int *cor_temp=new int[loci_number];
            for (int ln=0;ln<loci_number;ln++)
            {
                //cor_temp[loci_number]=ln;
                cor_temp[ln]=ln;
            }
            if (loci_number!=0)
            {
                cell_label=cell_label+1;
                int cellstage=Visual_range(x,y,3);
                shuffle(cor_temp, cor_temp+loci_number,RNG);
                //            gsl_ran_shuffle(r7, cor_temp, loci_number, sizeof (int));
                int loci=cor_temp_3[cor_temp[0]];
                cell_temp(1,1)=cor_small_1(1,loci);
                cell_temp(1,5)=cor_small_1(2,loci);
                int cell_index=cell_temp(1,15);
                Visual_range((int)cell_temp(1,1),(int)cell_temp(1,5),1)=1;
                Visual_range((int)cell_temp(1,1),(int)cell_temp(1,5),2)=cell_index;
                Visual_range((int)cell_temp(1,1),(int)cell_temp(1,5),3)=cellstage;
                Visual_range((int)cell_temp(1,1),(int)cell_temp(1,5),4)=cell_label;
                cell_temp(1,14)=1;
                cell_array(i,20)=0;
                cell_temp(1,20)=0;
                cell_temp(1,16)=0;
                cell_array(i,16)=0;
                
                cell_array(i,15)=r_label;
                cell_array(i,13)=cell_array(i,13)+1;
                
                cell_temp(1,Col)=cell_array(i,13);
                cell_array(i,Col)=cell_array(i,13);
            }
            else
            {
                if (cell_array(i,9)==1)
                {
                    cell_array(i,22)=0;
                    Visual_range(cell_array(i,1),cell_array(i,5),all)=0;
                    cell_temp(1,22)=0;
                }
                else
                {
                    if (utralsmall==1)
                    {
                        cell_temp(1,1)=cell_array(i,1);
                        cell_temp(1,5)=cell_array(i,5);
                        cell_array(i,14)=2;
                        cell_temp(1,14)=2;
                        Visual_range((int)cell_temp(1,1),(int)cell_temp(1,5),3)=2;
                        Visual_range((int)cell_temp(1,1),(int)cell_temp(1,5),4)=cell_label;
                        cell_array(i,20)=0;
                        cell_temp(1,20)=0;
                        cell_temp(1,16)=0;
                        cell_array(i,16)=0;
                        cell_array(i,13)=cell_array(i,13)+1;
                        cell_temp(1,Col)=cell_array(i,13);
                        cell_array(i,Col)=cell_array(i,13);
                        
                    }
                }
            }
            delete[] cor_temp;
            cor_temp = NULL;
            delete[] cor_temp_3;
            cor_temp_3 = NULL;
            break;
        }
        case 2:
        {
            int cor_temp_length=0;
            for (int s=1;s<=8;s++)
            {
                int x1=cor_small_1(1,s);
                int y1=cor_small_1(2,s);
                if (Visual_range(x1,y1,1)==0)
                {
                    cor_temp_length=cor_temp_length+1;
                }
            }
            if (cor_temp_length==0)
            {
                cell_array(i,22)=0;
            }
            else
            {
                cell_array(i,16)=cell_array(i,16)+deltah;
            }
            break;
        }
    }
    
    int celltypes=(int)cell_array(i,9);
    switch (celltypes)
    {
        case 1:
        {
            if (cell_array(i,10) > max_growth_rate_r)
            {
                cell_array(i,10) = max_growth_rate_r;
            }
            break;
        }
        case 2:
        {
            if(cell_array(i,10) > max_growth_rate_K)
            {
                cell_array(i,10) = max_growth_rate_K;
            }
            break;
        }
    }
    
    int celltypes_temp=(int)cell_temp(i,9);
    switch (celltypes_temp)
    {
        case 1:
        {
            if(cell_temp(1,10) > max_growth_rate_r)
            {
                cell_temp(1,10) = max_growth_rate_r;
            }
            break;
        }
        case 2:
        {
            if(cell_temp(1,10) > max_growth_rate_K)
            {
                cell_temp(1,10) = max_growth_rate_K;
            }
            break;
        }
    }
    
    if (cell_temp(1,1)!=0 && cell_temp(1,5)!=0)
    {
        cell_trace_temp.resize(2,150);
        cell_trace_temp(1,Range(2,4))=cell_array(i,Range(29,Col));
        cell_trace_temp(2,Range(2,4))=cell_temp(1,Range(29,Col));
        cell_trace_temp(1,1)=cell_temp(1,15);
        cell_trace_temp(2,1)=cell_temp(1,15);
        
        int current_size_trace=cell_trace.rows();
        int GN=cell_trace_temp(1,4);
        for (int rows=current_size_trace;rows>=1;rows--)
        {
            if(cell_trace_temp(1,3) == cell_trace(rows,2))
            {
                cell_trace_temp(1,Range(5,150))=cell_trace(rows,Range(5,150));
                cell_trace_temp(2,Range(5,150))=cell_trace(rows,Range(5,150));
                break;
            }
        }
        
        cell_trace_temp(1,GN+5)=1;
        cell_trace_temp(2,GN+5)=2;
        
        
        cell_trace.resizeAndPreserve(current_size_trace+2,150);
        cell_trace(current_size_trace+1,all)=cell_trace_temp(1,all);
        cell_trace(current_size_trace+2,all)=cell_trace_temp(2,all);
        
        
        
        int current_size=cell_array.rows();
        
        cell_array.resizeAndPreserve(current_size+1,Col);
        cell_array(current_size+1,all)=cell_temp(1,all);
        
    }
    gsl_rng_free(r7);
}

#endif /* free_living_division_hpp */
