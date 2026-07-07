//
//  division.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef division_hpp
#define division_hpp

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
#include "stateless_rng.hpp"
#include "cell_store.hpp"
#include "visual_range.hpp"
#include "int_matrix.hpp"
#include <chrono>

using std::chrono::high_resolution_clock;
using namespace blitz;

namespace division_detail
{
inline void append_cell(CellStore &cells, const CellRowBuffer &cell_temp, int Col)
{
    cell_store_append_row_from_array(cells, cell_temp, 1, Col);
}
}

template <typename CellArray>
inline void division(int i, double max_growth_rate_r, double max_growth_rate_K, CellArray &cell_array, VisualRange &Visual_range, CellRowBuffer &cell_temp,int &cell_label, double &deltah,int utralsmall,int Col, long rng_time_step)
{
    int row = i - 1;
    auto &x1_values = cell_array.x1();
    auto &x2_values = cell_array.x2();
    auto &x3_values = cell_array.x3();
    auto &x4_values = cell_array.x4();
    auto &y1_values = cell_array.y1();
    auto &y2_values = cell_array.y2();
    auto &y3_values = cell_array.y3();
    auto &y4_values = cell_array.y4();
    auto &types = cell_array.type();
    auto &growth_rates = cell_array.growth_rate();
    auto &density_growth_rates = cell_array.density_growth_rate();
    auto &migration_rate_bases = cell_array.migration_rate_base();
    auto &random_labels = cell_array.random_label();
    auto &stages = cell_array.stage();
    auto &ids = cell_array.id();
    auto &division_elapsed = cell_array.division_elapsed();
    auto &death_times = cell_array.death_time();
    auto &death_elapsed = cell_array.death_elapsed();
    auto &migration_elapsed = cell_array.migration_elapsed();
    auto &migration_intervals = cell_array.migration_interval();
    auto &viability = cell_array.viability();
    auto &migration_follow_flags = cell_array.migration_follow_flag();
    auto &migration_active = cell_array.migration_active();
    auto &migration_duration = cell_array.migration_duration();
    auto &migration_passed = cell_array.migration_passed();
    auto &migration_rates = cell_array.migration_rate();
    long cell_rng_id = (long)ids[row];
    if (cell_rng_id == 0)
    {
        cell_rng_id = i;
    }
    long rng_event = 100;
    int pro_loci[8]={0};
    int pro_loci1[4]={0};
    int pro_loci2[4]={0};
    int pro_loci1_new[4]={0};
    int pro_loci2_new[4]={0};
    IntMatrix cor_big_1(2, 16);
    IntMatrix cor_big_1_change_shape(2, 16);
    IntMatrix cor_small_1(2, 8);
    cell_temp.resize(1, Col);
    cell_temp=0;
    int x=(int)x1_values[row];
    int y=(int)y1_values[row];
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
    cor_big_1_change_shape.set_col(1, B, b);
    cor_big_1_change_shape.set_col(2, B, c);
    cor_big_1_change_shape.set_col(3, B, d);
    cor_big_1_change_shape.set_col(4, B, e);
    cor_big_1_change_shape.set_col(5, C, e);
    cor_big_1_change_shape.set_col(6, D, e);
    cor_big_1_change_shape.set_col(7, E, e);
    cor_big_1_change_shape.set_col(8, E, d);
    cor_big_1_change_shape.set_col(9, E, c);
    cor_big_1_change_shape.set_col(10, E, b);
    cor_big_1_change_shape.set_col(11, D, b);
    cor_big_1_change_shape.set_col(12, C, b);
    cor_big_1_change_shape.set_col(13, C, c);
    cor_big_1_change_shape.set_col(14, C, d);
    cor_big_1_change_shape.set_col(15, D, d);
    cor_big_1_change_shape.set_col(16, D, c);
    cor_small_1.set_col(1, B, b);
    cor_small_1.set_col(2, B, c);
    cor_small_1.set_col(3, B, d);
    cor_small_1.set_col(4, C, d);
    cor_small_1.set_col(5, D, d);
    cor_small_1.set_col(6, D, c);
    cor_small_1.set_col(7, D, b);
    cor_small_1.set_col(8, C, b);
    int cor_temp1[16]={0};
    
    if (stages[row]==0)
    {
        double growth_rate_inherent=growth_rates[row];
        double X1=growth_rate_inherent*(1-0.05);
        double X2=growth_rate_inherent*(1+0.05);
        growth_rates[row]=(X2-X1)*stateless_uniform(cell_rng_id, rng_time_step, rng_event++)+X1;
        cell_temp(1,9)=types[row];
        cell_temp(1,10)=growth_rates[row];
        cell_temp(1,11)=density_growth_rates[row];
        cell_temp(1,12)=migration_rate_bases[row];
        cell_temp(1,13)=random_labels[row];
        cell_temp(1,15)=ids[row];
        cell_temp(1,18)=0;
        cell_temp(1,19)=death_elapsed[row];
        cell_temp(1,21)=migration_intervals[row];
        cell_temp(1,22)=viability[row];
        cell_temp(1,23)=0;
        cell_temp(1,24)=0;
        division_elapsed[row]=0;
        cell_temp(1,16)=0;
        cell_temp(1,25)=migration_active[row];
        cell_temp(1,26)=migration_duration[row];
        cell_temp(1,27)=migration_passed[row];
        cell_temp(1,28)=migration_rates[row];
        int cor_temp_length=0;
        for (int s=1; s<=16; s++)
        {
            int x1=cor_big_1(1,s);
            int y1=cor_big_1(2,s);
            if (Visual_range.occupied(x1,y1)==0 && Visual_range.occupied(x1,y1+1)==0 && Visual_range.occupied(x1+1,y1)==0 && Visual_range.occupied(x1+1,y1+1)==0)
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
            stateless_shuffle(cor_big_temp_1,cor_big_temp_1+size_big, cell_rng_id, rng_time_step, rng_event++);
            int loci=cor_big_temp_1[0];
            cell_temp(1,1)=cor_big_1(1,loci);
            cell_temp(1,2)=cor_big_1(1,loci);
            cell_temp(1,3)=cor_big_1(1,loci)+1;
            cell_temp(1,4)=cor_big_1(1,loci)+1;
            cell_temp(1,5)=cor_big_1(2,loci);
            cell_temp(1,6)=cor_big_1(2,loci)+1;
            cell_temp(1,7)=cor_big_1(2,loci)+1;
            cell_temp(1,8)=cor_big_1(2,loci);
            cell_temp(1,14)=stages[row];
            cell_temp(1,22)=viability[row];
            cell_temp(1,24)=migration_follow_flags[row];
            migration_elapsed[row]=0;
            cell_temp(1,20)=0;
            cell_temp(1,16)=0;
            division_elapsed[row]=0;
            long cell_index= (int)ids[row];
            Visual_range.occupied((int)cell_temp(1,1),(int)cell_temp(1,5))=1;
            Visual_range.occupied((int)cell_temp(1,2),(int)cell_temp(1,6))=1;
            Visual_range.occupied((int)cell_temp(1,3),(int)cell_temp(1,7))=1;
            Visual_range.occupied((int)cell_temp(1,4),(int)cell_temp(1,8))=1;
            Visual_range.density_label((int)cell_temp(1,1),(int)cell_temp(1,5))=cell_index;
            Visual_range.density_label((int)cell_temp(1,2),(int)cell_temp(1,6))=cell_index;
            Visual_range.density_label((int)cell_temp(1,3),(int)cell_temp(1,7))=cell_index;
            Visual_range.density_label((int)cell_temp(1,4),(int)cell_temp(1,8))=cell_index;
            Visual_range.stage((int)cell_temp(1,1),(int)cell_temp(1,5))=(int)stages[row];
            Visual_range.stage((int)cell_temp(1,2),(int)cell_temp(1,6))=(int)stages[row];
            Visual_range.stage((int)cell_temp(1,3),(int)cell_temp(1,7))=(int)stages[row];
            Visual_range.stage((int)cell_temp(1,4),(int)cell_temp(1,8))=(int)stages[row];
            Visual_range.cell_label((int)cell_temp(1,1),(int)cell_temp(1,5))=cell_label;
            Visual_range.cell_label((int)cell_temp(1,2),(int)cell_temp(1,6))=cell_label;
            Visual_range.cell_label((int)cell_temp(1,3),(int)cell_temp(1,7))=cell_label;
            Visual_range.cell_label((int)cell_temp(1,4),(int)cell_temp(1,8))=cell_label;
        }
        else
        {
            if  (Visual_range.occupied(cor_big_1_change_shape(1,1),cor_big_1_change_shape(2,1))==0 && Visual_range.occupied(cor_big_1_change_shape(1,2),cor_big_1_change_shape(2,2))==0 && Visual_range.occupied(cor_big_1_change_shape(1,12),cor_big_1_change_shape(2,12))==0)
            {
                pro_loci[0]=1;
                pro_loci1[0]=1;
            }
            if (Visual_range.occupied(cor_big_1_change_shape(1,2),cor_big_1_change_shape(2,2))==0 && Visual_range.occupied(cor_big_1_change_shape(1,3),cor_big_1_change_shape(2,3))==0)
            {
                pro_loci[1]=2;
                pro_loci2[0]=2;
            }
            if (Visual_range.occupied(cor_big_1_change_shape(1,3),cor_big_1_change_shape(2,3))==0 && Visual_range.occupied(cor_big_1_change_shape(1,4),cor_big_1_change_shape(2,4))==0 && Visual_range.occupied(cor_big_1_change_shape(1,5),cor_big_1_change_shape(2,5))==0)
            {
                pro_loci[2]=3;
                pro_loci1[1]=3;
            }
            if (Visual_range.occupied(cor_big_1_change_shape(1,5),cor_big_1_change_shape(2,5))==0 && Visual_range.occupied(cor_big_1_change_shape(1,6),cor_big_1_change_shape(2,6))==0)
            {
                pro_loci[3]=4;
                pro_loci2[1]=4;
            }
            if (Visual_range.occupied(cor_big_1_change_shape(1,6),cor_big_1_change_shape(2,6))==0 && Visual_range.occupied(cor_big_1_change_shape(1,7),cor_big_1_change_shape(2,7))==0 && Visual_range.occupied(cor_big_1_change_shape(1,8),cor_big_1_change_shape(2,8))==0)
            {
                pro_loci[4]=5;
                pro_loci1[2]=5;
            }
            if (Visual_range.occupied(cor_big_1_change_shape(1,8),cor_big_1_change_shape(2,8))==0 && Visual_range.occupied(cor_big_1_change_shape(1,9),cor_big_1_change_shape(2,9))==0)
            {
                pro_loci[5]=6;
                pro_loci2[2]=6;
            }
            if (Visual_range.occupied(cor_big_1_change_shape(1,9),cor_big_1_change_shape(2,9))==0 && Visual_range.occupied(cor_big_1_change_shape(1,10),cor_big_1_change_shape(2,10))==0 && Visual_range.occupied(cor_big_1_change_shape(1,11),cor_big_1_change_shape(2,11))==0)
            {
                pro_loci[6]=7;
                pro_loci1[3]=7;
            }
            if (Visual_range.occupied(cor_big_1_change_shape(1,11),cor_big_1_change_shape(2,11))==0 && Visual_range.occupied(cor_big_1_change_shape(1,12),cor_big_1_change_shape(2,12))==0)
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
                migration_elapsed[row]=0;
                cell_temp(1,20)=0;
                cell_temp(1,16)=0;
                division_elapsed[row]=0;
                long cell_label_1=Visual_range.cell_label(x,y);
                long cellstage=Visual_range.stage(x,y);
                cell_label=cell_label+1;
                if (pro_loci_number1 > 1 && pro_loci_number2 >1)
                {
                    int random_pro_loci[2]={1,2};
                    stateless_shuffle(random_pro_loci,random_pro_loci+2, cell_rng_id, rng_time_step, rng_event++);
                    if (random_pro_loci[0]==1)
                    {
                        int sizep1=num1;
                        int random_pro_loci1[sizep1];
                        for(int aa=0;aa<sizep1;aa++)
                        {
                            random_pro_loci1[aa]=pro_loci1_new[aa];
                        }
                        stateless_shuffle(random_pro_loci1,random_pro_loci1+sizep1, cell_rng_id, rng_time_step, rng_event++);
                        int proloci1=random_pro_loci1[0];
                        int proloci2=random_pro_loci1[1];
                        if (proloci1==1)
                        {
                            Visual_range.clear_square(x,y);
                            for (int cor_cell=1; cor_cell<=4;cor_cell++)
                            {
                                int cor_cell_y=cor_cell+4;
                                cell_array.column(cor_cell)[row]=cell_array.column(cor_cell)[row]-1;
                                cell_array.column(cor_cell_y)[row]=cell_array.column(cor_cell_y)[row]-1;
                            }
                            Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                            Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                            Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                            Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                        }
                        else if (proloci1==3)
                        {
                            Visual_range.clear_square(x,y);
                            for (int cor_cell=1; cor_cell<=4;cor_cell++)
                            {
                                int cor_cell_y=cor_cell+4;
                                cell_array.column(cor_cell)[row]=cell_array.column(cor_cell)[row]-1;
                                cell_array.column(cor_cell_y)[row]=cell_array.column(cor_cell_y)[row]+1;
                            }
                            Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                            Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                            Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                            Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                        }
                        else if (proloci1==5)
                        {
                            Visual_range.clear_square(x,y);
                            for (int cor_cell=1; cor_cell<=4;cor_cell++)
                            {
                                int cor_cell_y=cor_cell+4;
                                cell_array.column(cor_cell)[row]=cell_array.column(cor_cell)[row]+1;
                                cell_array.column(cor_cell_y)[row]=cell_array.column(cor_cell_y)[row]+1;
                            }
                            Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                            Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                            Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                            Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                        }
                        else if (proloci1==7)
                        {
                            Visual_range.clear_square(x,y);
                            for (int cor_cell=1; cor_cell<=4;cor_cell++)
                            {
                                int cor_cell_y=cor_cell+4;
                                cell_array.column(cor_cell)[row]=cell_array.column(cor_cell)[row]+1;
                                cell_array.column(cor_cell_y)[row]=cell_array.column(cor_cell_y)[row]-1;
                            }
                            Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                            Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                            Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                            Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                        }
                        if (proloci2==1)
                        {
                            Visual_range.clear_square(x,y);
                            cell_temp(1,1)=x-1;
                            cell_temp(1,5)=y-1;
                            cell_temp(1,2)=cell_temp(1,1);
                            cell_temp(1,3)=cell_temp(1,1)+1;
                            cell_temp(1,4)=cell_temp(1,1)+1;
                            cell_temp(1,6)=cell_temp(1,5)+1;
                            cell_temp(1,7)=cell_temp(1,5)+1;
                            cell_temp(1,8)=cell_temp(1,5);
                            cell_temp(1,14)=stages[row];
                            cell_temp(1,22)=viability[row];
                            cell_temp(1,24)=migration_follow_flags[row];
                            Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                            Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                            Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                            Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                        }
                        else if (proloci2==3)
                        {
                            Visual_range.clear_square(x,y);
                            cell_temp(1,1)=x-1;
                            cell_temp(1,5)=y+1;
                            cell_temp(1,2)=cell_temp(1,1);
                            cell_temp(1,3)=cell_temp(1,1)+1;
                            cell_temp(1,4)=cell_temp(1,1)+1;
                            cell_temp(1,6)=cell_temp(1,5)+1;
                            cell_temp(1,7)=cell_temp(1,5)+1;
                            cell_temp(1,8)=cell_temp(1,5);
                            cell_temp(1,14)=stages[row];
                            cell_temp(1,22)=viability[row];
                            cell_temp(1,24)=migration_follow_flags[row];
                            Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                            Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                            Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                            Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                        }
                        else if (proloci2==5)
                        {
                            Visual_range.clear_square(x,y);
                            cell_temp(1,1)=x+1;
                            cell_temp(1,5)=y+1;
                            cell_temp(1,2)=cell_temp(1,1);
                            cell_temp(1,3)=cell_temp(1,1)+1;
                            cell_temp(1,4)=cell_temp(1,1)+1;
                            cell_temp(1,6)=cell_temp(1,5)+1;
                            cell_temp(1,7)=cell_temp(1,5)+1;
                            cell_temp(1,8)=cell_temp(1,5);
                            cell_temp(1,14)=stages[row];
                            cell_temp(1,22)=viability[row];
                            cell_temp(1,24)=migration_follow_flags[row];
                            Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                            Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                            Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                            Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                        }
                        else if (proloci2==7)
                        {
                            Visual_range.clear_square(x,y);
                            cell_temp(1,1)=x+1;
                            cell_temp(1,5)=y-1;
                            cell_temp(1,2)=cell_temp(1,1);
                            cell_temp(1,3)=cell_temp(1,1)+1;
                            cell_temp(1,4)=cell_temp(1,1)+1;
                            cell_temp(1,6)=cell_temp(1,5)+1;
                            cell_temp(1,7)=cell_temp(1,5)+1;
                            cell_temp(1,8)=cell_temp(1,5);
                            cell_temp(1,14)=stages[row];
                            cell_temp(1,22)=viability[row];
                            cell_temp(1,24)=migration_follow_flags[row];
                            Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                            Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                            Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                            Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                        }
                    }
                    else if (random_pro_loci[0]==2)
                    {
                        int sizep2=num2;
                        int random_pro_loci2[sizep2];
                        for(int aa=0;aa<sizep2;aa++)
                        {
                            random_pro_loci2[aa]=pro_loci2_new[aa];
                        }
                        stateless_shuffle(random_pro_loci2,random_pro_loci2+sizep2, cell_rng_id, rng_time_step, rng_event++);
                        int proloci1=random_pro_loci2[0];
                        int proloci2=random_pro_loci2[1];
                        if (proloci1==2)
                        {
                            Visual_range.clear_square(x,y);
                            for (int cor_cell=1; cor_cell<=4;cor_cell++)
                            {
                                cell_array.column(cor_cell)[row]=cell_array.column(cor_cell)[row]-1;
                            }
                            Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                            Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                            Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                            Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                        }
                        else if (proloci1==4)
                        {
                            Visual_range.clear_square(x,y);
                            for (int cor_cell_y=5; cor_cell_y<=8;cor_cell_y++)
                            {
                                cell_array.column(cor_cell_y)[row]=cell_array.column(cor_cell_y)[row]+1;
                            }
                            Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                            Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                            Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                            Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                        }
                        else if (proloci1==6)
                        {
                            Visual_range.clear_square(x,y);
                            for (int cor_cell=1; cor_cell<=4;cor_cell++)
                            {
                                cell_array.column(cor_cell)[row]=cell_array.column(cor_cell)[row]+1;
                            }
                            Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                            Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                            Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                            Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                        }
                        else if (proloci1==8)
                        {
                            Visual_range.clear_square(x,y);
                            for (int cor_cell_y=5; cor_cell_y<=8;cor_cell_y++)
                            {
                                cell_array.column(cor_cell_y)[row]=cell_array.column(cor_cell_y)[row]-1;
                            }
                            Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                            Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                            Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                            Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                        }
                        if (proloci2==2)
                        {
                            Visual_range.clear_square(x,y);
                            cell_temp(1,1)=x-1;
                            cell_temp(1,5)=y;
                            cell_temp(1,2)=cell_temp(1,1);
                            cell_temp(1,3)=cell_temp(1,1)+1;
                            cell_temp(1,4)=cell_temp(1,1)+1;
                            cell_temp(1,6)=cell_temp(1,5)+1;
                            cell_temp(1,7)=cell_temp(1,5)+1;
                            cell_temp(1,8)=cell_temp(1,5);
                            cell_temp(1,14)=stages[row];
                            cell_temp(1,22)=viability[row];
                            cell_temp(1,24)=migration_follow_flags[row];
                            Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                            Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                            Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                            Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                        }
                        else if (proloci2==4)
                        {
                            Visual_range.clear_square(x,y);
                            cell_temp(1,1)=x;
                            cell_temp(1,5)=y+1;
                            cell_temp(1,2)=cell_temp(1,1);
                            cell_temp(1,3)=cell_temp(1,1)+1;
                            cell_temp(1,4)=cell_temp(1,1)+1;
                            cell_temp(1,6)=cell_temp(1,5)+1;
                            cell_temp(1,7)=cell_temp(1,5)+1;
                            cell_temp(1,8)=cell_temp(1,5);
                            cell_temp(1,14)=stages[row];
                            cell_temp(1,22)=viability[row];
                            cell_temp(1,24)=migration_follow_flags[row];
                            Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                            Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                            Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                            Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                        }
                        else if (proloci2==6)
                        {
                            Visual_range.clear_square(x,y);
                            cell_temp(1,1)=x+1;
                            cell_temp(1,5)=y;
                            cell_temp(1,2)=cell_temp(1,1);
                            cell_temp(1,3)=cell_temp(1,1)+1;
                            cell_temp(1,4)=cell_temp(1,1)+1;
                            cell_temp(1,6)=cell_temp(1,5)+1;
                            cell_temp(1,7)=cell_temp(1,5)+1;
                            cell_temp(1,8)=cell_temp(1,5);
                            cell_temp(1,14)=stages[row];
                            cell_temp(1,22)=viability[row];
                            cell_temp(1,24)=migration_follow_flags[row];
                            Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                            Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                            Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                            Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                        }
                        else if (proloci2==8)
                        {
                            Visual_range.clear_square(x,y);
                            cell_temp(1,1)=x;
                            cell_temp(1,5)=y-1;
                            cell_temp(1,2)=cell_temp(1,1);
                            cell_temp(1,3)=cell_temp(1,1)+1;
                            cell_temp(1,4)=cell_temp(1,1)+1;
                            cell_temp(1,6)=cell_temp(1,5)+1;
                            cell_temp(1,7)=cell_temp(1,5)+1;
                            cell_temp(1,8)=cell_temp(1,5);
                            cell_temp(1,14)=stages[row];
                            cell_temp(1,22)=viability[row];
                            cell_temp(1,24)=migration_follow_flags[row];
                            Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                            Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                            Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                            Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
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
                    stateless_shuffle(random_pro_loci1,random_pro_loci1+sizep1, cell_rng_id, rng_time_step, rng_event++);
                    int proloci1=random_pro_loci1[0];
                    int proloci2=random_pro_loci1[1];
                    if (proloci1==1)
                    {
                        Visual_range.clear_square(x,y);
                        for (int cor_cell=1; cor_cell<=4;cor_cell++)
                        {
                            int cor_cell_y=cor_cell+4;
                            cell_array.column(cor_cell)[row]=cell_array.column(cor_cell)[row]-1;
                            cell_array.column(cor_cell_y)[row]=cell_array.column(cor_cell_y)[row]-1;
                        }
                        Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                        Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                        Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                        Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                    }
                    else if (proloci1==3)
                    {
                        Visual_range.clear_square(x,y);
                        for (int cor_cell=1; cor_cell<=4;cor_cell++)
                        {
                            int cor_cell_y=cor_cell+4;
                            cell_array.column(cor_cell)[row]=cell_array.column(cor_cell)[row]-1;
                            cell_array.column(cor_cell_y)[row]=cell_array.column(cor_cell_y)[row]+1;
                        }
                        Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                        Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                        Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                        Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                    }
                    else if (proloci1==5)
                    {
                        Visual_range.clear_square(x,y);
                        for (int cor_cell=1; cor_cell<=4;cor_cell++)
                        {
                            int cor_cell_y=cor_cell+4;
                            cell_array.column(cor_cell)[row]=cell_array.column(cor_cell)[row]+1;
                            cell_array.column(cor_cell_y)[row]=cell_array.column(cor_cell_y)[row]+1;
                        }
                        Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                        Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                        Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                        Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                    }
                    else if (proloci1==7)
                    {
                        Visual_range.clear_square(x,y);
                        for (int cor_cell=1; cor_cell<=4;cor_cell++)
                        {
                            int cor_cell_y=cor_cell+4;
                            cell_array.column(cor_cell)[row]=cell_array.column(cor_cell)[row]+1;
                            cell_array.column(cor_cell_y)[row]=cell_array.column(cor_cell_y)[row]-1;
                        }
                        Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                        Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                        Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                        Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                    }
                    if (proloci2==1)
                    {
                        Visual_range.clear_square(x,y);
                        cell_temp(1,1)=x-1;
                        cell_temp(1,5)=y-1;
                        cell_temp(1,2)=cell_temp(1,1);
                        cell_temp(1,3)=cell_temp(1,1)+1;
                        cell_temp(1,4)=cell_temp(1,1)+1;
                        cell_temp(1,6)=cell_temp(1,5)+1;
                        cell_temp(1,7)=cell_temp(1,5)+1;
                        cell_temp(1,8)=cell_temp(1,5);
                        cell_temp(1,14)=stages[row];
                        cell_temp(1,22)=viability[row];
                        cell_temp(1,24)=migration_follow_flags[row];
                        Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                        Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                        Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                        Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                    }
                    else if (proloci2==3)
                    {
                        Visual_range.clear_square(x,y);
                        cell_temp(1,1)=x-1;
                        cell_temp(1,5)=y+1;
                        cell_temp(1,2)=cell_temp(1,1);
                        cell_temp(1,3)=cell_temp(1,1)+1;
                        cell_temp(1,4)=cell_temp(1,1)+1;
                        cell_temp(1,6)=cell_temp(1,5)+1;
                        cell_temp(1,7)=cell_temp(1,5)+1;
                        cell_temp(1,8)=cell_temp(1,5);
                        cell_temp(1,14)=stages[row];
                        cell_temp(1,22)=viability[row];
                        cell_temp(1,24)=migration_follow_flags[row];
                        Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                        Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                        Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                        Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                    }
                    else if (proloci2==5)
                    {
                        Visual_range.clear_square(x,y);
                        cell_temp(1,1)=x+1;
                        cell_temp(1,5)=y+1;
                        cell_temp(1,2)=cell_temp(1,1);
                        cell_temp(1,3)=cell_temp(1,1)+1;
                        cell_temp(1,4)=cell_temp(1,1)+1;
                        cell_temp(1,6)=cell_temp(1,5)+1;
                        cell_temp(1,7)=cell_temp(1,5)+1;
                        cell_temp(1,8)=cell_temp(1,5);
                        cell_temp(1,14)=stages[row];
                        cell_temp(1,22)=viability[row];
                        cell_temp(1,24)=migration_follow_flags[row];
                        Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                        Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                        Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                        Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                    }
                    else if (proloci2==7)
                    {
                        Visual_range.clear_square(x,y);
                        cell_temp(1,1)=x+1;
                        cell_temp(1,5)=y-1;
                        cell_temp(1,2)=cell_temp(1,1);
                        cell_temp(1,3)=cell_temp(1,1)+1;
                        cell_temp(1,4)=cell_temp(1,1)+1;
                        cell_temp(1,6)=cell_temp(1,5)+1;
                        cell_temp(1,7)=cell_temp(1,5)+1;
                        cell_temp(1,8)=cell_temp(1,5);
                        cell_temp(1,14)=stages[row];
                        cell_temp(1,22)=viability[row];
                        cell_temp(1,24)=migration_follow_flags[row];
                        Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                        Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                        Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                        Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
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
                    stateless_shuffle(random_pro_loci2,random_pro_loci2+sizep2, cell_rng_id, rng_time_step, rng_event++);
                    int proloci1=random_pro_loci2[0];
                    int proloci2=random_pro_loci2[1];
                    if (proloci1==2)
                    {
                        Visual_range.clear_square(x,y);
                        for (int cor_cell=1; cor_cell<=4;cor_cell++)
                        {
                            cell_array.column(cor_cell)[row]=cell_array.column(cor_cell)[row]-1;
                        }
                        Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                        Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                        Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                        Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                    }
                    else if (proloci1==4)
                    {
                        Visual_range.clear_square(x,y);
                        for (int cor_cell_y=5; cor_cell_y<=8;cor_cell_y++)
                        {
                            cell_array.column(cor_cell_y)[row]=cell_array.column(cor_cell_y)[row]+1;
                        }
                        Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                        Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                        Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                        Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                    }
                    else if (proloci1==6)
                    {
                        Visual_range.clear_square(x,y);
                        for (int cor_cell=1; cor_cell<=4;cor_cell++)
                        {
                            cell_array.column(cor_cell)[row]=cell_array.column(cor_cell)[row]+1;
                        }
                        Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                        Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                        Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                        Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                    }
                    else if (proloci1==8)
                    {
                        Visual_range.clear_square(x,y);
                        for (int cor_cell_y=5; cor_cell_y<=8;cor_cell_y++)
                        {
                            cell_array.column(cor_cell_y)[row]=cell_array.column(cor_cell_y)[row]-1;
                        }
                        Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                        Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                        Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                        Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                    }
                    if (proloci2==2)
                    {
                        Visual_range.clear_square(x,y);
                        cell_temp(1,1)=x-1;
                        cell_temp(1,5)=y;
                        cell_temp(1,2)=cell_temp(1,1);
                        cell_temp(1,3)=cell_temp(1,1)+1;
                        cell_temp(1,4)=cell_temp(1,1)+1;
                        cell_temp(1,6)=cell_temp(1,5)+1;
                        cell_temp(1,7)=cell_temp(1,5)+1;
                        cell_temp(1,8)=cell_temp(1,5);
                        cell_temp(1,14)=stages[row];
                        cell_temp(1,22)=viability[row];
                        cell_temp(1,24)=migration_follow_flags[row];
                        Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                        Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                        Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                        Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                        
                    }
                    else if (proloci2==4)
                    {
                        Visual_range.clear_square(x,y);
                        cell_temp(1,1)=x;
                        cell_temp(1,5)=y+1;
                        cell_temp(1,2)=cell_temp(1,1);
                        cell_temp(1,3)=cell_temp(1,1)+1;
                        cell_temp(1,4)=cell_temp(1,1)+1;
                        cell_temp(1,6)=cell_temp(1,5)+1;
                        cell_temp(1,7)=cell_temp(1,5)+1;
                        cell_temp(1,8)=cell_temp(1,5);
                        cell_temp(1,14)=stages[row];
                        cell_temp(1,22)=viability[row];
                        cell_temp(1,24)=migration_follow_flags[row];
                        Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                        Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                        Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                        Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                    }
                    else if (proloci2==6)
                    {
                        Visual_range.clear_square(x,y);
                        cell_temp(1,1)=x+1;
                        cell_temp(1,5)=y;
                        cell_temp(1,2)=cell_temp(1,1);
                        cell_temp(1,3)=cell_temp(1,1)+1;
                        cell_temp(1,4)=cell_temp(1,1)+1;
                        cell_temp(1,6)=cell_temp(1,5)+1;
                        cell_temp(1,7)=cell_temp(1,5)+1;
                        cell_temp(1,8)=cell_temp(1,5);
                        cell_temp(1,14)=stages[row];
                        cell_temp(1,22)=viability[row];
                        cell_temp(1,24)=migration_follow_flags[row];
                        Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                        Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                        Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                        Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                    }
                    else if (proloci2==8)
                    {
                        Visual_range.clear_square(x,y);
                        cell_temp(1,1)=x;
                        cell_temp(1,5)=y-1;
                        cell_temp(1,2)=cell_temp(1,1);
                        cell_temp(1,3)=cell_temp(1,1)+1;
                        cell_temp(1,4)=cell_temp(1,1)+1;
                        cell_temp(1,6)=cell_temp(1,5)+1;
                        cell_temp(1,7)=cell_temp(1,5)+1;
                        cell_temp(1,8)=cell_temp(1,5);
                        cell_temp(1,14)=stages[row];
                        cell_temp(1,22)=viability[row];
                        cell_temp(1,24)=migration_follow_flags[row];
                        Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                        Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                        Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                        Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                    }
                }
                else if (pro_loci_number1 ==1 && pro_loci_number2 ==1)
                {
                    if (abs(pro_loci_number2-pro_loci_number1)>=2 && pro_loci_number2-pro_loci_number1!=7)
                    {
                        int random_pro_loci[2]={1,2};
                        stateless_shuffle(random_pro_loci,random_pro_loci+2, cell_rng_id, rng_time_step, rng_event++);
                        if (random_pro_loci[0]==1)
                        {
                            int proloci1=pro_loci1_new[0];
                            int proloci2=pro_loci2_new[0];
                            if (proloci1==1)
                            {
                                Visual_range.clear_square(x,y);
                                for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                {
                                    int cor_cell_y=cor_cell+4;
                                    cell_array.column(cor_cell)[row]=cell_array.column(cor_cell)[row]-1;
                                    cell_array.column(cor_cell_y)[row]=cell_array.column(cor_cell_y)[row]-1;
                                }
                                Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                                Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                                Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                                Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                            }
                            else if (proloci1==3)
                            {
                                Visual_range.clear_square(x,y);
                                for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                {
                                    int cor_cell_y=cor_cell+4;
                                    cell_array.column(cor_cell)[row]=cell_array.column(cor_cell)[row]-1;
                                    cell_array.column(cor_cell_y)[row]=cell_array.column(cor_cell_y)[row]+1;
                                }
                                Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                                Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                                Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                                Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                            }
                            else if (proloci1==5)
                            {
                                Visual_range.clear_square(x,y);
                                for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                {
                                    int cor_cell_y=cor_cell+4;
                                    cell_array.column(cor_cell)[row]=cell_array.column(cor_cell)[row]+1;
                                    cell_array.column(cor_cell_y)[row]=cell_array.column(cor_cell_y)[row]+1;
                                }
                                Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                                Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                                Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                                Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                            }
                            else if (proloci1==7)
                            {
                                Visual_range.clear_square(x,y);
                                for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                {
                                    int cor_cell_y=cor_cell+4;
                                    cell_array.column(cor_cell)[row]=cell_array.column(cor_cell)[row]+1;
                                    cell_array.column(cor_cell_y)[row]=cell_array.column(cor_cell_y)[row]-1;
                                }
                                Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                                Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                                Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                                Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                            }
                            if (proloci2==2)
                            {
                                Visual_range.clear_square(x,y);
                                cell_temp(1,1)=x-1;
                                cell_temp(1,5)=y;
                                cell_temp(1,2)=cell_temp(1,1);
                                cell_temp(1,3)=cell_temp(1,1)+1;
                                cell_temp(1,4)=cell_temp(1,1)+1;
                                cell_temp(1,6)=cell_temp(1,5)+1;
                                cell_temp(1,7)=cell_temp(1,5)+1;
                                cell_temp(1,8)=cell_temp(1,5);
                                cell_temp(1,14)=stages[row];
                                cell_temp(1,22)=viability[row];
                                cell_temp(1,24)=migration_follow_flags[row];
                                Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                                Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                                Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                                Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                                
                            }
                            else if (proloci2==4)
                            {
                                Visual_range.clear_square(x,y);
                                cell_temp(1,1)=x;
                                cell_temp(1,5)=y+1;
                                cell_temp(1,2)=cell_temp(1,1);
                                cell_temp(1,3)=cell_temp(1,1)+1;
                                cell_temp(1,4)=cell_temp(1,1)+1;
                                cell_temp(1,6)=cell_temp(1,5)+1;
                                cell_temp(1,7)=cell_temp(1,5)+1;
                                cell_temp(1,8)=cell_temp(1,5);
                                cell_temp(1,14)=stages[row];
                                cell_temp(1,22)=viability[row];
                                cell_temp(1,24)=migration_follow_flags[row];
                                Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                                Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                                Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                                Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                            }
                            else if (proloci2==6)
                            {
                                Visual_range.clear_square(x,y);
                                cell_temp(1,1)=x+1;
                                cell_temp(1,5)=y;
                                cell_temp(1,2)=cell_temp(1,1);
                                cell_temp(1,3)=cell_temp(1,1)+1;
                                cell_temp(1,4)=cell_temp(1,1)+1;
                                cell_temp(1,6)=cell_temp(1,5)+1;
                                cell_temp(1,7)=cell_temp(1,5)+1;
                                cell_temp(1,8)=cell_temp(1,5);
                                cell_temp(1,14)=stages[row];
                                cell_temp(1,22)=viability[row];
                                cell_temp(1,24)=migration_follow_flags[row];
                                Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                                Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                                Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                                Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                            }
                            else if (proloci2==8)
                            {
                                Visual_range.clear_square(x,y);
                                cell_temp(1,1)=x;
                                cell_temp(1,5)=y-1;
                                cell_temp(1,2)=cell_temp(1,1);
                                cell_temp(1,3)=cell_temp(1,1)+1;
                                cell_temp(1,4)=cell_temp(1,1)+1;
                                cell_temp(1,6)=cell_temp(1,5)+1;
                                cell_temp(1,7)=cell_temp(1,5)+1;
                                cell_temp(1,8)=cell_temp(1,5);
                                cell_temp(1,14)=stages[row];
                                cell_temp(1,22)=viability[row];
                                cell_temp(1,24)=migration_follow_flags[row];
                                Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                                Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                                Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                                Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                            }
                        }
                        else if (random_pro_loci[0]==2)
                        {
                            int proloci1=pro_loci2_new[0];
                            int proloci2=pro_loci1_new[0];
                            if (proloci1==2)
                            {
                                Visual_range.clear_square(x,y);
                                for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                {
                                    cell_array.column(cor_cell)[row]=cell_array.column(cor_cell)[row]-1;
                                }
                                Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                                Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                                Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                                Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                            }
                            else if (proloci1==4)
                            {
                                Visual_range.clear_square(x,y);
                                for (int cor_cell_y=5; cor_cell_y<=8;cor_cell_y++)
                                {
                                    cell_array.column(cor_cell_y)[row]=cell_array.column(cor_cell_y)[row]+1;
                                }
                                Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                                Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                                Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                                Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                            }
                            else if (proloci1==6)
                            {
                                Visual_range.clear_square(x,y);
                                for (int cor_cell=1; cor_cell<=4;cor_cell++)
                                {
                                    cell_array.column(cor_cell)[row]=cell_array.column(cor_cell)[row]+1;
                                }
                                Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                                Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                                Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                                Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                            }
                            else if (proloci1==8)
                            {
                                Visual_range.clear_square(x,y);
                                for (int cor_cell_y=5; cor_cell_y<=8;cor_cell_y++)
                                {
                                    cell_array.column(cor_cell_y)[row]=cell_array.column(cor_cell_y)[row]-1;
                                }
                                Visual_range.set_square_occupied(x1_values[row],y1_values[row],1);
                                Visual_range.set_square_density_label(x1_values[row],y1_values[row],(int)ids[row]);
                                Visual_range.set_square_stage(x1_values[row],y1_values[row],cellstage);
                                Visual_range.set_square_cell_label(x1_values[row],y1_values[row],cell_label_1);
                            }
                            if (proloci2==1)
                            {
                                Visual_range.clear_square(x,y);
                                cell_temp(1,1)=x-1;
                                cell_temp(1,5)=y-1;
                                cell_temp(1,2)=cell_temp(1,1);
                                cell_temp(1,3)=cell_temp(1,1)+1;
                                cell_temp(1,4)=cell_temp(1,1)+1;
                                cell_temp(1,6)=cell_temp(1,5)+1;
                                cell_temp(1,7)=cell_temp(1,5)+1;
                                cell_temp(1,8)=cell_temp(1,5);
                                cell_temp(1,14)=stages[row];
                                cell_temp(1,22)=viability[row];
                                cell_temp(1,24)=migration_follow_flags[row];
                                Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                                Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                                Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                                Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                            }
                            else if (proloci2==3)
                            {
                                Visual_range.clear_square(x,y);
                                cell_temp(1,1)=x-1;
                                cell_temp(1,5)=y+1;
                                cell_temp(1,2)=cell_temp(1,1);
                                cell_temp(1,3)=cell_temp(1,1)+1;
                                cell_temp(1,4)=cell_temp(1,1)+1;
                                cell_temp(1,6)=cell_temp(1,5)+1;
                                cell_temp(1,7)=cell_temp(1,5)+1;
                                cell_temp(1,8)=cell_temp(1,5);
                                cell_temp(1,14)=stages[row];
                                cell_temp(1,22)=viability[row];
                                cell_temp(1,24)=migration_follow_flags[row];
                                Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                                Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                                Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                                Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                            }
                            else if (proloci2==5)
                            {
                                Visual_range.clear_square(x,y);
                                cell_temp(1,1)=x+1;
                                cell_temp(1,5)=y+1;
                                cell_temp(1,2)=cell_temp(1,1);
                                cell_temp(1,3)=cell_temp(1,1)+1;
                                cell_temp(1,4)=cell_temp(1,1)+1;
                                cell_temp(1,6)=cell_temp(1,5)+1;
                                cell_temp(1,7)=cell_temp(1,5)+1;
                                cell_temp(1,8)=cell_temp(1,5);
                                cell_temp(1,14)=stages[row];
                                cell_temp(1,22)=viability[row];
                                cell_temp(1,24)=migration_follow_flags[row];
                                Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                                Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                                Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                                Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                            }
                            else if (proloci2==7)
                            {
                                Visual_range.clear_square(x,y);
                                cell_temp(1,1)=x+1;
                                cell_temp(1,5)=y-1;
                                cell_temp(1,2)=cell_temp(1,1);
                                cell_temp(1,3)=cell_temp(1,1)+1;
                                cell_temp(1,4)=cell_temp(1,1)+1;
                                cell_temp(1,6)=cell_temp(1,5)+1;
                                cell_temp(1,7)=cell_temp(1,5)+1;
                                cell_temp(1,8)=cell_temp(1,5);
                                cell_temp(1,14)=stages[row];
                                cell_temp(1,22)=viability[row];
                                cell_temp(1,24)=migration_follow_flags[row];
                                Visual_range.set_square_occupied(cell_temp(1,1),cell_temp(1,5),1);
                                Visual_range.set_square_density_label(cell_temp(1,1),cell_temp(1,5),(int)ids[row]);
                                Visual_range.set_square_stage(cell_temp(1,1),cell_temp(1,5),cellstage);
                                Visual_range.set_square_cell_label(cell_temp(1,1),cell_temp(1,5),cell_label);
                            }
                        }
                    }
                }
                else
                {
                    IntMatrix cor_pro_1_1(2, 16);
                    cor_pro_1_1=0;
                    int ssss=1;
                    for (int sss=1;sss<=16;sss++)
                    {
                        if (Visual_range.occupied(cor_big_1_change_shape(1,sss),cor_big_1_change_shape(2,sss))==0)
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
                    IntMatrix cor_pro_1(2, cor_pro_1_1_nozero_length + 4);
                    cor_pro_1=0;
                    int nzl=1;
                    for (int aa=1;aa<=16;aa++ )
                    {
                        if (cor_pro_1_1(1,aa)!=0)
                        {
                            cor_pro_1.copy_col_from(nzl, cor_pro_1_1, aa);
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
                    stateless_shuffle(random_pro_loci_1,random_pro_loci_1+cor_pro_1_length, cell_rng_id, rng_time_step, rng_event++);
                    Visual_range.clear_square(x,y);
                    cell_temp(1,1)=(double)cor_pro_1(1,random_pro_loci_1[0]);
                    cell_temp(1,5)=(double)cor_pro_1(2,random_pro_loci_1[0]);
                    cell_temp(1,14)=1;
                    cell_temp(1,15)=ids[row];
                    x1_values[row]=(double)cor_pro_1(1,random_pro_loci_1[1]);
                    y1_values[row]=(double)cor_pro_1(2,random_pro_loci_1[1]);
                    if (cell_temp(1,1)!=0 && cell_temp(1,5)!=0)
                    {
                        Visual_range.occupied((int)cell_temp(1,1),(int)cell_temp(1,5))=1;
                        Visual_range.density_label((int)cell_temp(1,1),(int)cell_temp(1,5))=(int)ids[row];
                        Visual_range.stage((int)cell_temp(1,1),(int)cell_temp(1,5))=cellstage;
                        Visual_range.cell_label((int)cell_temp(1,1),(int)cell_temp(1,5))=cell_label;
                        Visual_range.occupied((int)x1_values[row],(int)y1_values[row])=1;
                        Visual_range.density_label((int)x1_values[row],(int)y1_values[row])=(int)ids[row];
                        Visual_range.stage((int)x1_values[row],(int)y1_values[row])=cellstage;
                        Visual_range.cell_label((int)x1_values[row],(int)y1_values[row])=cell_label_1;
                    }
                    stages[row]=1;
                    x2_values[row]=0;
                    x3_values[row]=0;
                    x4_values[row]=0;
                    y2_values[row]=0;
                    y3_values[row]=0;
                    y4_values[row]=0;
                    delete[] random_pro_loci_1;
                    random_pro_loci_1 = NULL;
                }
            }
            else
            {
                migration_elapsed[row]=0;
                cell_temp(1,20)=0;
                cell_temp(1,16)=0;
                division_elapsed[row]=0;
                long cell_label_1=Visual_range.cell_label(x,y);
                long cellstage=Visual_range.stage(x,y);
                cell_label=cell_label+1;
                IntMatrix pro_loci_small_1(2, 16);
                pro_loci_small_1=0;
                int ddddd=1;
                for (int loci22=1;loci22<=16;loci22++)
                {
                    if (Visual_range.occupied(cor_big_1_change_shape(1,loci22),cor_big_1_change_shape(2,loci22))==0)
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
                IntMatrix pro_loci_small(2, pro_loci_small_1_nozero_length + 4);
                pro_loci_small=0;
                int nzl=1;
                for (int aa=1;aa<=16;aa++ )
                {
                    if (pro_loci_small_1(1,aa)!=0)
                    {
                        pro_loci_small.copy_col_from(nzl, pro_loci_small_1, aa);
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
                stateless_shuffle(random_pro_loci_1,random_pro_loci_1+pro_loci_small_length, cell_rng_id, rng_time_step, rng_event++);
                Visual_range.clear_square(x,y);
                cell_temp(1,1)=(double)pro_loci_small(1,random_pro_loci_1[0]);
                cell_temp(1,5)=(double)pro_loci_small(2,random_pro_loci_1[0]);
                cell_temp(1,14)=1;
                cell_temp(1,15)=ids[row];
                x1_values[row]=(double)pro_loci_small(1,random_pro_loci_1[1]);
                y1_values[row]=(double)pro_loci_small(2,random_pro_loci_1[1]);
                if (cell_temp(1,1)!=0 && cell_temp(1,5)!=0)
                {
                    Visual_range.occupied((int)cell_temp(1,1),(int)cell_temp(1,5))=1;
                    Visual_range.density_label((int)cell_temp(1,1),(int)cell_temp(1,5))=(int)ids[row];
                    Visual_range.stage((int)cell_temp(1,1),(int)cell_temp(1,5))=cellstage;
                    Visual_range.cell_label((int)cell_temp(1,1),(int)cell_temp(1,5))=cell_label;
                    Visual_range.occupied((int)x1_values[row],(int)y1_values[row])=1;
                    Visual_range.density_label((int)x1_values[row],(int)y1_values[row])=(int)ids[row];
                    Visual_range.stage((int)x1_values[row],(int)y1_values[row])=cellstage;
                    Visual_range.cell_label((int)x1_values[row],(int)y1_values[row])=cell_label_1;
                }
                stages[row]=1;
                x2_values[row]=0;
                x3_values[row]=0;
                x4_values[row]=0;
                y2_values[row]=0;
                y3_values[row]=0;
                y4_values[row]=0;
                delete[] random_pro_loci_1;
                random_pro_loci_1 = NULL;
            }
            delete[] pro_loci_new;
            pro_loci_new = NULL;
        }
    }
    else if (stages[row]==1)
    {
        double growth_rate_inherent=growth_rates[row];
        double X1=growth_rate_inherent*(1-0.05);
        double X2=growth_rate_inherent*(1+0.05);
        growth_rates[row]=(X2-X1)*stateless_uniform(cell_rng_id, rng_time_step, rng_event++)+X1;
        cell_temp(1,9)=types[row];
        cell_temp(1,10)=growth_rates[row];
        cell_temp(1,11)=density_growth_rates[row];
        cell_temp(1,12)=migration_rate_bases[row];
        cell_temp(1,13)=random_labels[row];
        cell_temp(1,15)=ids[row];
        cell_temp(1,18)=death_times[row];
        cell_temp(1,19)=death_elapsed[row];
        cell_temp(1,21)=migration_intervals[row];
        cell_temp(1,22)=viability[row];
        cell_temp(1,25)=migration_active[row];
        cell_temp(1,26)=migration_duration[row];
        cell_temp(1,27)=migration_passed[row];
        cell_temp(1,28)=migration_rates[row];
        IntMatrix cor_temp_2(1, 8);
        cor_temp_2=0;
        int cor_temp_length=1;
        for (int s=1;s<=8;s++)
        {
            int x1=cor_small_1(1,s);
            int y1=cor_small_1(2,s);
            if (Visual_range.occupied(x1,y1)==0)
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
            long cellstage=Visual_range.stage(x,y);
            stateless_shuffle(cor_temp, cor_temp+loci_number, cell_rng_id, rng_time_step, rng_event++);
            int loci=cor_temp_3[cor_temp[0]];
            cell_temp(1,1)=cor_small_1(1,loci);
            cell_temp(1,5)=cor_small_1(2,loci);
            long cell_index=ids[row];
            Visual_range.occupied((int)cell_temp(1,1),(int)cell_temp(1,5))=1;
            Visual_range.density_label((int)cell_temp(1,1),(int)cell_temp(1,5))=cell_index;
            Visual_range.stage((int)cell_temp(1,1),(int)cell_temp(1,5))=cellstage;
            Visual_range.cell_label((int)cell_temp(1,1),(int)cell_temp(1,5))=cell_label;
            cell_temp(1,14)=1;
            migration_elapsed[row]=0;
            cell_temp(1,20)=0;
            cell_temp(1,16)=0;
            division_elapsed[row]=0;
        }
        else
        {
            if (types[row]==1)
            {
                viability[row]=0;
                Visual_range.clear_site(x1_values[row],y1_values[row]);
                cell_temp(1,22)=0;
            }
            else
            {
                if (utralsmall==1)
                {
                cell_temp(1,1)=x1_values[row];
                cell_temp(1,5)=y1_values[row];
                stages[row]=2;
                cell_temp(1,14)=2;
                Visual_range.stage((int)cell_temp(1,1),(int)cell_temp(1,5))=2;
                Visual_range.cell_label((int)cell_temp(1,1),(int)cell_temp(1,5))=cell_label;
                migration_elapsed[row]=0;
                cell_temp(1,20)=0;
                cell_temp(1,16)=0;
                division_elapsed[row]=0;
                }
            }
        }
        delete[] cor_temp;
        cor_temp = NULL;
        delete[] cor_temp_3;
        cor_temp_3 = NULL;
    }
    else if (stages[row]==2)
    {
        int cor_temp_length=0;
        for (int s=1;s<=8;s++)
        {
            int x1=cor_small_1(1,s);
            int y1=cor_small_1(2,s);
            if (Visual_range.occupied(x1,y1)==0)
            {
                cor_temp_length=cor_temp_length+1;
            }
        }
        if (cor_temp_length==0)
        {
            viability[row]=0;
        }
        else
        {
            division_elapsed[row]=division_elapsed[row]+deltah;
        }
    }
    if (types[row]==1)
    {
        if (growth_rates[row] > max_growth_rate_r)
        {
            growth_rates[row] = max_growth_rate_r;
        }
    }
    else if(types[row]==2)
    {
        if(growth_rates[row] > max_growth_rate_K)
        {
            growth_rates[row] = max_growth_rate_K;
        }
    }
    if (cell_temp(1,9)==1)
    {
        if(cell_temp(1,10) > max_growth_rate_r)
        {
            cell_temp(1,10) = max_growth_rate_r;
        }
    }
    else if (cell_temp(1,9)==2)
    {
        if(cell_temp(1,10) > max_growth_rate_K)
        {
            cell_temp(1,10) = max_growth_rate_K;
        }
    }
    if (cell_temp(1,1)!=0 && cell_temp(1,5)!=0)
    {
        division_detail::append_cell(cell_array, cell_temp, Col);
        
    }
}
#endif /* division_hpp */
