//
//  sort_time_per_generation.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef sort_time_per_generation_hpp
#define sort_time_per_generation_hpp

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "deltah_calculation.hpp"
using namespace blitz;
void sort_time_per_generation(Array<float, 2> &cell_array, Array<float, 2> cell_array1,int Col)
{
    Range all = Range::all();
    int C_17= cell_array.rows();
    double cell_array_17[C_17];
    for (int CN=0; CN<C_17; CN++)
    {
        cell_array_17[CN]=cell_array(CN+1,17);
    }
    struct array{
                double data;
                int index;
            };
    int num;
    num=C_17;
    array *p=new array[num];
    for (int i=0; i<num; i++)
    {
        p[i].data = cell_array_17[i];
         p[i].index = i+1;
        
    }
    qsort(p, sizeof(cell_array_17)/sizeof(cell_array_17[0]), sizeof(p[0]), compare);
    int index_17[C_17];
    for (int in=0; in<C_17; in++)
    {
        index_17[in]=p[in].index;
    }
    cell_array1.resize(C_17, Col);
    for (int CNx=0; CNx<C_17; CNx++)
    {
        int ind=index_17[CNx];
        cell_array1(CNx+1,all)=cell_array(ind,all);
    }
    cell_array(all,all)=0;
    cell_array(all,all)=cell_array1(all,all);
    delete[] p;
    p=NULL;
}
#endif /* sort_time_per_generation_hpp */
