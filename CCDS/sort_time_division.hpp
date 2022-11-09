//
//  sort_time_division.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef sort_time_division_hpp
#define sort_time_division_hpp

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "deltah_calculation.hpp"
using namespace blitz;
void sort_time_division(Array<float, 2> &cell_array, Array<float, 2> cell_array1)
{
    Range all = Range::all();
    int C_16= cell_array.rows();
    double cell_array_16[C_16];
    for (int CN=0; CN<C_16; CN++)
    {
        cell_array_16[CN]=cell_array(CN+1,16);
    }
    struct array{
        double data;
        int index;
    };
    int num;
    num=C_16;
    array *p=new array[num];
    for (int i=0; i<num; i++)
    {
        p[i].data = cell_array_16[i];
        p[i].index = i+1;
    }
    qsort(p, sizeof(cell_array_16)/sizeof(cell_array_16[0]), sizeof(p[0]), compare);
    int index_16[C_16];
    for (int in=0; in<C_16; in++)
    {
        index_16[in]=p[in].index;
    }
    cell_array1.resize(C_16, 31);
    for (int CNx=0; CNx<C_16; CNx++)
    {
        int ind=index_16[CNx];
        cell_array1(CNx+1,all)=cell_array(ind,all);
    }
    cell_array(all,all)=0;
    cell_array(all,all)=cell_array1(all,all);
    delete[] p;
    p=NULL;
}
#endif /* sort_time_division_hpp */
