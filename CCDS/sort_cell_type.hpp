//
//  sort_cell_type.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef sort_cell_type_hpp
#define sort_cell_type_hpp

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "deltah_calculation.hpp"
using namespace blitz;
void sort_cell_type(Array<float, 2> &cell_array, Array<float, 2> cell_array1,int Col)
{
    Range all = Range::all();
    int C_9= cell_array.rows();
    double cell_array_9[C_9];
    for (int CN=0; CN<C_9; CN++)
    {
        cell_array_9[CN]=cell_array(CN+1,9);
    }
    
    struct array{
        double data;
        int index;
    };
    int num;
    num=C_9;
    array *p=new array[num];
    for (int i=0; i<num; i++)
    {
        p[i].data = cell_array_9[i];
        p[i].index = i+1;
        
    }
    qsort(p, sizeof(cell_array_9)/sizeof(cell_array_9[0]), sizeof(p[0]), compare);
    int index_9[C_9];
    for (int in=0; in<C_9; in++)
    {
        index_9[in]=p[in].index;
    }
    cell_array1.resize(C_9, Col);
    int number =1;
    for (int CNx=C_9; CNx>0; CNx--)
    {
        int ind=index_9[CNx-1];
        cell_array1(number,all)=cell_array(ind,all);
        number = number + 1;
    }
    cell_array(all,all)=0;
    cell_array(all,all)=cell_array1(all,all);
    delete[] p;
    p=NULL;
}
#endif /* sort_cell_type_hpp */
