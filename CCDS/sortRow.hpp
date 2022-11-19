//
//  sortRow.hpp
//  CCDS
//
//  Created by Tao Lee on 11/16/22.
//  Copyright © 2022 Tao Lee. All rights reserved.
//

#ifndef sortRow_hpp
#define sortRow_hpp


#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "deltah_calculation.hpp"
using namespace blitz;
using namespace std;
void sortRow(Array<float, 2> &cell_array, Array<float, 2> cell_array1, int Col, int colnum_to_sort)
{
    Range all = Range::all();
    int C_16= cell_array.rows();
    double cell_array_16[C_16];
    
    struct array{
        double data;
        int index;
    };
    
    array *p=new array[C_16];
    
    
    for (int CN=0; CN<C_16; CN++)
    {
        cell_array_16[CN]=cell_array(CN+1,colnum_to_sort);
        p[CN].data = cell_array(CN+1,colnum_to_sort);
        p[CN].index = CN+1;
    }
    
    qsort(p, sizeof(cell_array_16)/sizeof(cell_array_16[0]), sizeof(p[0]), compare);
    
    cell_array1.resize(C_16, Col);
    
    for (int CNx=0; CNx<C_16; CNx++)
    {
        cell_array1(CNx+1,all)=cell_array(p[CNx].index,all);
    }
    
    cell_array(all,all)=0;
    cell_array(all,all)=cell_array1(all,all);
    delete[] p;
    p=NULL;
}

#endif /* sortRow_hpp */
