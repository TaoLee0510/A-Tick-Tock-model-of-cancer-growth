//
//  out_initiation_visualrange.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef out_initiation_visualrange_hpp
#define out_initiation_visualrange_hpp

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
using namespace blitz;
Array<int,3> outer_initiation_visualrange(Array<float,2> cell_array0,int N0,int Vx,int Vy,int &cell_label)
{
    Range all = Range::all();
    Array<int,3> Visual_range(Vx,Vy,4,FortranArray<3>());
    Visual_range(all,all,all)=0;
    N0=cell_array0.rows();
    for (int x=1; x<=N0; x++)
    {
        int x1 = cell_array0(x,1);
        int y1 = cell_array0(x,5);
        int x2 = cell_array0(x,2);
        int y2 = cell_array0(x,6);
        int x3 = cell_array0(x,3);
        int y3 = cell_array0(x,7);
        int x4 = cell_array0(x,4);
        int y4 = cell_array0(x,8);
        int cell_array_index=cell_array0(x,15);
        int cell_array_stage=cell_array0(x,14);
        Visual_range(x1,y1,1)=1;
        Visual_range(x1,y1,2)=cell_array_index;
        Visual_range(x1,y1,3)=cell_array_stage;
        Visual_range(x1,y1,4)=cell_label;
        Visual_range(x2,y2,1)=1;
        Visual_range(x2,y2,2)=cell_array_index;
        Visual_range(x2,y2,3)=cell_array_stage;
        Visual_range(x2,y2,4)=cell_label;
        Visual_range(x3,y3,1)=1;
        Visual_range(x3,y3,2)=cell_array_index;
        Visual_range(x3,y3,3)=cell_array_stage;
        Visual_range(x3,y3,4)=cell_label;
        Visual_range(x4,y4,1)=1;
        Visual_range(x4,y4,2)=cell_array_index;
        Visual_range(x4,y4,3)=cell_array_stage;
        Visual_range(x4,y4,4)=cell_label;
        cell_label=cell_label+1;
    }
    return Visual_range;
}
#endif /* out_initiation_visualrange_hpp */
