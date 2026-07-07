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
#include "cell_columns.hpp"
#include "cell_store.hpp"
#include "visual_range.hpp"
using namespace blitz;

inline void set_outer_visual_range_row(VisualRange &Visual_range, int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4, int cell_array_index, int cell_array_stage, int cell_label)
{
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
}

inline VisualRange outer_initiation_visualrange(const CellStore &cells,int N0,int Vx,int Vy,int &cell_label)
{
    VisualRange Visual_range(Vx,Vy);
    N0=cells.rows();
    for (int x=1; x<=N0; x++)
    {
        int row = x - 1;
        int x1 = cells.x1()[row];
        int y1 = cells.y1()[row];
        int x2 = cells.x2()[row];
        int y2 = cells.y2()[row];
        int x3 = cells.x3()[row];
        int y3 = cells.y3()[row];
        int x4 = cells.x4()[row];
        int y4 = cells.y4()[row];
        int cell_array_index=cells.id()[row];
        int cell_array_stage=cells.stage()[row];
        set_outer_visual_range_row(Visual_range, x1, y1, x2, y2, x3, y3, x4, y4, cell_array_index, cell_array_stage, cell_label);
        cell_label=cell_label+1;
    }
    return Visual_range;
}
#endif /* out_initiation_visualrange_hpp */
