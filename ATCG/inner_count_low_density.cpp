//
//  inner_count_low_density.hpp
//  CCDS
//
//  Created by Taolee on 2019/4/13.
//  Copyright © 2019 Tao Lee. All rights reserved.
//

#include "inner_count_low_density.hpp"

#include <stdio.h>
#include "visual_range.hpp"
int inner_count_low_density(int Visual_range_x, int Visual_range_y, const VisualRange &Visual_range, int N01, double R1)
{
    for (int x=1;x<=Visual_range_x;x++)
    {
        for(int y=1;y<=Visual_range_y;y++)
        {
            if (pow((x-Visual_range_x/2),2)+pow((y-Visual_range_y/2),2)<=pow(R1,2))
            {
                if (Visual_range.occupied(x,y)==0)
                {
                    N01=N01+1;
                }
            }
        }
    }
    return N01;
}
