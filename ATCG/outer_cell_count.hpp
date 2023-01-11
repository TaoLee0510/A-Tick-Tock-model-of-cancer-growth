//
//  outer_cell_count.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef outer_cell_count_hpp
#define outer_cell_count_hpp

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>

using namespace blitz;
void outer_cell_count(int Visual_range_x,int Visual_range_y,int &N0,double R0, double R1)
{
    double r0=R0/2;
    double r1=R1/2;
    for (int x=1;x<=Visual_range_x/2;x++)
    {
        for(int y=1;y<=Visual_range_y/2;y++)
        {
            if (pow((x-Visual_range_x/4),2)+pow((y-Visual_range_y/4),2)<=pow(r0,2) && pow((x-Visual_range_x/4),2)+pow((y-Visual_range_y/4),2)>=pow(r1,2))
            {
                N0=N0+1;
            }
        }
    }
}
#endif /* outer_cell_count_hpp */
