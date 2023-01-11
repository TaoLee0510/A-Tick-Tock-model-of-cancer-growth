//
//  outer_corr.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef outer_corr_hpp
#define outer_corr_hpp

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>

using namespace blitz;
void outer_corr(int Visual_range_x,int Visual_range_y,double R0, double R1,Array<int,2> &A)
{
    double r0=R0/2;
    double r1=R1/2;
    for (int x=1;x<=Visual_range_x/2;x++)
    {
        for(int y=1;y<=Visual_range_y/2;y++)
        {
            if (pow((x-Visual_range_x/4),2)+pow((y-Visual_range_y/4),2)<=pow(r0,2) && pow((x-Visual_range_x/4),2)+pow((y-Visual_range_y/4),2)>=pow(r1,2))
            {
                A(x,y)=1;
            }
            else
            {
                A(x,y)=0;
            }
        }
    }
}
#endif /* outer_corr_hpp */
