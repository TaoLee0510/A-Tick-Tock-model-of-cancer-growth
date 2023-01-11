//
//  deltah_recalculation.hpp
//  CCDS
//
//  Created by Taolee on 2019/1/13.
//  Copyright © 2019 Tao Lee. All rights reserved.
//

#ifndef deltah_recalculation_hpp
#define deltah_recalculation_hpp

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>


void deltah_recalculation(double &deltah, Array<double,2> cell_array, int &MMR, int DDM)
{
    if (DDM==1)
    {
        int C_12= cell_array.rows();
        double cell_array_12[C_12];
        for (int CN=0; CN<C_12; CN++)
        {
            cell_array_12[CN]=cell_array(CN+1,12);
        }
        qsort(cell_array_12,C_12,sizeof(cell_array_12[0]),cmp);
        double max_mig_r=cell_array_12[C_12-1];
        double aT=1;
        MMR=(int)ceil(max_mig_r);
        deltah = aT/(double)MMR;
        if (deltah>0.1)
        {
            deltah=0.1;
        }
    }
}
#endif /* deltah_recalculation_hpp */
