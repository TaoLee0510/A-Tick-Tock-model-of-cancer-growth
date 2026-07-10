//
//  CellDivision.hpp
//  CCDS
//
//  Created by Tao Lee on 12/7/22.
//  Copyright © 2022 Tao Lee. All rights reserved.
//

#include "CellDivision.hpp"

#include <stdio.h>
#include <omp.h>
#include <iostream>
#include <time.h>
#include <memory>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <getopt.h>
#include "visual_range.hpp"
#include "outer_corr.hpp"
#include "outer_cell_count.hpp"
#include "deltah_calculation.hpp"
#include "outer_initiation_array.hpp"
#include "out_initiation_visualrange.hpp"
#include "inner_count.hpp"
#include "inner_initiation_array.hpp"
#include "density_growth_rate_calculation_1.hpp"
#include "stage_convert.hpp"
#include "death_judgement.hpp"
#include "save_data_free_living.hpp"
#include "random_migration.hpp"
#include "migration.hpp"
#include "division.hpp"
#include "free_living_division.hpp"
#include "migrate_activation.hpp"
#include "density_calculation.hpp"
#include "deltah_recalculation.hpp"
#include <omp.h>

void CellDivision(int i, double max_growth_rate_r, double max_growth_rate_K, CellStore &cell_array, VisualRange &Visual_range, CellRowBuffer &cell_temp,int &cell_label, double &deltah,int utralsmall,int Col,double deathjudge,int borderx,int bordery, long rng_time_step)
{
    if(cell_array.x1()[i - 1]==0 && cell_array.y1()[i - 1] ==0)
    {
        cell_array.viability()[i - 1]=0;
    }
    if (cell_array.x1()[i - 1]>=100 && cell_array.y1()[i - 1] >=100 && cell_array.x1()[i - 1]<=borderx && cell_array.y1()[i - 1]<=bordery)
    {
//        if (cell_array.density_growth_rate()[i - 1]>deathjudge)
        if (cell_array.division_time()[i - 1]>0)
        {
            if (cell_array.division_elapsed()[i - 1]>=cell_array.division_time()[i - 1])
            {
                division(i, max_growth_rate_r, max_growth_rate_K, cell_array, Visual_range, cell_temp,cell_label,deltah,utralsmall, Col, rng_time_step);
            }
        }
    }
}
