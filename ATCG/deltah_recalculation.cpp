//
//  deltah_recalculation.hpp
//  CCDS
//
//  Created by Taolee on 2019/1/13.
//  Copyright © 2019 Tao Lee. All rights reserved.
//

#include "deltah_recalculation.hpp"

#include <stdio.h>
#include <algorithm>
#include <cmath>
#include "cell_columns.hpp"
#include "cell_store.hpp"

void set_deltah_from_max_migration_rate(double &deltah, int &MMR, double max_mig_r)
{
    double aT=1;
    MMR=(int)ceil(max_mig_r);
    if (MMR <= 0)
    {
        MMR = 1;
    }
    deltah = aT/(double)MMR;
    if (deltah>0.1)
    {
        deltah=0.1;
    }
}

void deltah_recalculation(double &deltah, const CellStore &cells, int &MMR, int DDM)
{
    if (DDM==1)
    {
        const CellStore::Column &migration_rates = cells.migration_rate_base();
        if (migration_rates.empty())
        {
            return;
        }
        double max_mig_r=*std::max_element(migration_rates.begin(), migration_rates.end());
        set_deltah_from_max_migration_rate(deltah, MMR, max_mig_r);
    }
}
