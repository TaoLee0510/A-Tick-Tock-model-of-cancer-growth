//
//  density_growth_rate_calculation_1.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#include "density_growth_rate_calculation_1.hpp"

#include <algorithm>
#include <stdio.h>
#include <vector>
#include <ctime>
#include "cell_columns.hpp"
#include "cell_store.hpp"
#include "visual_range.hpp"
#include "stateless_rng.hpp"


long unique_nonzero_count(long *values, int count)
{
    std::vector<long> unique_values(values, values + count);
    std::sort(unique_values.begin(), unique_values.end());
    unique_values.erase(std::unique(unique_values.begin(), unique_values.end()), unique_values.end());

    long result = unique_values.size();
    if (!unique_values.empty() && unique_values[0] == 0)
    {
        result = result - 1;
    }
    return result;
}

bool is_r_density_label(long label, int N00, int N01)
{
    return (label <= N00 / 2 && label >= 1) || (label <= N00 + (N01 / 2) && label > N00);
}

DensityGrowthCounts density_growth_neighborhood_counts(int x1, int y1, int N00, int N01, const VisualRange &Visual_range)
{
    long subcell_r[36]={0};
    long subcell_K[36]={0};
    int r_cells_number=0;
    int K_cells_number=0;
    long stage_sum=0;

    for(int cxx=1;cxx<=6;cxx++)
    {
        for(int cyy=1;cyy<=6;cyy++)
        {
            int visual_x = x1 - 3 + cxx;
            int visual_y = y1 - 3 + cyy;
            long density_label = Visual_range.density_label(visual_x, visual_y);
            long cell_label = Visual_range.cell_label(visual_x, visual_y);
            if(is_r_density_label(density_label, N00, N01))
            {
                subcell_r[r_cells_number]=cell_label;
                r_cells_number=r_cells_number+1;
            }
            else
            {
                subcell_K[K_cells_number]=cell_label;
                K_cells_number=K_cells_number+1;
            }

            if (Visual_range.stage(visual_x, visual_y)==2)
            {
                stage_sum=stage_sum+1;
            }
        }
    }

    DensityGrowthCounts counts;
    counts.rc = unique_nonzero_count(subcell_r, 36);
    counts.kc = unique_nonzero_count(subcell_K, 36) + stage_sum;
    counts.cells_number = counts.rc + counts.kc;
    return counts;
}

double calculate_density_growth_rate(int cell_type, double inherent_growth_rate, const DensityGrowthCounts &counts, double r_limit, double K_limit, double alpha, double beta, double carrying_capacity_r, double carrying_capacity_K)
{
    switch (cell_type)
    {
        case 1:
        {
            if (counts.cells_number < r_limit)
            {
                return inherent_growth_rate;
            }
            double crowding = counts.rc + counts.kc + alpha * counts.kc - r_limit;
            return inherent_growth_rate - (inherent_growth_rate * 2 * crowding) / carrying_capacity_r;
        }
        case 2:
        {
            if (counts.cells_number < K_limit)
            {
                return inherent_growth_rate;
            }
            double crowding = beta * counts.rc + counts.rc + counts.kc - K_limit;
            return inherent_growth_rate - (inherent_growth_rate * 2 * crowding) / carrying_capacity_K;
        }
        default:
        {
            return inherent_growth_rate;
        }
    }
}

void density_growth_rate_calculation_1(int Visual_range_x, int Visual_range_y, int N00, int N01, double r_limit, double K_limit, double lambda_r, double lambda_K, double alpha, double beta, double carrying_capacity_r, double carrying_capacity_K, double Cr, double CK, double death_time_range_r, double death_time_range_K, CellStore &cells,const VisualRange &Visual_range, long rng_time_step)
{
    (void)lambda_r;
    (void)lambda_K;
    (void)Cr;
    (void)CK;
    (void)death_time_range_r;
    (void)death_time_range_K;

    const CellStore::Column &ids = cells.id();
    const CellStore::Column &xs = cells.x1();
    const CellStore::Column &ys = cells.y1();
    const CellStore::Column &types = cells.type();
    const CellStore::Column &growth_rates = cells.growth_rate();
    CellStore::Column &density_growth_rates = cells.density_growth_rate();
    CellStore::Column &division_times = cells.division_time();
    CellStore::Column &migration_intervals = cells.migration_interval();
    const CellStore::Column &migration_rates = cells.migration_rate();

    int C0= cells.rows();
    for (int i=1; i<=C0; i++)
    {
        int row = i - 1;
        long cell_rng_id = (long)ids[row];
        if (cell_rng_id == 0)
        {
            cell_rng_id = i;
        }
        long rng_event = 100;
        int x1=(int)xs[row];
        int y1=(int)ys[row];
        int cell_type=(int)types[row];
        if (x1<100 || y1<100 || x1>Visual_range_x+100 || y1>Visual_range_y+100)
        {
            continue;
        }

        switch (cell_type)
        {
            case 1:
            {
                DensityGrowthCounts counts = density_growth_neighborhood_counts(x1, y1, N00, N01, Visual_range);
                double growth_rate_inherent_r=growth_rates[row];
                density_growth_rates[row]=calculate_density_growth_rate(cell_type, growth_rate_inherent_r, counts, r_limit, K_limit, alpha, beta, carrying_capacity_r, carrying_capacity_K);
                if (density_growth_rates[row]>0)
                {
                    double expected_division_time=24/density_growth_rates[row];
                    double undividing_time=0.9*expected_division_time;
                    double diving_time_range=0.1*expected_division_time;
                    double probability_of_division=1/diving_time_range;
                    double expected_dividing_time=undividing_time+stateless_geometric(cell_rng_id, rng_time_step, rng_event++, probability_of_division);;
                    division_times[row]=expected_dividing_time;
                }
                migration_intervals[row]=1/migration_rates[row];
                break;
            }
            case 2:
            {
                DensityGrowthCounts counts = density_growth_neighborhood_counts(x1, y1, N00, N01, Visual_range);
                double growth_rate_inherent_K=growth_rates[row];
                density_growth_rates[row]=calculate_density_growth_rate(cell_type, growth_rate_inherent_K, counts, r_limit, K_limit, alpha, beta, carrying_capacity_r, carrying_capacity_K);
                if (density_growth_rates[row]>0)
                {
                    double expected_division_time=24/density_growth_rates[row];
                    double undividing_time=0.9*expected_division_time;
                    double diving_time_range=0.1*expected_division_time;
                    double probability_of_division=1/diving_time_range;
                    double expected_dividing_time=undividing_time+stateless_geometric(cell_rng_id, rng_time_step, rng_event++, probability_of_division);;
                    division_times[row]=expected_dividing_time;
                }
                migration_intervals[row]=1/migration_rates[row];
                break;
            }
        }
    }
}
