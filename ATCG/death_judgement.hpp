//
//  death_judgement.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef death_judgement_hpp
#define death_judgement_hpp

#include <stdio.h>
#include <omp.h>
#include <algorithm>
#include <random>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "cell_columns.hpp"
#include "cell_store.hpp"
#include "visual_range.hpp"
#include "density_growth_rate_calculation_1.hpp"
#include "stateless_rng.hpp"
using namespace blitz;

template <typename CellArray>
inline bool death_judgement_live_cell(const CellArray &cell_array, int row)
{
    int idx = row - 1;
    return cell_array.viability()[idx]==1 && cell_array.x1()[idx]!=0 && cell_array.y1()[idx]!=0;
}

template <typename CellArray>
inline void clear_dead_cell_visual(int site, CellArray &cell_array, VisualRange &Visual_range, int C)
{
    int idx = site - 1;
    auto &x1_values = cell_array.x1();
    auto &x2_values = cell_array.x2();
    auto &x3_values = cell_array.x3();
    auto &x4_values = cell_array.x4();
    auto &y1_values = cell_array.y1();
    auto &y2_values = cell_array.y2();
    auto &y3_values = cell_array.y3();
    auto &y4_values = cell_array.y4();
    auto &stages = cell_array.stage();
    int stage = (int)stages[idx];
    if(stage==0)
    {
        Visual_range.clear_site((int)x1_values[idx],(int)y1_values[idx]);
        Visual_range.clear_site((int)x2_values[idx],(int)y2_values[idx]);
        Visual_range.clear_site((int)x3_values[idx],(int)y3_values[idx]);
        Visual_range.clear_site((int)x4_values[idx],(int)y4_values[idx]);
    }
    else if(stage==1)
    {
        Visual_range.clear_site((int)x1_values[idx],(int)y1_values[idx]);
    }
    else if(stage==2)
    {
        double usx=x1_values[idx];
        double usy=y1_values[idx];
        x1_values[idx]=1;
        y1_values[idx]=1;
        for (int us=1;us<=C;us++)
        {
            int us_idx = us - 1;
            if(x1_values[us_idx]!=0 && y1_values[us_idx]!=0 && stages[idx]==2)
            {
                if(x1_values[us_idx]==usx && y1_values[us_idx]==usy)
                {
                    Visual_range.stage((int)x1_values[us_idx],(int)y1_values[us_idx])=1;
                    stages[us_idx]=1;
                }
            }
        }
    }
}

inline void compact_after_death_judgement(CellStore &cells, VisualRange &Visual_range, int, int nthreads, int C)
{
    int current_size=cells.rows();
    omp_set_num_threads(nthreads);
    int sum =0;
    #pragma omp parallel for schedule(dynamic) reduction(+:sum)
    {
        for (int CN=1; CN<=current_size; ++CN)
        {
            if (death_judgement_live_cell(cells, CN))
            {
                sum=sum+1;
            }
        }
    }
    CellStore compacted(cells.column_count());
    compacted.reserve(sum);
    for (int site=1; site<= current_size; ++site)
    {
        if (death_judgement_live_cell(cells, site))
        {
            compacted.append_row_from(cells, site);
        }
        else
        {
            clear_dead_cell_visual(site, cells, Visual_range, C);
        }
    }
    cells = compacted;
}

template <typename CellArray>
inline void death_judgement(int Visual_range_x, int Visual_range_y, int N00, int N01, double r_limit, double K_limit, double lambda_r, double lambda_K, double alpha, double beta, double carrying_capacity_r, double carrying_capacity_K, double Cr, double CK, double death_time_range_r, double death_time_range_K, double deltah, double &h, CellArray &cell_array, VisualRange &Visual_range, double deathjudge, int Col,int nthreads,long rng_time_step)
{
    int C= cell_array.rows();
    auto &x1_values = cell_array.x1();
    auto &y1_values = cell_array.y1();
    auto &types = cell_array.type();
    auto &growth_rates = cell_array.growth_rate();
    auto &density_growth_rates = cell_array.density_growth_rate();
    auto &ids = cell_array.id();
    auto &division_times = cell_array.division_time();
    auto &death_times = cell_array.death_time();
    auto &death_elapsed = cell_array.death_elapsed();
    auto &migration_intervals = cell_array.migration_interval();
    auto &viability = cell_array.viability();
    auto &migration_rates = cell_array.migration_rate();
//    omp_set_num_threads(nthreads);
//    #pragma omp parallel for schedule(dynamic)
//    {
        for (int rows=1; rows<=C; ++rows)
        {
            int row = rows - 1;
            long cell_rng_id = (long)ids[row];
            if (cell_rng_id == 0)
            {
                cell_rng_id = rows;
            }
            long rng_event = 100;
            if (density_growth_rates[row]>deathjudge)
            {
                if (types[row]==1)
                {
                    if (x1_values[row]>=100 && y1_values[row] >=100 && x1_values[row]<=Visual_range_x+100 && y1_values[row]<=Visual_range_y+100)
                    {
                        DensityGrowthCounts counts = density_growth_neighborhood_counts((int)x1_values[row], (int)y1_values[row], N00, N01, Visual_range);
                        long rc=counts.rc;
                        long kc=counts.kc;
                        long cells_number=counts.cells_number;
                        double growth_rate_inherent_r=growth_rates[row];
                        if (cells_number>=r_limit)
                        {
                            density_growth_rates[row]=growth_rate_inherent_r-((growth_rate_inherent_r*2*(rc+kc+alpha*kc-r_limit))/carrying_capacity_r);
                        }
                        else
                        {
                            density_growth_rates[row]=growth_rate_inherent_r;
                        }
                        if (density_growth_rates[row]>deathjudge)
                        {
                            double expected_division_time=24/density_growth_rates[row];
                            double undividing_time=0.9*expected_division_time;
                            double diving_time_range=0.1*expected_division_time;
                            double probability_of_division=1/diving_time_range;
                            double expected_dividing_time=undividing_time+stateless_geometric(cell_rng_id, rng_time_step, rng_event++, probability_of_division);;
                            division_times[row]=expected_dividing_time;
                        }
                        else
                        {
                            if (death_times[row]==0)
                            {
                                double probability_to_death=1/death_time_range_r;
                                death_times[row]=stateless_geometric(cell_rng_id, rng_time_step, rng_event++, probability_to_death);
                            }
                            death_elapsed[row]=death_elapsed[row]+deltah;
                            division_times[row]=0;
                        }
                        migration_intervals[row]=1/migration_rates[row];
                    }
                    migration_intervals[row]=1/migration_rates[row];
                }
                else
                {
                    if (x1_values[row]>=100 && y1_values[row] >=100 && x1_values[row]<=Visual_range_x+100 && y1_values[row]<=Visual_range_y+100)
                    {
                        DensityGrowthCounts counts = density_growth_neighborhood_counts((int)x1_values[row], (int)y1_values[row], N00, N01, Visual_range);
                        long rc=counts.rc;
                        long kc=counts.kc;
                        long cells_number=counts.cells_number;
                        double growth_rate_inherent_K=growth_rates[row];
                        if (cells_number>=K_limit)
                        {
                            density_growth_rates[row]=growth_rate_inherent_K-((growth_rate_inherent_K*2*(beta*rc+rc+kc-K_limit))/carrying_capacity_K);
                        }
                        else
                        {
                            density_growth_rates[row]=growth_rate_inherent_K;
                        }
                        if (density_growth_rates[row]>deathjudge)
                        {
                            double expected_division_time=24/density_growth_rates[row];
                            double undividing_time=0.9*expected_division_time;
                            double diving_time_range=0.1*expected_division_time;
                            double probability_of_division=1/diving_time_range;
                            double expected_dividing_time=undividing_time+stateless_geometric(cell_rng_id, rng_time_step, rng_event++, probability_of_division);;
                            division_times[row]=expected_dividing_time;
                        }
                        else
                        {
                            if (death_times[row]==0)
                            {
                                double probability_to_death=1/death_time_range_K;
                                death_times[row]=stateless_geometric(cell_rng_id, rng_time_step, rng_event++, probability_to_death);
                            }
                            death_elapsed[row]=death_elapsed[row]+deltah;
                            division_times[row]=0;
                        }
                        migration_intervals[row]=1/migration_rates[row];
                    }
                    migration_intervals[row]=1/migration_rates[row];
                }
            }
            else if (density_growth_rates[row]<=deathjudge)
            {
                if(h==0)
                {
                    double probability_to_death=0;
                    if (types[row]==1)
                    {
                        probability_to_death=1/death_time_range_r;
                    }
                    else if (types[row]==2)
                    {
                        probability_to_death=1/death_time_range_K;
                    }
                    death_times[row]=stateless_geometric(cell_rng_id, rng_time_step, rng_event++, probability_to_death);
                }
                else
                {
                    double D_time_1=1.5*(24/growth_rates[row]);
                    double D_time_2=0.9*death_times[row];
                    double D_time = 0;
                    if (D_time_1<=D_time_2)
                    {
                        D_time = D_time_1;
                    }
                    else
                    {
                        D_time = D_time_2;
                    }
                    if (death_elapsed[row]<=D_time)
                    {
                        if (types[row]==1)
                        {
                            if (x1_values[row]>=100 && y1_values[row] >=100 && x1_values[row]<=Visual_range_x+100 && y1_values[row]<=Visual_range_y+100)
                            {
                                DensityGrowthCounts counts = density_growth_neighborhood_counts((int)x1_values[row], (int)y1_values[row], N00, N01, Visual_range);
                                long rc=counts.rc;
                                long kc=counts.kc;
                                long cells_number=counts.cells_number;
                                double growth_rate_inherent_r=growth_rates[row];
                                if (cells_number>=r_limit)
                                {
                                    density_growth_rates[row]=growth_rate_inherent_r-((growth_rate_inherent_r*2*(rc+kc+alpha*kc-r_limit))/carrying_capacity_r);
                                }
                                else
                                {
                                    density_growth_rates[row]=growth_rate_inherent_r;
                                }
                                if (density_growth_rates[row]>deathjudge)
                                {
                                    double expected_division_time=24/density_growth_rates[row];
                                    double undividing_time=0.9*expected_division_time;
                                    double diving_time_range=0.1*expected_division_time;
                                    double probability_of_division=1/diving_time_range;
                                    double expected_dividing_time=undividing_time+stateless_geometric(cell_rng_id, rng_time_step, rng_event++, probability_of_division);;
                                    division_times[row]=expected_dividing_time;
                                    death_times[row]=0;
                                    death_elapsed[row]=0;
                                }
                                else
                                {
                                    death_elapsed[row]=death_elapsed[row]+deltah;
                                    division_times[row]=0;
                                }
                                migration_intervals[row]=1/migration_rates[row];
                            }
                            migration_intervals[row]=1/migration_rates[row];
                        }
                        else
                        {
                            if (x1_values[row]>=100 && y1_values[row] >=100 && x1_values[row]<=Visual_range_x+100 && y1_values[row]<=Visual_range_y+100)
                            {
                                DensityGrowthCounts counts = density_growth_neighborhood_counts((int)x1_values[row], (int)y1_values[row], N00, N01, Visual_range);
                                long rc=counts.rc;
                                long kc=counts.kc;
                                long cells_number=counts.cells_number;
                                double growth_rate_inherent_K=growth_rates[row];
                                if (cells_number>=K_limit)
                                {
                                    density_growth_rates[row]=growth_rate_inherent_K-((growth_rate_inherent_K*2*(beta*rc+rc+kc-K_limit))/carrying_capacity_K);
                                }
                                else
                                {
                                    density_growth_rates[row]=growth_rate_inherent_K;
                                }
                                if (density_growth_rates[row]>deathjudge)
                                {
                                    double expected_division_time=24/density_growth_rates[row];
                                    double undividing_time=0.9*expected_division_time;
                                    double diving_time_range=0.1*expected_division_time;
                                    double probability_of_division=1/diving_time_range;
                                    double expected_dividing_time=undividing_time+stateless_geometric(cell_rng_id, rng_time_step, rng_event++, probability_of_division);;
                                    division_times[row]=expected_dividing_time;
                                    death_times[row]=0;
                                    death_elapsed[row]=0;
                                }
                                else
                                {
                                    death_elapsed[row]=death_elapsed[row]+deltah;
                                    division_times[row]=0;
                                }
                                migration_intervals[row]=1/migration_rates[row];
                            }
                            migration_intervals[row]=1/migration_rates[row];
                        }
                    }
                    else
                    {
                        if(death_times[row]<=death_elapsed[row])
                        {
                            viability[row]=0;
                        }
                        else
                        {
                            death_elapsed[row]=death_elapsed[row]+deltah;
                        }
                    }
                }
            }
        }
//    }
    compact_after_death_judgement(cell_array, Visual_range, Col, nthreads, C);
}
#endif /* death_judgement_hpp */
