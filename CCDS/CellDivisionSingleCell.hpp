//
//  CellDivisionSingleCell.hpp
//  CCDS
//
//  Created by Tao Lee on 12/7/22.
//  Copyright © 2022 Tao Lee. All rights reserved.
//

#ifndef CellDivisionSingleCell_hpp
#define CellDivisionSingleCell_hpp

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
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_matrix.h>
#define BZ_THREADSAFE
#define BZ_THREADSAFE_USE_OPENMP
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "random_uniform.hpp"
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
#include "free_living_division.hpp"
#include "migrate_activation.hpp"
#include "density_calculation.hpp"
#include "deltah_recalculation.hpp"
#include "sortRow.hpp"
#include <omp.h>

void CellDivisionSingleCell(int i, double max_growth_rate_r, double max_growth_rate_K, Array<float, 2> &cell_array, Array<float,2> cell_array_temp, Array<int, 3> &Visual_range, Array<int,2> cor_big_1, Array<int, 2> cor_big_1_change_shape, Array<int, 2> cor_small_1, Array<int, 2> proliferation_loci, Array<float, 2> cell_temp,int &cell_label, double &deltah,int utralsmall, double beta_distribution_alpha_for_normal_migration,double beta_distribution_beta_for_normal_migration,double migration_rate_K_mean,double uniup_K, double unilow_K,double sigmahatK,double muhatK,int &K_label,Array<int, 3> sub_visual,double beta_distribution_alpha, double beta_distribution_beta, double migration_rate_r_mean,double migration_rate_r_mean_quia,double beta_distribution_expected_for_normal_migration,Array<float,2> &cell_trace,Array<float,2> cell_trace_temp, int &cell_index,int &generation,int &r_label,int Col,double K_formation_rate,double deathjudge,int Visual_range_x,int Visual_range_y)
{
    if(cell_array(i,1)==0 && cell_array(i,5) ==0)
    {
        cell_array(i,22)=0;
    }
    if (cell_array(i,1)>=100 && cell_array(i,5) >=100 && cell_array(i,1)<=Visual_range_x+100 && cell_array(i,5)<=Visual_range_y+100)
    {
        if (cell_array(i,11)>deathjudge)
        {
            if (cell_array(i,16)>=cell_array(i,17))
            {
                free_living_division(i, max_growth_rate_r, max_growth_rate_K, cell_array, cell_array_temp, Visual_range, cor_big_1, cor_big_1_change_shape, cor_small_1, proliferation_loci, cell_temp,cell_label,deltah,utralsmall,beta_distribution_alpha_for_normal_migration, beta_distribution_beta_for_normal_migration,migration_rate_K_mean, uniup_K,unilow_K,sigmahatK,muhatK,K_label,sub_visual,beta_distribution_alpha, beta_distribution_beta, migration_rate_r_mean, migration_rate_r_mean_quia, beta_distribution_expected_for_normal_migration,cell_trace,cell_trace_temp,cell_index,generation,r_label,Col, K_formation_rate);
            }
        }
    }
}




#endif /* CellDivisionSingleCell_hpp */
