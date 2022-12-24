//
//  random_migration.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef random_migration_hpp
#define random_migration_hpp

#include <stdio.h>
#include <random>
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
#include "deltah_calculation.hpp"
using namespace blitz;
void random_migration(int i, double deltah,Array<float, 2> &cell_array, Array<int, 3> &Visual_range, Array<int,2> cor_big, Array<int, 2> area_square, Array<int, 2> sub_area_square, Array<int, 2> cor_small, Array<int, 2> area_square_s, Array<int, 2>  sub_area_square_s,double &migration_judgement)
{
    //    std::random_device r;
    //    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    //    std::mt19937 RNG(seed);
    Range all = Range::all();
    const gsl_rng_type *T5;
    gsl_rng *r5;
    gsl_rng_env_setup();
    T5 = gsl_rng_ranlxs0;
    gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    r5 = gsl_rng_alloc(T5);
    int x1=cell_array(i,1);
    int y1=cell_array(i,5);
    if (cell_array(i,14)==0)
    {
        cor_big.resize(4,4);
        cor_big=0;
        cor_big(all,all)=Visual_range(Range(x1-1,x1+2),Range(y1-1,y1+2),1);
        int direction[8]={0};
        if (cor_big(2,1)==0 && cor_big(1,1)==0 && cor_big(1,2)==0)
        {
            direction[0]=1;
        }
        if (cor_big(1,2)==0 && cor_big(1,3)==0)
        {
            direction[1]=2;
        }
        if (cor_big(1,3)==0 && cor_big(1,4)==0 && cor_big(2,4)==0)
        {
            direction[2]=3;
        }
        if (cor_big(2,4)==0 && cor_big(3,4)==0)
        {
            direction[3]=4;
        }
        if (cor_big(3,4)==0 && cor_big(4,4)==0 && cor_big(4,3)==0)
        {
            direction[4]=5;
        }
        if (cor_big(4,3)==0 && cor_big(4,2)==0)
        {
            direction[5]=6;
        }
        if (cor_big(4,2)==0 && cor_big(4,1)==0 && cor_big(3,1)==0)
        {
            direction[6]=7;
        }
        if (cor_big(3,1)==0 && cor_big(2,1)==0)
        {
            direction[7]=8;
        }
        int mloci=0;
        for (int mlo=0; mlo<8; mlo++)
        {
            if (direction[mlo]!=0)
            {
                mloci++;
            }
        }
        
        if (mloci>0)
        {
            int *direction1=new int[mloci];
            int new_loci=0;
            for (int loci=0; loci<8; loci++)
            {
                if (direction[loci]!=0)
                {
                    direction1[new_loci]=direction[loci];
                    new_loci++;
                }
            }
            int order=0;
            long length_dir=new_loci;
            if (length_dir==0)
            {
                order=direction1[0];
            }
            else
            {
                //                shuffle(direction1, direction1+new_loci,RNG);
                gsl_ran_shuffle(r5, direction1, new_loci, sizeof (int));
                order=direction1[0];
            }
            
            switch (order)
            {
                case 1:
                {
                    int cell_label_1=Visual_range(x1,y1,4);
                    int cellstage=Visual_range(x1,y1,3);
                    Visual_range(Range(x1,x1+1),Range(y1,y1+1),all)=0;
                    for (int cor_cell=1; cor_cell<=4; cor_cell++)
                    {
                        int cor_cell_y=cor_cell+4;
                        cell_array(i,cor_cell)=cell_array(i,cor_cell)-1;
                        cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)-1;
                    }
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=cell_array(i,15);
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                    cell_array(i,23)=1;
                    cell_array(i,20)=0;
                    break;
                }
                case 2:
                {
                    int cell_label_1=Visual_range(x1,y1,4);
                    int cellstage=Visual_range(x1,y1,3);
                    Visual_range(Range(x1,x1+1),Range(y1,y1+1),all)=0;
                    for (int cor_cell=1; cor_cell<=4; cor_cell++)
                    {
                        cell_array(i,cor_cell)=cell_array(i,cor_cell)-1;
                    }
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=cell_array(i,15);
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                    cell_array(i,23)=2;
                    cell_array(i,20)=0;
                    break;
                }
                case 3:
                {
                    int cell_label_1=Visual_range(x1,y1,4);
                    int cellstage=Visual_range(x1,y1,3);
                    Visual_range(Range(x1,x1+1),Range(y1,y1+1),all)=0;
                    for (int cor_cell=1; cor_cell<=4; cor_cell++)
                    {
                        int cor_cell_y=cor_cell+4;
                        cell_array(i,cor_cell)=cell_array(i,cor_cell)-1;
                        cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)+1;
                    }
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=cell_array(i,15);
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                    cell_array(i,23)=3;
                    cell_array(i,20)=0;
                    break;
                }
                case 4:
                {
                    int cell_label_1=Visual_range(x1,y1,4);
                    int cellstage=Visual_range(x1,y1,3);
                    Visual_range(Range(x1,x1+1),Range(y1,y1+1),all)=0;
                    for (int cor_cell_y=5; cor_cell_y<=8; cor_cell_y++)
                    {
                        cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)+1;
                    }
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=cell_array(i,15);
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                    cell_array(i,23)=4;
                    cell_array(i,20)=0;
                    break;
                }
                case 5:
                {
                    int cell_label_1=Visual_range(x1,y1,4);
                    int cellstage=Visual_range(x1,y1,3);
                    Visual_range(Range(x1,x1+1),Range(y1,y1+1),all)=0;
                    for (int cor_cell=1; cor_cell<=4; cor_cell++)
                    {
                        int cor_cell_y=cor_cell+4;
                        cell_array(i,cor_cell)=cell_array(i,cor_cell)+1;
                        cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)+1;
                    }
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=cell_array(i,15);
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                    cell_array(i,23)=5;
                    cell_array(i,20)=0;
                    break;
                }
                case 6:
                {
                    int cell_label_1=Visual_range(x1,y1,4);
                    int cellstage=Visual_range(x1,y1,3);
                    Visual_range(Range(x1,x1+1),Range(y1,y1+1),all)=0;
                    for (int cor_cell=1; cor_cell<=4; cor_cell++)
                    {
                        cell_array(i,cor_cell)=cell_array(i,cor_cell)+1;
                    }
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=cell_array(i,15);
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                    cell_array(i,23)=6;
                    cell_array(i,20)=0;
                    break;
                }
                case 7:
                {
                    int cell_label_1=Visual_range(x1,y1,4);
                    int cellstage=Visual_range(x1,y1,3);
                    Visual_range(Range(x1,x1+1),Range(y1,y1+1),all)=0;
                    for (int cor_cell=1; cor_cell<=4; cor_cell++)
                    {
                        int cor_cell_y=cor_cell+4;
                        cell_array(i,cor_cell)=cell_array(i,cor_cell)+1;
                        cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)-1;
                    }
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=cell_array(i,15);
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                    cell_array(i,23)=7;
                    cell_array(i,20)=0;
                    break;
                }
                case 8:
                {
                    int cell_label_1=Visual_range(x1,y1,4);
                    int cellstage=Visual_range(x1,y1,3);
                    Visual_range(Range(x1,x1+1),Range(y1,y1+1),all)=0;
                    for (int cor_cell_y=5; cor_cell_y<=8; cor_cell_y++)
                    {
                        cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)-1;
                    }
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=cell_array(i,15);
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                    cell_array(i,23)=8;
                    cell_array(i,20)=0;
                    break;
                }
                    
            }
            delete[] direction1;
            direction1 = NULL;
            cell_array(i,20)=0;
            cell_array(i,24)=1;
        }
    }
    else //small stage
    {
        cor_small.resize(3, 3);
        cor_small=0;
        cor_small(all,all)=Visual_range(Range(x1-1,x1+1),Range(y1-1,y1+1),1);
        int direction[8]={0};
        if (cor_small(1,1)==0)
        {
            direction[0]=1;
        }
        if (cor_small(1,2)==0)
        {
            direction[1]=2;
        }
        if (cor_small(1,3)==0)
        {
            direction[2]=3;
        }
        if (cor_small(2,3)==0)
        {
            direction[3]=4;
        }
        if (cor_small(3,3)==0)
        {
            direction[4]=5;
        }
        if (cor_small(3,2)==0)
        {
            direction[5]=6;
        }
        if (cor_small(3,1)==0)
        {
            direction[6]=7;
        }
        if (cor_small(2,1)==0)
        {
            direction[7]=8;
        }
        int mloci=0;
        for (int mlo=0; mlo<8; mlo++)
        {
            if (direction[mlo]!=0)
            {
                mloci++;
            }
        }
        if (mloci>0)
        {
            int *direction1=new int[mloci];
            int new_loci=0;
            for (int loci=0; loci<8; loci++)
            {
                if (direction[loci]!=0)
                {
                    direction1[new_loci]=direction[loci];
                    new_loci++;
                }
            }
            int order=0;
            long length_dir=new_loci;
            if (length_dir==0)
            {
                order=direction1[0];
            }
            else if (length_dir>0)
            {
                //                shuffle(direction1, direction1+new_loci,RNG);
                gsl_ran_shuffle(r5, direction1, new_loci, sizeof (int));
                order=direction1[0];
            }
            switch (order)
            {
                case 1:
                {
                    int cell_label_1=Visual_range(x1,y1,4);
                    int cellstage=Visual_range(x1,y1,3);
                    Visual_range(x1,y1,all)=0;
                    Visual_range(x1-1,y1-1,1)=1;
                    Visual_range(x1-1,y1-1,2)=cell_array(i,15);
                    Visual_range(x1-1,y1-1,3)=cellstage;
                    Visual_range(x1-1,y1-1,4)=cell_label_1;
                    cell_array(i,1)=cell_array(i,1)-1;
                    cell_array(i,5)=cell_array(i,5)-1;
                    cell_array(i,23)=1;
                    cell_array(i,20)=0;
                    break;
                }
                case 2:
                {
                    int cell_label_1=Visual_range(x1,y1,4);
                    int cellstage=Visual_range(x1,y1,3);
                    Visual_range(x1,y1,all)=0;
                    Visual_range(x1-1,y1,1)=1;
                    Visual_range(x1-1,y1,2)=cell_array(i,15);
                    Visual_range(x1-1,y1,3)=cellstage;
                    Visual_range(x1-1,y1,4)=cell_label_1;
                    cell_array(i,1)=cell_array(i,1)-1;
                    cell_array(i,23)=2;
                    cell_array(i,20)=0;
                    break;
                }
                case 3:
                {
                    int cell_label_1=Visual_range(x1,y1,4);
                    int cellstage=Visual_range(x1,y1,3);
                    Visual_range(x1,y1,all)=0;
                    Visual_range(x1-1,y1+1,1)=1;
                    Visual_range(x1-1,y1+1,2)=cell_array(i,15);
                    Visual_range(x1-1,y1+1,3)=cellstage;
                    Visual_range(x1-1,y1+1,4)=cell_label_1;
                    cell_array(i,1)=cell_array(i,1)-1;
                    cell_array(i,5)=cell_array(i,5)+1;
                    cell_array(i,23)=3;
                    cell_array(i,20)=0;
                    break;
                }
                case 4:
                {
                    int cell_label_1=Visual_range(x1,y1,4);
                    int cellstage=Visual_range(x1,y1,3);
                    Visual_range(x1,y1,all)=0;
                    Visual_range(x1,y1+1,1)=1;
                    Visual_range(x1,y1+1,2)=cell_array(i,15);
                    Visual_range(x1,y1+1,3)=cellstage;
                    Visual_range(x1,y1+1,4)=cell_label_1;
                    cell_array(i,5)=cell_array(i,5)+1;
                    cell_array(i,23)=4;
                    cell_array(i,20)=0;
                    break;
                }
                case 5:
                {
                    int cell_label_1=Visual_range(x1,y1,4);
                    int cellstage=Visual_range(x1,y1,3);
                    Visual_range(x1,y1,all)=0;
                    Visual_range(x1+1,y1+1,1)=1;
                    Visual_range(x1+1,y1+1,2)=cell_array(i,15);
                    Visual_range(x1+1,y1+1,3)=cellstage;
                    Visual_range(x1+1,y1+1,4)=cell_label_1;
                    cell_array(i,1)=cell_array(i,1)+1;
                    cell_array(i,5)=cell_array(i,5)+1;
                    cell_array(i,23)=5;
                    cell_array(i,20)=0;
                    break;
                }
                case 6:
                {
                    int cell_label_1=Visual_range(x1,y1,4);
                    int cellstage=Visual_range(x1,y1,3);
                    Visual_range(x1,y1,all)=0;
                    Visual_range(x1+1,y1,1)=1;
                    Visual_range(x1+1,y1,2)=cell_array(i,15);
                    Visual_range(x1+1,y1,3)=cellstage;
                    Visual_range(x1+1,y1,4)=cell_label_1;
                    cell_array(i,1)=cell_array(i,1)+1;
                    cell_array(i,23)=6;
                    cell_array(i,20)=0;
                    break;
                }
                case 7:
                {
                    int cell_label_1=Visual_range(x1,y1,4);
                    int cellstage=Visual_range(x1,y1,3);
                    Visual_range(x1,y1,all)=0;
                    Visual_range(x1+1,y1-1,1)=1;
                    Visual_range(x1+1,y1-1,2)=cell_array(i,15);
                    Visual_range(x1+1,y1-1,3)=cellstage;
                    Visual_range(x1+1,y1-1,4)=cell_label_1;
                    cell_array(i,1)=cell_array(i,1)+1;
                    cell_array(i,5)=cell_array(i,5)-1;
                    cell_array(i,23)=7;
                    cell_array(i,20)=0;
                    break;
                }
                case 8:
                {
                    int cell_label_1=Visual_range(x1,y1,4);
                    int cellstage=Visual_range(x1,y1,3);
                    Visual_range(x1,y1,all)=0;
                    Visual_range(x1,y1-1,1)=1;
                    Visual_range(x1,y1-1,2)=cell_array(i,15);
                    Visual_range(x1,y1-1,3)=cellstage;
                    Visual_range(x1,y1-1,4)=cell_label_1;
                    cell_array(i,5)=cell_array(i,5)-1;
                    cell_array(i,23)=8;
                    cell_array(i,20)=0;
                    break;
                }
                    
            }
            delete[] direction1;
            direction1 = NULL;
            cell_array(i,20)=0;
            cell_array(i,24)=1;
        }
    }
    migration_judgement=migration_judgement+0.0001;
    gsl_rng_free(r5);
}

#endif /* random_migration_hpp */
