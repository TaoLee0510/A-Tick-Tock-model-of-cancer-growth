//
//  density_growth_rate_calculation_1.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef density_growth_rate_calculation_1_hpp
#define density_growth_rate_calculation_1_hpp

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include <blitz/blitz.h>
#include <blitz/array.h>

using namespace blitz;
void density_growth_rate_calculation_1(int Visual_range_x, int Visual_range_y, int N00, int N01, double r_limit, double K_limit, double lambda_r, double lambda_K, float alpha, float beta, double carrying_capacity_r, double carrying_capacity_K, double Cr, double CK, double death_time_range_r, double death_time_range_K, Array<float, 2> &cell_array,Array<int, 3> &sub_visual,Array<int,3> Visual_range)
{
    Range all = Range::all();
    int C0= cell_array.rows();
    const gsl_rng_type *T3;
    gsl_rng *r3;
    gsl_rng_env_setup();
    T3 = gsl_rng_ranlxs0;
    gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    r3 = gsl_rng_alloc(T3);
    for (int i=1; i<=C0; i++)
    {
        sub_visual.resize(6,6,4);
        sub_visual=0;
        int cell_type=(int)cell_array(i,9);
        switch (cell_type)
        {
            case 1:
            {
                if (cell_array(i,1)>=100 && cell_array(i,5) >=100 && cell_array(i,1)<=Visual_range_x+100 && cell_array(i,5)<=Visual_range_y+100)
                {
                    sub_visual(all,all,all)=Visual_range(Range(cell_array(i,1)-2,cell_array(i,1)+3),Range(cell_array(i,5)-2,cell_array(i,5)+3),all);
                    int cell_count[36];
                    int cc=0;
                    for (int cx=0; cx<6; cx++)
                    {
                        for(int cy=0; cy<6; cy++)
                        {
                            cell_count[cc]=sub_visual(cx+1,cy+1,4);
                            cc++;
                        }
                    }
                    vector<int> mycellcount (cell_count, cell_count+36);
                    sort(mycellcount.begin(),mycellcount.end());
                    mycellcount.erase(unique(mycellcount.begin(), mycellcount.end()), mycellcount.end());
                    long cells_number=0;
                    cells_number = mycellcount.size();
                    if (mycellcount[0]==0)
                    {
                        cells_number=cells_number-1;
                    }
                    int subcell_r[36]={0};
                    int subcell_K[36]={0};
                    int r_cells_number=0;
                    int K_cells_number=0;
                    for(int cxx=1;cxx<=6;cxx++)
                    {
                        for(int cyy=1;cyy<=6;cyy++)
                        {
                            if((sub_visual(cxx,cyy,2)<=N00/2 && sub_visual(cxx,cyy,2)>=1)||(sub_visual(cxx,cyy,2)<=N00+(N01/2) && sub_visual(cxx,cyy,2) >N00))
                            {
                                subcell_r[r_cells_number]=sub_visual(cxx,cyy,4);
                                r_cells_number=r_cells_number+1;
                            }
                            else
                            {
                                subcell_K[K_cells_number]=sub_visual(cxx,cyy,4);
                                K_cells_number=K_cells_number+1;
                                
                            }
                        }
                    }
                    vector<int> myrcellcount (subcell_r, subcell_r+36);
                    sort(myrcellcount.begin(),myrcellcount.end());
                    myrcellcount.erase(unique(myrcellcount.begin(), myrcellcount.end()), myrcellcount.end());
                    vector<int> myKcellcount (subcell_K, subcell_K+36);
                    sort(myKcellcount.begin(),myKcellcount.end());
                    myKcellcount.erase(unique(myKcellcount.begin(), myKcellcount.end()), myKcellcount.end());
                    long rc=myrcellcount.size();
                    long kc1=myKcellcount.size();
                    if (myrcellcount[0]==0)
                    {
                        rc=rc-1;
                    }
                    if (myKcellcount[0]==0)
                    {
                        kc1=kc1-1;
                    }
                    long kc_usmall=0;
                    for (int cx=0; cx<6; cx++)
                    {
                        for(int cy=0; cy<6; cy++)
                        {
                            if (sub_visual(cx+1,cy+1,3)==2)
                            {
                                kc_usmall=kc_usmall+1;
                            }
                        }
                    }
                    long kc=kc1+kc_usmall;
                    cells_number=kc+rc;
                    double growth_rate_inherent_r=cell_array(i,10);
                    if (cells_number>=r_limit)
                    {
                        cell_array(i,11)=growth_rate_inherent_r-((growth_rate_inherent_r*2*(rc+kc+alpha*kc-r_limit))/carrying_capacity_r);
                    }
                    if (cell_array(i,11)>0)
                    {
                        double expected_division_time=24/cell_array(i,11);
                        double undividing_time=0.9*expected_division_time;
                        double diving_time_range=0.1*expected_division_time;
                        double probability_of_division=1/diving_time_range;
                        double expected_dividing_time=undividing_time+gsl_ran_geometric(r3, probability_of_division);;
                        cell_array(i,17)=expected_dividing_time;
                    }
                    cell_array(i,21)=1/cell_array(i,28);
                }
                break;
            }
            case 2:
            {
                if (cell_array(i,1)>=100 && cell_array(i,5) >=100 && cell_array(i,1)<=Visual_range_x+100 && cell_array(i,5)<=Visual_range_y+100)
                {
                    sub_visual(all,all,all)=Visual_range(Range(cell_array(i,1)-2,cell_array(i,1)+3),Range(cell_array(i,5)-2,cell_array(i,5)+3),all);
                    int cell_count[36];
                    int cc=0;
                    for (int cx=0; cx<6; cx++)
                    {
                        for(int cy=0; cy<6; cy++)
                        {
                            cell_count[cc]=sub_visual(cx+1,cy+1,4);
                            cc++;
                        }
                    }
                    int cs=0;
                    for (int csx=0; csx<6; csx++)
                    {
                        for(int csy=0; csy<6; csy++)
                        {
                            if(sub_visual(csx+1,csy+1,3)==2)
                            {
                                cs=cs+1;
                            }
                        }
                    }
                    int stage_sum=cs;
                    vector<int> mycellcount (cell_count, cell_count+36);
                    sort(mycellcount.begin(),mycellcount.end());
                    mycellcount.erase(unique(mycellcount.begin(), mycellcount.end()), mycellcount.end());
                    long cells_number=0;
                    cells_number = mycellcount.size();
                    if (mycellcount[0]==0)
                    {
                        cells_number=cells_number-1;
                    }
                    int subcell_r[36]={0};
                    int subcell_K[36]={0};
                    int r_cells_number=0;
                    int K_cells_number=0;
                    for(int cxx=1;cxx<=6;cxx++)
                    {
                        for(int cyy=1;cyy<=6;cyy++)
                        {
                            if((sub_visual(cxx,cyy,2)<=N00/2 && sub_visual(cxx,cyy,2)>=1)||(sub_visual(cxx,cyy,2)<=N00+(N01/2) && sub_visual(cxx,cyy,2) >N00))
                            {
                                subcell_r[r_cells_number]=sub_visual(cxx,cyy,4);
                                r_cells_number=r_cells_number+1;
                            }
                            else
                            {
                                subcell_K[K_cells_number]=sub_visual(cxx,cyy,4);
                                K_cells_number=K_cells_number+1;
                                
                            }
                        }
                    }
                    vector<int> myrcellcount (subcell_r, subcell_r+36);
                    sort(myrcellcount.begin(),myrcellcount.end());
                    myrcellcount.erase(unique(myrcellcount.begin(), myrcellcount.end()), myrcellcount.end());
                    vector<int> myKcellcount (subcell_K, subcell_K+36);
                    sort(myKcellcount.begin(),myKcellcount.end());
                    myKcellcount.erase(unique(myKcellcount.begin(), myKcellcount.end()), myKcellcount.end());
                    long rc=myrcellcount.size();
                    long kc=myKcellcount.size();
                    if (myrcellcount[0]==0)
                    {
                        rc=rc-1;
                    }
                    if (myKcellcount[0]==0)
                    {
                        kc=kc-1;
                    }
                    kc=kc+stage_sum;
                    cells_number=kc+rc;
                    double growth_rate_inherent_K=cell_array(i,10);
                    if (cells_number>=K_limit)
                    {
                        cell_array(i,11)=growth_rate_inherent_K-((growth_rate_inherent_K*2*(beta*rc+rc+kc-K_limit))/carrying_capacity_K);
                    }
                    if (cell_array(i,11)>0)
                    {
                        double expected_division_time=24/cell_array(i,11);
                        double undividing_time=0.9*expected_division_time;
                        double diving_time_range=0.1*expected_division_time;
                        double probability_of_division=1/diving_time_range;
                        double expected_dividing_time=undividing_time+gsl_ran_geometric(r3, probability_of_division);;
                        cell_array(i,17)=expected_dividing_time;
                    }
                    cell_array(i,21)=1/cell_array(i,28);
                }
                break;
            }
        }
    }
    gsl_rng_free(r3);
}
#endif /* density_growth_rate_calculation_1_hpp */
