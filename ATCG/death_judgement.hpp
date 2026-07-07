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
#include "stateless_rng.hpp"
using namespace blitz;

template <typename CellArray>
inline bool death_judgement_live_cell(const CellArray &cell_array, int row)
{
    return cell_array(row,cell_col::kViability)==1 && cell_array(row,cell_col::kX1)!=0 && cell_array(row,cell_col::kY1)!=0;
}

template <typename CellArray>
inline void clear_dead_cell_visual(int site, CellArray &cell_array, Array<long,3> &Visual_range, int C)
{
    Range all = Range::all();
    int stage = (int)cell_array(site,cell_col::kStage);
    if(stage==0)
    {
        Visual_range((int)cell_array(site,cell_col::kX1),(int)cell_array(site,cell_col::kY1),all)=0;
        Visual_range((int)cell_array(site,cell_col::kX2),(int)cell_array(site,cell_col::kY2),all)=0;
        Visual_range((int)cell_array(site,cell_col::kX3),(int)cell_array(site,cell_col::kY3),all)=0;
        Visual_range((int)cell_array(site,cell_col::kX4),(int)cell_array(site,cell_col::kY4),all)=0;
    }
    else if(stage==1)
    {
        Visual_range((int)cell_array(site,cell_col::kX1),(int)cell_array(site,cell_col::kY1),all)=0;
    }
    else if(stage==2)
    {
        double usx=cell_array(site,cell_col::kX1);
        double usy=cell_array(site,cell_col::kY1);
        cell_array(site,cell_col::kX1)=1;
        cell_array(site,cell_col::kY1)=1;
        for (int us=1;us<=C;us++)
        {
            if(cell_array(us,cell_col::kX1)!=0 && cell_array(us,cell_col::kY1)!=0 && cell_array(site,cell_col::kStage)==2)
            {
                if(cell_array(us,cell_col::kX1)==usx && cell_array(us,cell_col::kY1)==usy)
                {
                    Visual_range((int)cell_array(us,cell_col::kX1),(int)cell_array(us,cell_col::kY1),3)=1;
                    cell_array(us,cell_col::kStage)=1;
                }
            }
        }
    }
}

inline void compact_after_death_judgement(Array<double, 2> &cell_array, Array<long,3> &Visual_range, int Col, int nthreads, int C)
{
    Range all = Range::all();
    int current_size=cell_array.rows();
    omp_set_num_threads(nthreads);
    int sum =0;
    #pragma omp parallel for schedule(dynamic) reduction(+:sum)
    {
        for (int CN=1; CN<=current_size; ++CN)
        {
            if (death_judgement_live_cell(cell_array, CN))
            {
                sum=sum+1;
            }
        }
    }
    Array<double,2> cell_array_temp(sum,Col,FortranArray<2>());
    cell_array_temp=0;
    int site1=1;
    for (int site=1; site<= current_size; ++site)
    {
        if (death_judgement_live_cell(cell_array, site))
        {
            cell_array_temp(site1,all)=cell_array(site,all);
            site1++;
        }
        else
        {
            clear_dead_cell_visual(site, cell_array, Visual_range, C);
        }
    }
    cell_array.resize(sum,Col);
    cell_array=0;
    cell_array(all,all)=cell_array_temp(all,all);
}

inline void compact_after_death_judgement(CellStore &cells, Array<long,3> &Visual_range, int, int nthreads, int C)
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
inline void death_judgement(int Visual_range_x, int Visual_range_y, int N00, int N01, double r_limit, double K_limit, double lambda_r, double lambda_K, double alpha, double beta, double carrying_capacity_r, double carrying_capacity_K, double Cr, double CK, double death_time_range_r, double death_time_range_K, double deltah, double &h, CellArray &cell_array, Array<long, 3> &sub_visual, Array<long,3> &Visual_range, double deathjudge, int Col,int nthreads,long rng_time_step)
{
    Range all = Range::all();
    int C= cell_array.rows();
//    omp_set_num_threads(nthreads);
//    #pragma omp parallel for schedule(dynamic)
//    {
        for (int rows=1; rows<=C; ++rows)
        {
            long cell_rng_id = (long)cell_array(rows,15);
            if (cell_rng_id == 0)
            {
                cell_rng_id = rows;
            }
            long rng_event = 100;
            if (cell_array(rows,11)>deathjudge)
            {
                if (cell_array(rows,9)==1)
                {
                    if (cell_array(rows,1)>=100 && cell_array(rows,5) >=100 && cell_array(rows,1)<=Visual_range_x+100 && cell_array(rows,5)<=Visual_range_y+100)
                    {
                        sub_visual.resize(6,6,4);
                        sub_visual=0;
                        sub_visual(all,all,all)=Visual_range(Range(cell_array(rows,1)-2,cell_array(rows,1)+3),Range(cell_array(rows,5)-2,cell_array(rows,5)+3),all);
                        long cell_count[36];
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
                        long subcell_r[36]={0};
                        long subcell_K[36]={0};
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
                        double growth_rate_inherent_r=cell_array(rows,10);
                        if (cells_number>=r_limit)
                        {
                            cell_array(rows,11)=growth_rate_inherent_r-((growth_rate_inherent_r*2*(rc+kc+alpha*kc-r_limit))/carrying_capacity_r);
                        }
                        else
                        {
                            cell_array(rows,11)=growth_rate_inherent_r;
                        }
                        if (cell_array(rows,11)>deathjudge)
                        {
                            double expected_division_time=24/cell_array(rows,11);
                            double undividing_time=0.9*expected_division_time;
                            double diving_time_range=0.1*expected_division_time;
                            double probability_of_division=1/diving_time_range;
                            double expected_dividing_time=undividing_time+stateless_geometric(cell_rng_id, rng_time_step, rng_event++, probability_of_division);;
                            cell_array(rows,17)=expected_dividing_time;
                        }
                        else
                        {
                            if (cell_array(rows,18)==0)
                            {
                                double probability_to_death=1/death_time_range_r;
                                cell_array(rows,18)=stateless_geometric(cell_rng_id, rng_time_step, rng_event++, probability_to_death);
                            }
                            cell_array(rows,19)=cell_array(rows,19)+deltah;
                            cell_array(rows,17)=0;
                        }
                        cell_array(rows,21)=1/cell_array(rows,28);
                    }
                    cell_array(rows,21)=1/cell_array(rows,28);
                }
                else
                {
                    if (cell_array(rows,1)>=100 && cell_array(rows,5) >=100 && cell_array(rows,1)<=Visual_range_x+100 && cell_array(rows,5)<=Visual_range_y+100)
                    {
                        sub_visual.resize(6,6,4);
                        sub_visual=0;
                        sub_visual(all,all,all)=Visual_range(Range(cell_array(rows,1)-2,cell_array(rows,1)+3),Range(cell_array(rows,5)-2,cell_array(rows,5)+3),all);
                        long cell_count[36];
                        int cc=0;
                        for (int cx=0; cx<6; cx++)
                        {
                            for(int cy=0; cy<6; cy++)
                            {
                                cell_count[cc]=sub_visual(cx+1,cy+1,4);
                                cc++;
                            }
                        }
                        long cell_stage[36];
                        int cs=0;
                        for (int csx=0; csx<6; csx++)
                        {
                            for(int csy=0; csy<6; csy++)
                            {
                                if(sub_visual(csx+1,csy+1,3)==2)
                                {
                                    cell_stage[cs]=sub_visual(csx+1,csy+1,3);
                                    cs++;
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
                        long subcell_r[36]={0};
                        long subcell_K[36]={0};
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
                        double growth_rate_inherent_K=cell_array(rows,10);
                        if (cells_number>=K_limit)
                        {
                            cell_array(rows,11)=growth_rate_inherent_K-((growth_rate_inherent_K*2*(beta*rc+rc+kc-K_limit))/carrying_capacity_K);
                        }
                        else
                        {
                            cell_array(rows,11)=growth_rate_inherent_K;
                        }
                        if (cell_array(rows,11)>deathjudge)
                        {
                            double expected_division_time=24/cell_array(rows,11);
                            double undividing_time=0.9*expected_division_time;
                            double diving_time_range=0.1*expected_division_time;
                            double probability_of_division=1/diving_time_range;
                            double expected_dividing_time=undividing_time+stateless_geometric(cell_rng_id, rng_time_step, rng_event++, probability_of_division);;
                            cell_array(rows,17)=expected_dividing_time;
                        }
                        else
                        {
                            if (cell_array(rows,18)==0)
                            {
                                double probability_to_death=1/death_time_range_K;
                                cell_array(rows,18)=stateless_geometric(cell_rng_id, rng_time_step, rng_event++, probability_to_death);
                            }
                            cell_array(rows,19)=cell_array(rows,19)+deltah;
                            cell_array(rows,17)=0;
                        }
                        cell_array(rows,21)=1/cell_array(rows,28);
                    }
                    cell_array(rows,21)=1/cell_array(rows,28);
                }
            }
            else if (cell_array(rows,11)<=deathjudge)
            {
                if(h==0)
                {
                    double probability_to_death=0;
                    if (cell_array(rows,9)==1)
                    {
                        probability_to_death=1/death_time_range_r;
                    }
                    else if (cell_array(rows,9)==2)
                    {
                        probability_to_death=1/death_time_range_K;
                    }
                    cell_array(rows,18)=stateless_geometric(cell_rng_id, rng_time_step, rng_event++, probability_to_death);
                }
                else
                {
                    double D_time_1=1.5*(24/cell_array(rows,10));
                    double D_time_2=0.9*cell_array(rows,18);
                    double D_time = 0;
                    if (D_time_1<=D_time_2)
                    {
                        D_time = D_time_1;
                    }
                    else
                    {
                        D_time = D_time_2;
                    }
                    if (cell_array(rows,19)<=D_time)
                    {
                        if (cell_array(rows,9)==1)
                        {
                            if (cell_array(rows,1)>=100 && cell_array(rows,5) >=100 && cell_array(rows,1)<=Visual_range_x+100 && cell_array(rows,5)<=Visual_range_y+100)
                            {
                                sub_visual.resize(6,6,4);
                                sub_visual=0;
                                sub_visual(all,all,all)=Visual_range(Range(cell_array(rows,1)-2,cell_array(rows,1)+3),Range(cell_array(rows,5)-2,cell_array(rows,5)+3),all);
                                long cell_count[36];
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
                                long subcell_r[36]={0};
                                long subcell_K[36]={0};
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
                                double growth_rate_inherent_r=cell_array(rows,10);
                                if (cells_number>=r_limit)
                                {
                                    cell_array(rows,11)=growth_rate_inherent_r-((growth_rate_inherent_r*2*(rc+kc+alpha*kc-r_limit))/carrying_capacity_r);
                                }
                                else
                                {
                                    cell_array(rows,11)=growth_rate_inherent_r;
                                }
                                if (cell_array(rows,11)>deathjudge)
                                {
                                    double expected_division_time=24/cell_array(rows,11);
                                    double undividing_time=0.9*expected_division_time;
                                    double diving_time_range=0.1*expected_division_time;
                                    double probability_of_division=1/diving_time_range;
                                    double expected_dividing_time=undividing_time+stateless_geometric(cell_rng_id, rng_time_step, rng_event++, probability_of_division);;
                                    cell_array(rows,17)=expected_dividing_time;
                                    cell_array(rows,18)=0;
                                    cell_array(rows,19)=0;
                                }
                                else
                                {
                                    cell_array(rows,19)=cell_array(rows,19)+deltah;
                                    cell_array(rows,17)=0;
                                }
                                cell_array(rows,21)=1/cell_array(rows,28);
                            }
                            cell_array(rows,21)=1/cell_array(rows,28);
                        }
                        else
                        {
                            if (cell_array(rows,1)>=100 && cell_array(rows,5) >=100 && cell_array(rows,1)<=Visual_range_x+100 && cell_array(rows,5)<=Visual_range_y+100)
                            {
                                sub_visual.resize(6,6,4);
                                sub_visual=0;
                                sub_visual(all,all,all)=Visual_range(Range(cell_array(rows,1)-2,cell_array(rows,1)+3),Range(cell_array(rows,5)-2,cell_array(rows,5)+3),all);
                                long cell_count[36];
                                int cc=0;
                                for (int cx=0; cx<6; cx++)
                                {
                                    for(int cy=0; cy<6; cy++)
                                    {
                                        cell_count[cc]=sub_visual(cx+1,cy+1,4);
                                        cc++;
                                    }
                                }
                                long cell_stage[36];
                                int cs=0;
                                for (int csx=0; csx<6; csx++)
                                {
                                    for(int csy=0; csy<6; csy++)
                                    {
                                        if(sub_visual(csx+1,csy+1,3)==2)
                                        {
                                            cell_stage[cs]=sub_visual(csx+1,csy+1,3);
                                            cs++;
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
                                long subcell_r[36]={0};
                                long subcell_K[36]={0};
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
                                double growth_rate_inherent_K=cell_array(rows,10);
                                if (cells_number>=K_limit)
                                {
                                    cell_array(rows,11)=growth_rate_inherent_K-((growth_rate_inherent_K*2*(beta*rc+rc+kc-K_limit))/carrying_capacity_K);
                                }
                                else
                                {
                                    cell_array(rows,11)=growth_rate_inherent_K;
                                }
                                if (cell_array(rows,11)>deathjudge)
                                {
                                    double expected_division_time=24/cell_array(rows,11);
                                    double undividing_time=0.9*expected_division_time;
                                    double diving_time_range=0.1*expected_division_time;
                                    double probability_of_division=1/diving_time_range;
                                    double expected_dividing_time=undividing_time+stateless_geometric(cell_rng_id, rng_time_step, rng_event++, probability_of_division);;
                                    cell_array(rows,17)=expected_dividing_time;
                                    cell_array(rows,18)=0;
                                    cell_array(rows,19)=0;
                                }
                                else
                                {
                                    cell_array(rows,19)=cell_array(rows,19)+deltah;
                                    cell_array(rows,17)=0;
                                }
                                cell_array(rows,21)=1/cell_array(rows,28);
                            }
                            cell_array(rows,21)=1/cell_array(rows,28);
                        }
                    }
                    else
                    {
                        if(cell_array(rows,18)<=cell_array(rows,19))
                        {
                            cell_array(rows,22)=0;
                        }
                        else
                        {
                            cell_array(rows,19)=cell_array(rows,19)+deltah;
                        }
                    }
                }
            }
        }
//    }
    compact_after_death_judgement(cell_array, Visual_range, Col, nthreads, C);
}
#endif /* death_judgement_hpp */
