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
//#include <omp.h>
#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include <blitz/blitz.h>
#include <blitz/array.h>
using namespace blitz;
void death_judgement(int Visual_range_x, int Visual_range_y, int N00, int N01, double r_limit, double K_limit, double lambda_r, double lambda_K, double alpha, double beta, double carrying_capacity_r, double carrying_capacity_K, double Cr, double CK, double death_time_range_r, double death_time_range_K, double deltah, double &h, Array<double, 2> &cell_array, Array<double,2> &cell_array_temp, Array<long, 3> sub_visual, Array<long,3> &Visual_range, double deathjudge, int Col)
{
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 RNG(seed);
    Range all = Range::all();
    int C= cell_array.rows();
    const gsl_rng_type *T4;
    gsl_rng *r4;
    gsl_rng_env_setup();
    T4 = gsl_rng_ranlxs0;
    gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    r4 = gsl_rng_alloc(T4);
    for (int rows=1; rows<=C; rows++)
    {
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
                        double expected_dividing_time=undividing_time+gsl_ran_geometric(r4, probability_of_division);;
                        cell_array(rows,17)=expected_dividing_time;
                    }
                    else
                    {
                        if (cell_array(rows,18)==0)
                        {
                            double probability_to_death=1/death_time_range_r;
                            cell_array(rows,18)=gsl_ran_geometric(r4,probability_to_death);
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
                        double expected_dividing_time=undividing_time+gsl_ran_geometric(r4, probability_of_division);;
                        cell_array(rows,17)=expected_dividing_time;
                    }
                    else
                    {
                        if (cell_array(rows,18)==0)
                        {
                            double probability_to_death=1/death_time_range_K;
                            cell_array(rows,18)=gsl_ran_geometric(r4,probability_to_death);
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
                cell_array(rows,18)=gsl_ran_geometric(r4,probability_to_death);
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
                                double expected_dividing_time=undividing_time+gsl_ran_geometric(r4, probability_of_division);;
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
                                double expected_dividing_time=undividing_time+gsl_ran_geometric(r4, probability_of_division);;
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
    int current_size=cell_array.rows();
    int sum =0;
    for (int CN=1; CN<=current_size; CN++)
    {
        if (cell_array(CN,22)==1 && cell_array(CN,1)!=0 && cell_array(CN,5)!=0)
        {
            sum=sum+1;
        }
    }
    cell_array_temp.resize(sum,Col);
    cell_array_temp=0;
    int site1=1;
    for (int site=1; site<= current_size; site++)
    {
        if (cell_array(site,22)==1 && cell_array(site,1)!=0 && cell_array(site,5)!=0)
        {
            cell_array_temp(site1,all)=cell_array(site,all);
            site1++;
        }
        else
        {
            if(cell_array(site,14)==0)
            {
                Visual_range(cell_array(site,1),cell_array(site,5),all)=0;
                Visual_range(cell_array(site,2),cell_array(site,6),all)=0;
                Visual_range(cell_array(site,3),cell_array(site,7),all)=0;
                Visual_range(cell_array(site,4),cell_array(site,8),all)=0;
            }
            else if(cell_array(site,14)==1)
            {
                Visual_range(cell_array(site,1),cell_array(site,5),all)=0;
            }
            else if(cell_array(site,14)==2)
            {
                double usx=cell_array(site,1);
                double usy=cell_array(site,5);
                cell_array(site,1)=1;
                cell_array(site,5)=1;
                for (int us=1;us<=C;us++)
                {
                    if(cell_array(us,1)!=0 && cell_array(us,5)!=0 && cell_array(site,14)==2)
                    {
                        if(cell_array(us,1)==usx && cell_array(us,5)==usy)
                        {
                            Visual_range((int)cell_array(us,1),(int)cell_array(us,5),3)=1;
                            cell_array(us,14)=1;
                        }
                    }
                }
            }
        }
    }
    cell_array.resize(sum,Col);
    cell_array=0;
    cell_array(all,all)=cell_array_temp(all,all);
    gsl_rng_free(r4);
}
#endif /* death_judgement_hpp */
