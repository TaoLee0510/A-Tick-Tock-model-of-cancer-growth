//
//  migration.hpp
//  CCSCIM
//
//  Created by Tao Lee on 5/11/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef migration_hpp
#define migration_hpp

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
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "deltah_calculation.hpp"

using namespace blitz;
void migration(int i, double deltah, Array<float, 2> &cell_array, Array<int, 3> &Visual_range, Array<int,2> cor_big, Array<int, 2> area_square, Array<int, 2> sub_area_square, Array<int, 2> cor_small, Array<int, 2> area_square_s, Array<int, 2>  sub_area_square_s,double &migration_judgement)
{
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 RNG(seed);
    Range all = Range::all();
    const gsl_rng_type *T6;
    gsl_rng *r6;
    gsl_rng_env_setup();
    T6 = gsl_rng_ranlxs0;
    gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    r6 = gsl_rng_alloc(T6);
    int x1=cell_array(i,1);
    int y1=cell_array(i,5);
    if (cell_array(i,9)==1)
    {
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
                double density[8]={0};
                area_square.resize(10,10);
                area_square=0;
                area_square(all,all)=Visual_range(Range(x1-4,x1+5),Range(y1-4,y1+5),4);
                //////////////////////////////////////////////8 direction density calculation////////////////////////////////////
                for (int loci_for_mig=0; loci_for_mig<8;loci_for_mig++)
                {
                    if (loci_for_mig==0)
                    {
                        sub_area_square.resize(10,10);
                        sub_area_square=0;
                        sub_area_square(Range(1,5),Range(1,5))=area_square(Range(1,5),Range(1,5));
                        int cell_count[100]={0};
                        int cc=0;
                        for (int cx=0; cx<10; cx++)
                        {
                            for(int cy=0; cy<10; cy++)
                            {
                                cell_count[cc]=sub_area_square(cx+1,cy+1);
                                cc++;
                            }
                        }
                        vector<int> mycellcount (cell_count, cell_count+100);
                        sort(mycellcount.begin(),mycellcount.end());
                        mycellcount.erase(unique(mycellcount.begin(), mycellcount.end()), mycellcount.end());
                        long cells_number=0;
                        cells_number = mycellcount.size();
                        if (mycellcount[0]==0)
                        {
                            cells_number=cells_number-1;
                        }
                        density[0]=(double)cells_number/(double)25;
                    }
                    else if (loci_for_mig==2)
                    {
                        sub_area_square.resize(10,10);
                        sub_area_square=0;
                        sub_area_square(Range(1,5),Range(6,10))=area_square(Range(1,5),Range(6,10));
                        int cell_count[100]={0};
                        int cc=0;
                        for (int cx=0; cx<10; cx++)
                        {
                            for(int cy=0; cy<10; cy++)
                            {
                                cell_count[cc]=sub_area_square(cx+1,cy+1);
                                cc++;
                            }
                        }
                        vector<int> mycellcount (cell_count, cell_count+100);
                        sort(mycellcount.begin(),mycellcount.end());
                        mycellcount.erase(unique(mycellcount.begin(), mycellcount.end()), mycellcount.end());
                        long cells_number=0;
                        cells_number = mycellcount.size();
                        if (mycellcount[0]==0)
                        {
                            cells_number=cells_number-1;
                        }
                        density[2]=(double)cells_number/(double)25;
                    }
                    else if (loci_for_mig==4)
                    {
                        sub_area_square.resize(10,10);
                        sub_area_square=0;
                        sub_area_square(Range(6,10),Range(6,10))=area_square(Range(6,10),Range(6,10));
                        int cell_count[100]={0};
                        int cc=0;
                        for (int cx=0; cx<10; cx++)
                        {
                            for(int cy=0; cy<10; cy++)
                            {
                                cell_count[cc]=sub_area_square(cx+1,cy+1);
                                cc++;
                            }
                        }
                        vector<int> mycellcount (cell_count, cell_count+100);
                        sort(mycellcount.begin(),mycellcount.end());
                        mycellcount.erase(unique(mycellcount.begin(), mycellcount.end()), mycellcount.end());
                        long cells_number=0;
                        cells_number = mycellcount.size();
                        if (mycellcount[0]==0)
                        {
                            cells_number=cells_number-1;
                        }
                        density[4]=(double)cells_number/(double)25;
                    }
                    else if (loci_for_mig==6)
                    {
                        sub_area_square.resize(10,10);
                        sub_area_square=0;
                        sub_area_square(Range(6,10),Range(1,5))=area_square(Range(6,10),Range(1,5));
                        int cell_count[100]={0};
                        int cc=0;
                        for (int cx=0; cx<10; cx++)
                        {
                            for(int cy=0; cy<10; cy++)
                            {
                                cell_count[cc]=sub_area_square(cx+1,cy+1);
                                cc++;
                            }
                        }
                        vector<int> mycellcount (cell_count, cell_count+100);
                        sort(mycellcount.begin(),mycellcount.end());
                        mycellcount.erase(unique(mycellcount.begin(), mycellcount.end()), mycellcount.end());
                        long cells_number=0;
                        cells_number = mycellcount.size();
                        if (mycellcount[0]==0)
                        {
                            cells_number=cells_number-1;
                        }
                        density[6]=(double)cells_number/(double)25;
                    }
                    else if (loci_for_mig==1)
                    {
                        sub_area_square.resize(10,10);
                        sub_area_square=0;
                        int deltay=0;
                        for (int xss=1; xss<=5; xss++)
                        {
                            for (int yss=1+deltay; yss<=10-deltay; yss++)
                            {
                                sub_area_square(xss,yss)=area_square(xss,yss);
                            }
                            deltay++;
                        }
                        int cell_count[100]={0};
                        int cc=0;
                        for (int cx=0; cx<10; cx++)
                        {
                            for(int cy=0; cy<10; cy++)
                            {
                                cell_count[cc]=sub_area_square(cx+1,cy+1);
                                cc++;
                            }
                        }
                        vector<int> mycellcount (cell_count, cell_count+100);
                        sort(mycellcount.begin(),mycellcount.end());
                        mycellcount.erase(unique(mycellcount.begin(), mycellcount.end()), mycellcount.end());
                        long cells_number=0;
                        cells_number = mycellcount.size();
                        if (mycellcount[0]==0)
                        {
                            cells_number=cells_number-1;
                        }
                        density[1]=(double)cells_number/(double)30;
                    }
                    else if (loci_for_mig==3)
                    {
                        sub_area_square.resize(10,10);
                        sub_area_square=0;
                        int deltay=0;
                        for (int yss=10; yss>=6; yss--)
                        {
                            for (int xss=10-deltay; xss>=1+deltay; xss--)
                            {
                                sub_area_square(xss,yss)=area_square(xss,yss);
                            }
                            deltay++;
                        }
                        int cell_count[100]={0};
                        int cc=0;
                        for (int cx=0; cx<10; cx++)
                        {
                            for(int cy=0; cy<10; cy++)
                            {
                                cell_count[cc]=sub_area_square(cx+1,cy+1);
                                cc++;
                            }
                        }
                        vector<int> mycellcount (cell_count, cell_count+100);
                        sort(mycellcount.begin(),mycellcount.end());
                        mycellcount.erase(unique(mycellcount.begin(), mycellcount.end()), mycellcount.end());
                        long cells_number=0;
                        cells_number = mycellcount.size();
                        if (mycellcount[0]==0)
                        {
                            cells_number=cells_number-1;
                        }
                        density[3]=(double)cells_number/(double)30;
                    }
                    else if (loci_for_mig==5)
                    {
                        sub_area_square.resize(10,10);
                        sub_area_square=0;
                        int deltay=0;
                        for (int xss=10; xss>=6; xss--)
                        {
                            for (int yss=1+deltay; yss<=10-deltay; yss++)
                            {
                                sub_area_square(xss,yss)=area_square(xss,yss);
                            }
                            deltay++;
                        }
                        int cell_count[100]={0};
                        int cc=0;
                        for (int cx=0; cx<10; cx++)
                        {
                            for(int cy=0; cy<10; cy++)
                            {
                                cell_count[cc]=sub_area_square(cx+1,cy+1);
                                cc++;
                            }
                        }
                        vector<int> mycellcount (cell_count, cell_count+100);
                        sort(mycellcount.begin(),mycellcount.end());
                        mycellcount.erase(unique(mycellcount.begin(), mycellcount.end()), mycellcount.end());
                        long cells_number=0;
                        cells_number = mycellcount.size();
                        if (mycellcount[0]==0)
                        {
                            cells_number=cells_number-1;
                        }
                        density[5]=(double)cells_number/(double)30;
                    }
                    else if (loci_for_mig==7)
                    {
                        sub_area_square.resize(10,10);
                        sub_area_square=0;
                        int deltay=0;
                        for (int yss=1; yss<=5; yss++)
                        {
                            for (int xss=1+deltay; xss<=10-deltay; xss++)
                            {
                                sub_area_square(xss,yss)=area_square(xss,yss);
                            }
                            deltay++;
                        }
                        int cell_count[100]={0};
                        int cc=0;
                        for (int cx=0; cx<10; cx++)
                        {
                            for(int cy=0; cy<10; cy++)
                            {
                                cell_count[cc]=sub_area_square(cx+1,cy+1);
                                cc++;
                            }
                        }
                        vector<int> mycellcount (cell_count, cell_count+100);
                        sort(mycellcount.begin(),mycellcount.end());
                        mycellcount.erase(unique(mycellcount.begin(), mycellcount.end()), mycellcount.end());
                        long cells_number=0;
                        cells_number = mycellcount.size();
                        if (mycellcount[0]==0)
                        {
                            cells_number=cells_number-1;
                        }
                        density[7]=(double)cells_number/(double)30;
                    }
                }
                int order=0;
                double mean_density=0.6;
                ////////////////////////////////////////////////initial migration direction dudgement////////////////////////////////////////////
                if (cell_array(i,23)==0)
                {
                    int new_direction_number=0;
                    for (int dl=0;dl<8;dl++)
                    {
                        if (density[dl]<=mean_density)
                        {
                            new_direction_number=new_direction_number+1;
                        }
                    }
                    if (new_direction_number>0)
                    {
                        int *new_direction_for_migration=new int[new_direction_number];
                        int new_direction_number_for_migration=0;
                        for (int dl=0;dl<8;dl++)
                        {
                            if (density[dl]<=mean_density)
                            {
                                new_direction_for_migration[new_direction_number_for_migration]=direction[dl];
                                new_direction_number_for_migration=new_direction_number_for_migration+1;
                            }
                        }
                        int new_loci_for_migratio_number=0;
                        for (int density_nonzero=0; density_nonzero<new_direction_number_for_migration; density_nonzero++)
                        {
                            if (new_direction_for_migration[density_nonzero]>0)
                            {
                                new_loci_for_migratio_number=new_loci_for_migratio_number+1;
                            }
                        }
                        int *new_loci_for_migration=new int[new_loci_for_migratio_number];
                        int locinumber=0;
                        for (int loci_mig=0;loci_mig<new_direction_number_for_migration;loci_mig++)
                        {
                            if (new_direction_for_migration[loci_mig]>0)
                            {
                                new_loci_for_migration[locinumber]=new_direction_for_migration[loci_mig];
                                locinumber=locinumber+1;
                            }
                        }
                        if (locinumber>0)
                        {
                            shuffle(new_loci_for_migration,new_loci_for_migration+locinumber,RNG);
                            order=new_loci_for_migration[0];
                        }
                        delete[] new_direction_for_migration;
                        new_direction_for_migration = NULL;
                        delete[] new_loci_for_migration;
                        new_loci_for_migration = NULL;
                        cell_array(i,24)=1;
                    }
                }
                ///////////////////////////////////////////////////////following migration direction dudgement////////////////////////////////////////////
                else if (cell_array(i,23)==1)
                {
                    int pro_loci_left=0;
                    int pro_loci_right=0;
                    int pro_loci_mid=0;
                    for (int x=0; x<new_loci; x++)
                    {
                        if (direction1[x]==2)
                        {
                            pro_loci_left++;
                        }
                        else if (direction1[x]==8)
                        {
                            pro_loci_right++;
                        }
                        else if (direction1[x]==1)
                        {
                            pro_loci_mid++;
                        }
                    }
                    if (density[0]<=mean_density || density[1]<=mean_density || density[7]<=mean_density)
                    {
                        if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                        {
                            int order_loci[3]={8,1,2};
                            int n=3;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            int order_loci[2]={2,8};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                        {
                            int order_loci[2]={1,2};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                        {
                            int order_loci[2]={8,1};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            order=2;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                        {
                            order=1;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            order=8;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            shuffle(direction1,direction1+new_loci,RNG);
                            order=direction1[0];
                        }
                    }
                    else
                    {
                        int new_direction_number=0;
                        for (int dl=0;dl<8;dl++)
                        {
                            if (density[dl]<=mean_density)
                            {
                                new_direction_number=new_direction_number+1;
                            }
                        }
                        if (new_direction_number>0)
                        {
                            int *new_direction_for_migration=new int[new_direction_number];
                            int new_direction_number_for_migration=0;
                            for (int dl=0;dl<8;dl++)
                            {
                                if (density[dl]<=mean_density)
                                {
                                    new_direction_for_migration[new_direction_number_for_migration]=direction[dl];
                                    new_direction_number_for_migration=new_direction_number_for_migration+1;
                                }
                            }
                            int new_loci_for_migratio_number=0;
                            for (int density_nonzero=0; density_nonzero<new_direction_number_for_migration; density_nonzero++)
                            {
                                if (new_direction_for_migration[density_nonzero]>0)
                                {
                                    new_loci_for_migratio_number=new_loci_for_migratio_number+1;
                                }
                            }
                            int *new_loci_for_migration=new int[new_loci_for_migratio_number];
                            int locinumber=0;
                            for (int loci_mig=0;loci_mig<new_direction_number_for_migration;loci_mig++)
                            {
                                if (new_direction_for_migration[loci_mig]>0)
                                {
                                    new_loci_for_migration[locinumber]=new_direction_for_migration[loci_mig];
                                    locinumber=locinumber+1;
                                }
                            }
                            if (locinumber>0)
                            {
                                shuffle(new_loci_for_migration,new_loci_for_migration+locinumber,RNG);
                                order=new_loci_for_migration[0];
                            }
                            delete[] new_direction_for_migration;
                            new_direction_for_migration = NULL;
                            delete[] new_loci_for_migration;
                            new_loci_for_migration = NULL;
                        }
                    }
                }
                else if (cell_array(i,23)==8)
                {
                    int pro_loci_left=0;
                    int pro_loci_right=0;
                    int pro_loci_mid=0;
                    for (int x=0; x<new_loci; x++)
                    {
                        if (direction1[x]==7)
                        {
                            pro_loci_left++;
                        }
                        else if (direction1[x]==1)
                        {
                            pro_loci_right++;
                        }
                        else if (direction1[x]==8)
                        {
                            pro_loci_mid++;
                        }
                    }
                    if (density[7]<=mean_density || density[6]<=mean_density || density[0]<=mean_density)
                    {
                        if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                        {
                            int order_loci[3]={7,8,1};
                            int n=3;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            int order_loci[2]={7,1};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                        {
                            int order_loci[2]={7,8};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                        {
                            int order_loci[2]={8,1};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            order=7;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                        {
                            order=8;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            order=1;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            shuffle(direction1,direction1+new_loci,RNG);
                            order=direction1[0];
                        }
                    }
                    else
                    {
                        int new_direction_number=0;
                        for (int dl=0;dl<8;dl++)
                        {
                            if (density[dl]<=mean_density)
                            {
                                new_direction_number=new_direction_number+1;
                            }
                        }
                        if (new_direction_number>0)
                        {
                            int *new_direction_for_migration=new int[new_direction_number];
                            int new_direction_number_for_migration=0;
                            for (int dl=0;dl<8;dl++)
                            {
                                if (density[dl]<=mean_density)
                                {
                                    new_direction_for_migration[new_direction_number_for_migration]=direction[dl];
                                    new_direction_number_for_migration=new_direction_number_for_migration+1;
                                }
                            }
                            int new_loci_for_migratio_number=0;
                            for (int density_nonzero=0; density_nonzero<new_direction_number_for_migration; density_nonzero++)
                            {
                                if (new_direction_for_migration[density_nonzero]>0)
                                {
                                    new_loci_for_migratio_number=new_loci_for_migratio_number+1;
                                }
                            }
                            int *new_loci_for_migration=new int[new_loci_for_migratio_number];
                            int locinumber=0;
                            for (int loci_mig=0;loci_mig<new_direction_number_for_migration;loci_mig++)
                            {
                                if (new_direction_for_migration[loci_mig]>0)
                                {
                                    new_loci_for_migration[locinumber]=new_direction_for_migration[loci_mig];
                                    locinumber=locinumber+1;
                                }
                            }
                            if (locinumber>0)
                            {
                                shuffle(new_loci_for_migration,new_loci_for_migration+locinumber,RNG);
                                order=new_loci_for_migration[0];
                            }
                            delete[] new_direction_for_migration;
                            new_direction_for_migration = NULL;
                            delete[] new_loci_for_migration;
                            new_loci_for_migration = NULL;
                        }
                    }
                }
                else if (cell_array(i,23)==2)
                {
                    int pro_loci_left=0;
                    int pro_loci_right=0;
                    int pro_loci_mid=0;
                    for (int x=0; x<new_loci; x++)
                    {
                        if (direction1[x]==1)
                        {
                            pro_loci_left++;
                        }
                        else if (direction1[x]==3)
                        {
                            pro_loci_right++;
                        }
                        else if (direction1[x]==2)
                        {
                            pro_loci_mid++;
                        }
                    }
                    if (density[1]<=mean_density || density[0]<=mean_density || density[2]<=mean_density)
                    {
                        if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                        {
                            int order_loci[3]={1,2,3};
                            int n=3;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            int order_loci[2]={1,3};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                        {
                            int order_loci[2]={1,2};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                        {
                            int order_loci[2]={2,3};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            order=1;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                        {
                            order=2;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            order=3;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            shuffle(direction1,direction1+new_loci,RNG);
                            order=direction1[0];
                        }
                    }
                    else
                    {
                        int new_direction_number=0;
                        for (int dl=0;dl<8;dl++)
                        {
                            if (density[dl]<=mean_density)
                            {
                                new_direction_number=new_direction_number+1;
                            }
                        }
                        if (new_direction_number>0)
                        {
                            int *new_direction_for_migration=new int[new_direction_number];
                            int new_direction_number_for_migration=0;
                            for (int dl=0;dl<8;dl++)
                            {
                                if (density[dl]<=mean_density)
                                {
                                    new_direction_for_migration[new_direction_number_for_migration]=direction[dl];
                                    new_direction_number_for_migration=new_direction_number_for_migration+1;
                                }
                            }
                            int new_loci_for_migratio_number=0;
                            for (int density_nonzero=0; density_nonzero<new_direction_number_for_migration; density_nonzero++)
                            {
                                if (new_direction_for_migration[density_nonzero]>0)
                                {
                                    new_loci_for_migratio_number=new_loci_for_migratio_number+1;
                                }
                            }
                            int *new_loci_for_migration=new int[new_loci_for_migratio_number];
                            int locinumber=0;
                            for (int loci_mig=0;loci_mig<new_direction_number_for_migration;loci_mig++)
                            {
                                if (new_direction_for_migration[loci_mig]>0)
                                {
                                    new_loci_for_migration[locinumber]=new_direction_for_migration[loci_mig];
                                    locinumber=locinumber+1;
                                }
                            }
                            if (locinumber>0)
                            {
                                shuffle(new_loci_for_migration,new_loci_for_migration+locinumber,RNG);
                                order=new_loci_for_migration[0];
                            }
                            delete[] new_direction_for_migration;
                            new_direction_for_migration = NULL;
                            delete[] new_loci_for_migration;
                            new_loci_for_migration = NULL;
                        }
                    }
                }
                else if (cell_array(i,23)==3)
                {
                    int pro_loci_left=0;
                    int pro_loci_right=0;
                    int pro_loci_mid=0;
                    for (int x=0; x<new_loci; x++)
                    {
                        if (direction1[x]==2)
                        {
                            pro_loci_left++;
                        }
                        else if (direction1[x]==4)
                        {
                            pro_loci_right++;
                        }
                        else if (direction1[x]==3)
                        {
                            pro_loci_mid++;
                        }
                    }
                    if (density[2]<=mean_density || density[1]<=mean_density || density[3]<=mean_density)
                    {
                        if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                        {
                            int order_loci[3]={2,3,4};
                            int n=3;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            int order_loci[2]={2,4};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                        {
                            int order_loci[2]={3,3};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                        {
                            int order_loci[2]={3,4};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            order=2;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                        {
                            order=3;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            order=4;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            shuffle(direction1,direction1+new_loci,RNG);
                            order=direction1[0];
                        }
                    }
                    else
                    {
                        int new_direction_number=0;
                        for (int dl=0;dl<8;dl++)
                        {
                            if (density[dl]<=mean_density)
                            {
                                new_direction_number=new_direction_number+1;
                            }
                        }
                        if (new_direction_number>0)
                        {
                            int *new_direction_for_migration=new int[new_direction_number];
                            int new_direction_number_for_migration=0;
                            for (int dl=0;dl<8;dl++)
                            {
                                if (density[dl]<=mean_density)
                                {
                                    new_direction_for_migration[new_direction_number_for_migration]=direction[dl];
                                    new_direction_number_for_migration=new_direction_number_for_migration+1;
                                }
                            }
                            int new_loci_for_migratio_number=0;
                            for (int density_nonzero=0; density_nonzero<new_direction_number_for_migration; density_nonzero++)
                            {
                                if (new_direction_for_migration[density_nonzero]>0)
                                {
                                    new_loci_for_migratio_number=new_loci_for_migratio_number+1;
                                }
                            }
                            int *new_loci_for_migration=new int[new_loci_for_migratio_number];
                            int locinumber=0;
                            for (int loci_mig=0;loci_mig<new_direction_number_for_migration;loci_mig++)
                            {
                                if (new_direction_for_migration[loci_mig]>0)
                                {
                                    new_loci_for_migration[locinumber]=new_direction_for_migration[loci_mig];
                                    locinumber=locinumber+1;
                                }
                            }
                            if (locinumber>0)
                            {
                                shuffle(new_loci_for_migration,new_loci_for_migration+locinumber,RNG);
                                order=new_loci_for_migration[0];
                            }
                            delete[] new_direction_for_migration;
                            new_direction_for_migration = NULL;
                            delete[] new_loci_for_migration;
                            new_loci_for_migration = NULL;
                        }
                    }
                }
                else if (cell_array(i,23)==4)
                {
                    int pro_loci_left=0;
                    int pro_loci_right=0;
                    int pro_loci_mid=0;
                    for (int x=0; x<new_loci; x++)
                    {
                        if (direction1[x]==3)
                        {
                            pro_loci_left++;
                        }
                        else if (direction1[x]==5)
                        {
                            pro_loci_right++;
                        }
                        else if (direction1[x]==4)
                        {
                            pro_loci_mid++;
                        }
                    }
                    if (density[3]<=mean_density || density[2]<=mean_density || density[4]<=mean_density)
                    {
                        if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                        {
                            int order_loci[3]={3,4,5};
                            int n=3;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            int order_loci[2]={3,5};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                        {
                            int order_loci[2]={3,4};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                        {
                            int order_loci[2]={4,5};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            order=3;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                        {
                            order=4;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            order=5;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            shuffle(direction1,direction1+new_loci,RNG);
                            order=direction1[0];
                        }
                    }
                    else
                    {
                        int new_direction_number=0;
                        for (int dl=0;dl<8;dl++)
                        {
                            if (density[dl]<=mean_density)
                            {
                                new_direction_number=new_direction_number+1;
                            }
                        }
                        if (new_direction_number>0)
                        {
                            int *new_direction_for_migration=new int[new_direction_number];
                            int new_direction_number_for_migration=0;
                            for (int dl=0;dl<8;dl++)
                            {
                                if (density[dl]<=mean_density)
                                {
                                    new_direction_for_migration[new_direction_number_for_migration]=direction[dl];
                                    new_direction_number_for_migration=new_direction_number_for_migration+1;
                                }
                            }
                            int new_loci_for_migratio_number=0;
                            for (int density_nonzero=0; density_nonzero<new_direction_number_for_migration; density_nonzero++)
                            {
                                if (new_direction_for_migration[density_nonzero]>0)
                                {
                                    new_loci_for_migratio_number=new_loci_for_migratio_number+1;
                                }
                            }
                            int *new_loci_for_migration=new int[new_loci_for_migratio_number];
                            int locinumber=0;
                            for (int loci_mig=0;loci_mig<new_direction_number_for_migration;loci_mig++)
                            {
                                if (new_direction_for_migration[loci_mig]>0)
                                {
                                    new_loci_for_migration[locinumber]=new_direction_for_migration[loci_mig];
                                    locinumber=locinumber+1;
                                }
                            }
                            if (locinumber>0)
                            {
                                shuffle(new_loci_for_migration,new_loci_for_migration+locinumber,RNG);
                                order=new_loci_for_migration[0];
                            }
                            delete[] new_direction_for_migration;
                            new_direction_for_migration = NULL;
                            delete[] new_loci_for_migration;
                            new_loci_for_migration = NULL;
                        }
                    }
                }
                else if (cell_array(i,23)==5)
                {
                    int pro_loci_left=0;
                    int pro_loci_right=0;
                    int pro_loci_mid=0;
                    for (int x=0; x<new_loci; x++)
                    {
                        if (direction1[x]==4)
                        {
                            pro_loci_left++;
                        }
                        else if (direction1[x]==6)
                        {
                            pro_loci_right++;
                        }
                        else if (direction1[x]==5)
                        {
                            pro_loci_mid++;
                        }
                    }
                    if (density[4]<=mean_density || density[3]<=mean_density || density[5]<=mean_density)
                    {
                        if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                        {
                            int order_loci[3]={4,5,6};
                            int n=3;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            int order_loci[2]={4,6};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                        {
                            int order_loci[2]={4,5};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                        {
                            int order_loci[2]={5,6};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            order=4;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                        {
                            order=5;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            order=6;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            shuffle(direction1,direction1+new_loci,RNG);
                            order=direction1[0];
                        }
                    }
                    else
                    {
                        int new_direction_number=0;
                        for (int dl=0;dl<8;dl++)
                        {
                            if (density[dl]<=mean_density)
                            {
                                new_direction_number=new_direction_number+1;
                            }
                        }
                        if (new_direction_number>0)
                        {
                            int *new_direction_for_migration=new int[new_direction_number];
                            int new_direction_number_for_migration=0;
                            for (int dl=0;dl<8;dl++)
                            {
                                if (density[dl]<=mean_density)
                                {
                                    new_direction_for_migration[new_direction_number_for_migration]=direction[dl];
                                    new_direction_number_for_migration=new_direction_number_for_migration+1;
                                }
                            }
                            int new_loci_for_migratio_number=0;
                            for (int density_nonzero=0; density_nonzero<new_direction_number_for_migration; density_nonzero++)
                            {
                                if (new_direction_for_migration[density_nonzero]>0)
                                {
                                    new_loci_for_migratio_number=new_loci_for_migratio_number+1;
                                }
                            }
                            int *new_loci_for_migration=new int[new_loci_for_migratio_number];
                            int locinumber=0;
                            for (int loci_mig=0;loci_mig<new_direction_number_for_migration;loci_mig++)
                            {
                                if (new_direction_for_migration[loci_mig]>0)
                                {
                                    new_loci_for_migration[locinumber]=new_direction_for_migration[loci_mig];
                                    locinumber=locinumber+1;
                                }
                            }
                            if (locinumber>0)
                            {
                                shuffle(new_loci_for_migration,new_loci_for_migration+locinumber,RNG);
                                order=new_loci_for_migration[0];
                            }
                            delete[] new_direction_for_migration;
                            new_direction_for_migration = NULL;
                            delete[] new_loci_for_migration;
                            new_loci_for_migration = NULL;
                        }
                    }
                }
                else if (cell_array(i,23)==6)
                {
                    int pro_loci_left=0;
                    int pro_loci_right=0;
                    int pro_loci_mid=0;
                    for (int x=0; x<new_loci; x++)
                    {
                        if (direction1[x]==5)
                        {
                            pro_loci_left++;
                        }
                        else if (direction1[x]==7)
                        {
                            pro_loci_right++;
                        }
                        else if (direction1[x]==6)
                        {
                            pro_loci_mid++;
                        }
                    }
                    if (density[5]<=mean_density || density[4]<=mean_density || density[6]<=mean_density)
                    {
                        if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                        {
                            int order_loci[3]={5,6,7};
                            int n=3;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            int order_loci[2]={5,7};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                        {
                            int order_loci[2]={5,6};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                        {
                            int order_loci[2]={6,7};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            order=5;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                        {
                            order=6;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            order=7;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            shuffle(direction1,direction1+new_loci,RNG);
                            order=direction1[0];
                        }
                    }
                    else
                    {
                        int new_direction_number=0;
                        for (int dl=0;dl<8;dl++)
                        {
                            if (density[dl]<=mean_density)
                            {
                                new_direction_number=new_direction_number+1;
                            }
                        }
                        if (new_direction_number>0)
                        {
                            int *new_direction_for_migration=new int[new_direction_number];
                            int new_direction_number_for_migration=0;
                            for (int dl=0;dl<8;dl++)
                            {
                                if (density[dl]<=mean_density)
                                {
                                    new_direction_for_migration[new_direction_number_for_migration]=direction[dl];
                                    new_direction_number_for_migration=new_direction_number_for_migration+1;
                                }
                            }
                            int new_loci_for_migratio_number=0;
                            for (int density_nonzero=0; density_nonzero<new_direction_number_for_migration; density_nonzero++)
                            {
                                if (new_direction_for_migration[density_nonzero]>0)
                                {
                                    new_loci_for_migratio_number=new_loci_for_migratio_number+1;
                                }
                            }
                            int *new_loci_for_migration=new int[new_loci_for_migratio_number];
                            int locinumber=0;
                            for (int loci_mig=0;loci_mig<new_direction_number_for_migration;loci_mig++)
                            {
                                if (new_direction_for_migration[loci_mig]>0)
                                {
                                    new_loci_for_migration[locinumber]=new_direction_for_migration[loci_mig];
                                    locinumber=locinumber+1;
                                }
                            }
                            if (locinumber>0)
                            {
                                shuffle(new_loci_for_migration,new_loci_for_migration+locinumber,RNG);
                                order=new_loci_for_migration[0];
                            }
                            delete[] new_direction_for_migration;
                            new_direction_for_migration = NULL;
                            delete[] new_loci_for_migration;
                            new_loci_for_migration = NULL;
                        }
                    }
                }
                else if (cell_array(i,23)==7)
                {
                    int pro_loci_left=0;
                    int pro_loci_right=0;
                    int pro_loci_mid=0;
                    for (int x=0; x<new_loci; x++)
                    {
                        if (direction1[x]==6)
                        {
                            pro_loci_left++;
                        }
                        else if (direction1[x]==8)
                        {
                            pro_loci_right++;
                        }
                        else if (direction1[x]==7)
                        {
                            pro_loci_mid++;
                        }
                    }
                    if (density[6]<=mean_density || density[5]<=mean_density || density[7]<=mean_density)
                    {
                        if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                        {
                            int order_loci[3]={6,7,8};
                            int n=3;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            int order_loci[2]={6,8};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                        {
                            int order_loci[2]={6,7};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                        {
                            int order_loci[2]={7,8};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            order=6;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                        {
                            order=7;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            order=8;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            shuffle(direction1,direction1+new_loci,RNG);
                            order=direction1[0];
                        }
                    }
                    else
                    {
                        int new_direction_number=0;
                        for (int dl=0;dl<8;dl++)
                        {
                            if (density[dl]<=mean_density)
                            {
                                new_direction_number=new_direction_number+1;
                            }
                        }
                        if (new_direction_number>0)
                        {
                            int *new_direction_for_migration=new int[new_direction_number];
                            int new_direction_number_for_migration=0;
                            for (int dl=0;dl<8;dl++)
                            {
                                if (density[dl]<=mean_density)
                                {
                                    new_direction_for_migration[new_direction_number_for_migration]=direction[dl];
                                    new_direction_number_for_migration=new_direction_number_for_migration+1;
                                }
                            }
                            int new_loci_for_migratio_number=0;
                            for (int density_nonzero=0; density_nonzero<new_direction_number_for_migration; density_nonzero++)
                            {
                                if (new_direction_for_migration[density_nonzero]>0)
                                {
                                    new_loci_for_migratio_number=new_loci_for_migratio_number+1;
                                }
                            }
                            int *new_loci_for_migration=new int[new_loci_for_migratio_number];
                            int locinumber=0;
                            for (int loci_mig=0;loci_mig<new_direction_number_for_migration;loci_mig++)
                            {
                                if (new_direction_for_migration[loci_mig]>0)
                                {
                                    new_loci_for_migration[locinumber]=new_direction_for_migration[loci_mig];
                                    locinumber=locinumber+1;
                                }
                            }
                            if (locinumber>0)
                            {
                                shuffle(new_loci_for_migration,new_loci_for_migration+locinumber,RNG);
                                order=new_loci_for_migration[0];
                            }
                            delete[] new_direction_for_migration;
                            new_direction_for_migration = NULL;
                            delete[] new_loci_for_migration;
                            new_loci_for_migration = NULL;
                        }
                    }
                }
                ////////////////////////////////////////////// migration ////////////////////////////////////
                if (order == 1)
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
                }
                else if (order == 2)
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
                }
                else if (order == 3)
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
                    
                }
                else if (order == 4)
                {
                    int cell_label_1=Visual_range(x1,y1,4);
                    int cellstage=Visual_range(x1,y1,3);
                    Visual_range(Range(x1,x1+1),Range(y1,y1+1),all)=0;
                    for (int cor_cell_y=5;cor_cell_y<=8;cor_cell_y++)
                    {
                        cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)+1;
                    }
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=cell_array(i,15);
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                    cell_array(i,23)=4;
                    cell_array(i,20)=0;
                }
                else if (order == 5)
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
                }
                else if (order == 6)
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
                }
                else if (order == 7)
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
                }
                else if (order == 8)
                {
                    int cell_label_1=Visual_range(x1,y1,4);
                    int cellstage=Visual_range(x1,y1,3);
                    Visual_range(Range(x1,x1+1),Range(y1,y1+1),all)=0;
                    for (int cor_cell_y=5;cor_cell_y<=8;cor_cell_y++)
                    {
                        cell_array(i,cor_cell_y)=cell_array(i,cor_cell_y)-1;
                    }
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),1)=1;
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),2)=cell_array(i,15);
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),3)=cellstage;
                    Visual_range(Range(cell_array(i,1),cell_array(i,1)+1),Range(cell_array(i,5),cell_array(i,5)+1),4)=cell_label_1;
                    cell_array(i,23)=8;
                    cell_array(i,20)=0;
                }
                delete[] direction1;
                direction1 = NULL;
                cell_array(i,20)=0;
            }
        }
        else
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
                double density[8]={0};
                area_square.resize(9,9);
                area_square_s=0;
                area_square_s(all,all)=Visual_range(Range(x1-4,x1+4),Range(y1-4,y1+4),4);
                /////////////////////////////////////////////8 directions density dudgement/////////////////////////////////////
                for (int loci_for_mig=0; loci_for_mig<8;loci_for_mig++)
                {
                    if (loci_for_mig==0)
                    {
                        sub_area_square.resize(9,9);
                        sub_area_square_s=0;
                        for (int xss=1; xss<=5; xss++)
                        {
                            for (int yss=1; yss<=5; yss++)
                            {
                                sub_area_square_s(xss,yss)=area_square_s(xss,yss);
                            }
                        }
                        int cell_count[81];
                        int cc=0;
                        for (int cx=0; cx<9; cx++)
                        {
                            for(int cy=0; cy<9; cy++)
                            {
                                cell_count[cc]=sub_area_square_s(cx+1,cy+1);
                                cc++;
                            }
                        }
                        vector<int> mycellcount (cell_count, cell_count+100);
                        sort(mycellcount.begin(),mycellcount.end());
                        mycellcount.erase(unique(mycellcount.begin(), mycellcount.end()), mycellcount.end());
                        long cells_number=0;
                        cells_number = mycellcount.size();
                        if (mycellcount[0]==0)
                        {
                            cells_number=cells_number-1;
                        }
                        density[0]=(double)cells_number/(double)25;
                    }
                    else if (loci_for_mig==2)
                    {
                        sub_area_square.resize(9,9);
                        sub_area_square_s=0;
                        for (int xss=1; xss<=5; xss++)
                        {
                            for (int yss=5; yss<=9; yss++)
                            {
                                sub_area_square_s(xss,yss)=area_square_s(xss,yss);
                            }
                        }
                        int cell_count[81];
                        int cc=0;
                        for (int cx=0; cx<9; cx++)
                        {
                            for(int cy=0; cy<9; cy++)
                            {
                                cell_count[cc]=sub_area_square_s(cx+1,cy+1);
                                cc++;
                            }
                        }
                        vector<int> mycellcount (cell_count, cell_count+100);
                        sort(mycellcount.begin(),mycellcount.end());
                        mycellcount.erase(unique(mycellcount.begin(), mycellcount.end()), mycellcount.end());
                        long cells_number=0;
                        cells_number = mycellcount.size();
                        if (mycellcount[0]==0)
                        {
                            cells_number=cells_number-1;
                        }
                        density[2]=(double)cells_number/(double)25;
                    }
                    else if (loci_for_mig==4)
                    {
                        sub_area_square.resize(9,9);
                        sub_area_square_s=0;
                        for (int xss=5; xss<=9; xss++)
                        {
                            for (int yss=5; yss<=9; yss++)
                            {
                                sub_area_square_s(xss,yss)=area_square_s(xss,yss);
                            }
                        }
                        int cell_count[81];
                        int cc=0;
                        for (int cx=0; cx<9; cx++)
                        {
                            for(int cy=0; cy<9; cy++)
                            {
                                cell_count[cc]=sub_area_square_s(cx+1,cy+1);
                                cc++;
                            }
                        }
                        vector<int> mycellcount (cell_count, cell_count+100);
                        sort(mycellcount.begin(),mycellcount.end());
                        mycellcount.erase(unique(mycellcount.begin(), mycellcount.end()), mycellcount.end());
                        long cells_number=0;
                        cells_number = mycellcount.size();
                        if (mycellcount[0]==0)
                        {
                            cells_number=cells_number-1;
                        }
                        density[4]=(double)cells_number/(double)25;
                    }
                    else if (loci_for_mig==6)
                    {
                        sub_area_square.resize(9,9);
                        sub_area_square_s=0;
                        for (int xss=5; xss<=9; xss++)
                        {
                            for (int yss=1; yss<=5; yss++)
                            {
                                sub_area_square_s(xss,yss)=area_square_s(xss,yss);
                            }
                        }
                        int cell_count[81];
                        int cc=0;
                        for (int cx=0; cx<9; cx++)
                        {
                            for(int cy=0; cy<9; cy++)
                            {
                                cell_count[cc]=sub_area_square_s(cx+1,cy+1);
                                cc++;
                            }
                        }
                        vector<int> mycellcount (cell_count, cell_count+100);
                        sort(mycellcount.begin(),mycellcount.end());
                        mycellcount.erase(unique(mycellcount.begin(), mycellcount.end()), mycellcount.end());
                        long cells_number=0;
                        cells_number = mycellcount.size();
                        if (mycellcount[0]==0)
                        {
                            cells_number=cells_number-1;
                        }
                        density[6]=(double)cells_number/(double)25;
                    }
                    else if (loci_for_mig==1)
                    {
                        sub_area_square.resize(9,9);
                        sub_area_square_s=0;
                        int deltay=0;
                        for (int xss=1; xss<=5; xss++)
                        {
                            for (int yss=1+deltay; yss<=9-deltay; yss++)
                            {
                                sub_area_square_s(xss,yss)=area_square_s(xss,yss);
                            }
                            deltay=deltay+1;
                        }
                        int cell_count[81];
                        int cc=0;
                        for (int cx=0; cx<9; cx++)
                        {
                            for(int cy=0; cy<9; cy++)
                            {
                                cell_count[cc]=sub_area_square_s(cx+1,cy+1);
                                cc++;
                            }
                        }
                        vector<int> mycellcount (cell_count, cell_count+100);
                        sort(mycellcount.begin(),mycellcount.end());
                        mycellcount.erase(unique(mycellcount.begin(), mycellcount.end()), mycellcount.end());
                        long cells_number=0;
                        cells_number = mycellcount.size();
                        if (mycellcount[0]==0)
                        {
                            cells_number=cells_number-1;
                        }
                        density[1]=(double)cells_number/(double)25;
                    }
                    else if (loci_for_mig==3)
                    {
                        sub_area_square.resize(9,9);
                        sub_area_square_s=0;
                        int deltay=0;
                        for (int yss=9; yss>=5; yss--)
                        {
                            for (int xss=1+deltay; xss<=9-deltay; xss++)
                            {
                                sub_area_square_s(xss,yss)=area_square_s(xss,yss);
                            }
                            deltay=deltay+1;
                        }
                        int cell_count[81];
                        int cc=0;
                        for (int cx=0; cx<9; cx++)
                        {
                            for(int cy=0; cy<9; cy++)
                            {
                                cell_count[cc]=sub_area_square_s(cx+1,cy+1);
                                cc++;
                            }
                        }
                        vector<int> mycellcount (cell_count, cell_count+100);
                        sort(mycellcount.begin(),mycellcount.end());
                        mycellcount.erase(unique(mycellcount.begin(), mycellcount.end()), mycellcount.end());
                        long cells_number=0;
                        cells_number = mycellcount.size();
                        if (mycellcount[0]==0)
                        {
                            cells_number=cells_number-1;
                        }
                        density[3]=(double)cells_number/(double)25;
                    }
                    else if (loci_for_mig==5)
                    {
                        sub_area_square.resize(9,9);
                        sub_area_square_s=0;
                        int deltay=0;
                        for (int xss=9; xss>=5; xss--)
                        {
                            for (int yss=1+deltay; yss<=9-deltay; yss++)
                            {
                                sub_area_square_s(xss,yss)=area_square_s(xss,yss);
                            }
                            deltay=deltay+1;
                        }
                        int cell_count[81];
                        int cc=0;
                        for (int cx=0; cx<9; cx++)
                        {
                            for(int cy=0; cy<9; cy++)
                            {
                                cell_count[cc]=sub_area_square_s(cx+1,cy+1);
                                cc++;
                            }
                        }
                        vector<int> mycellcount (cell_count, cell_count+100);
                        sort(mycellcount.begin(),mycellcount.end());
                        mycellcount.erase(unique(mycellcount.begin(), mycellcount.end()), mycellcount.end());
                        long cells_number=0;
                        cells_number = mycellcount.size();
                        if (mycellcount[0]==0)
                        {
                            cells_number=cells_number-1;
                        }
                        density[5]=(double)cells_number/(double)25;
                    }
                    else if (loci_for_mig==7)
                    {
                        sub_area_square.resize(9,9);
                        sub_area_square_s=0;
                        int deltay=0;
                        for (int yss=1; yss<=5; yss++)
                        {
                            for (int xss=1+deltay; xss<=9-deltay; xss++)
                            {
                                sub_area_square_s(xss,yss)=area_square_s(xss,yss);
                            }
                            deltay=deltay+1;
                        }
                        int cell_count[81];
                        int cc=0;
                        for (int cx=0; cx<9; cx++)
                        {
                            for(int cy=0; cy<9; cy++)
                            {
                                cell_count[cc]=sub_area_square_s(cx+1,cy+1);
                                cc++;
                            }
                        }
                        vector<int> mycellcount (cell_count, cell_count+100);
                        sort(mycellcount.begin(),mycellcount.end());
                        mycellcount.erase(unique(mycellcount.begin(), mycellcount.end()), mycellcount.end());
                        long cells_number=0;
                        cells_number = mycellcount.size();
                        if (mycellcount[0]==0)
                        {
                            cells_number=cells_number-1;
                        }
                        density[7]=(double)cells_number/(double)25;
                    }
                }
                int order=0;
                double mean_density=0.6;
                ////////////////////////////////////////////////////initial migration direction dudgement/////////////////////////////////////////////////////////////
                if (cell_array(i,23)==0)
                {
                    int new_direction_number=0;
                    for (int dl=0;dl<8;dl++)
                    {
                        if (density[dl]<=mean_density)
                        {
                            new_direction_number=new_direction_number+1;
                        }
                    }
                    if (new_direction_number>0)
                    {
                        int *new_direction_for_migration=new int[new_direction_number];
                        int new_direction_number_for_migration=0;
                        for (int dl=0;dl<8;dl++)
                        {
                            if (density[dl]<=mean_density)
                            {
                                new_direction_for_migration[new_direction_number_for_migration]=direction[dl];
                                new_direction_number_for_migration=new_direction_number_for_migration+1;
                            }
                        }
                        int new_loci_for_migratio_number=0;
                        for (int density_nonzero=0; density_nonzero<new_direction_number_for_migration; density_nonzero++)
                        {
                            if (new_direction_for_migration[density_nonzero]>0)
                            {
                                new_loci_for_migratio_number=new_loci_for_migratio_number+1;
                            }
                        }
                        int *new_loci_for_migration=new int[new_loci_for_migratio_number];
                        int locinumber=0;
                        for (int loci_mig=0;loci_mig<new_direction_number_for_migration;loci_mig++)
                        {
                            if (new_direction_for_migration[loci_mig]>0)
                            {
                                new_loci_for_migration[locinumber]=new_direction_for_migration[loci_mig];
                                locinumber=locinumber+1;
                            }
                        }
                        if (locinumber>0)
                        {
                            shuffle(new_loci_for_migration,new_loci_for_migration+locinumber,RNG);
                            order=new_loci_for_migration[0];
                        }
                        delete[] new_direction_for_migration;
                        new_direction_for_migration = NULL;
                        delete[] new_loci_for_migration;
                        new_loci_for_migration = NULL;
                        cell_array(i,24)=1;
                    }
                }
                //////////////////////////////////////////////////////following migration direction dudgement////////////////////////////////////////////////////////////
                else if (cell_array(i,23)==1)
                {
                    int pro_loci_left=0;
                    int pro_loci_right=0;
                    int pro_loci_mid=0;
                    for (int x=0; x<new_loci; x++)
                    {
                        if (direction1[x]==2)
                        {
                            pro_loci_left++;
                        }
                        else if (direction1[x]==8)
                        {
                            pro_loci_right++;
                        }
                        else if (direction1[x]==1)
                        {
                            pro_loci_mid++;
                        }
                    }
                    if (density[0]<=mean_density || density[1]<=mean_density || density[7]<=mean_density)
                    {
                        if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                        {
                            int order_loci[3]={8,1,2};
                            int n=3;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            int order_loci[2]={2,8};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                        {
                            int order_loci[2]={1,2};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                        {
                            int order_loci[2]={8,1};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            order=2;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                        {
                            order=1;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            order=8;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            shuffle(direction1,direction1+new_loci,RNG);
                            order=direction1[0];
                        }
                    }
                    else
                    {
                        int new_direction_number=0;
                        for (int dl=0;dl<8;dl++)
                        {
                            if (density[dl]<=mean_density)
                            {
                                new_direction_number=new_direction_number+1;
                            }
                        }
                        if (new_direction_number>0)
                        {
                            int *new_direction_for_migration=new int[new_direction_number];
                            int new_direction_number_for_migration=0;
                            for (int dl=0;dl<8;dl++)
                            {
                                if (density[dl]<=mean_density)
                                {
                                    new_direction_for_migration[new_direction_number_for_migration]=direction[dl];
                                    new_direction_number_for_migration=new_direction_number_for_migration+1;
                                }
                            }
                            int new_loci_for_migratio_number=0;
                            for (int density_nonzero=0; density_nonzero<new_direction_number_for_migration; density_nonzero++)
                            {
                                if (new_direction_for_migration[density_nonzero]>0)
                                {
                                    new_loci_for_migratio_number=new_loci_for_migratio_number+1;
                                }
                            }
                            int *new_loci_for_migration=new int[new_loci_for_migratio_number];
                            int locinumber=0;
                            for (int loci_mig=0;loci_mig<new_direction_number_for_migration;loci_mig++)
                            {
                                if (new_direction_for_migration[loci_mig]>0)
                                {
                                    new_loci_for_migration[locinumber]=new_direction_for_migration[loci_mig];
                                    locinumber=locinumber+1;
                                }
                            }
                            if (locinumber>0)
                            {
                                shuffle(new_loci_for_migration,new_loci_for_migration+locinumber,RNG);
                                order=new_loci_for_migration[0];
                            }
                            delete[] new_direction_for_migration;
                            new_direction_for_migration = NULL;
                            delete[] new_loci_for_migration;
                            new_loci_for_migration = NULL;
                        }
                    }
                }
                else if (cell_array(i,23)==8)
                {
                    int pro_loci_left=0;
                    int pro_loci_right=0;
                    int pro_loci_mid=0;
                    for (int x=0; x<new_loci; x++)
                    {
                        if (direction1[x]==7)
                        {
                            pro_loci_left++;
                        }
                        else if (direction1[x]==1)
                        {
                            pro_loci_right++;
                        }
                        else if (direction1[x]==8)
                        {
                            pro_loci_mid++;
                        }
                    }
                    
                    if (density[7]<=mean_density || density[0]<=mean_density || density[6]<=mean_density)
                    {
                        if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                        {
                            int order_loci[3]={7,8,1};
                            int n=3;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            int order_loci[2]={7,1};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                        {
                            int order_loci[2]={7,8};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                        {
                            int order_loci[2]={8,1};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            order=7;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                        {
                            order=8;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            order=1;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            shuffle(direction1,direction1+new_loci,RNG);
                            order=direction1[0];
                        }
                        
                    }
                    else
                    {
                        int new_direction_number=0;
                        for (int dl=0;dl<8;dl++)
                        {
                            if (density[dl]<=mean_density)
                            {
                                new_direction_number=new_direction_number+1;
                            }
                        }
                        if (new_direction_number>0)
                        {
                            int *new_direction_for_migration=new int[new_direction_number];
                            int new_direction_number_for_migration=0;
                            for (int dl=0;dl<8;dl++)
                            {
                                if (density[dl]<=mean_density)
                                {
                                    new_direction_for_migration[new_direction_number_for_migration]=direction[dl];
                                    new_direction_number_for_migration=new_direction_number_for_migration+1;
                                }
                            }
                            int new_loci_for_migratio_number=0;
                            for (int density_nonzero=0; density_nonzero<new_direction_number_for_migration; density_nonzero++)
                            {
                                if (new_direction_for_migration[density_nonzero]>0)
                                {
                                    new_loci_for_migratio_number=new_loci_for_migratio_number+1;
                                }
                            }
                            int *new_loci_for_migration=new int[new_loci_for_migratio_number];
                            int locinumber=0;
                            for (int loci_mig=0;loci_mig<new_direction_number_for_migration;loci_mig++)
                            {
                                if (new_direction_for_migration[loci_mig]>0)
                                {
                                    new_loci_for_migration[locinumber]=new_direction_for_migration[loci_mig];
                                    locinumber=locinumber+1;
                                }
                            }
                            if (locinumber>0)
                            {
                                shuffle(new_loci_for_migration,new_loci_for_migration+locinumber,RNG);
                                order=new_loci_for_migration[0];
                            }
                            delete[] new_direction_for_migration;
                            new_direction_for_migration = NULL;
                            delete[] new_loci_for_migration;
                            new_loci_for_migration = NULL;
                        }
                    }
                }
                else if (cell_array(i,23)==2)
                {
                    int pro_loci_left=0;
                    int pro_loci_right=0;
                    int pro_loci_mid=0;
                    for (int x=0; x<new_loci; x++)
                    {
                        if (direction1[x]==1)
                        {
                            pro_loci_left++;
                        }
                        else if (direction1[x]==3)
                        {
                            pro_loci_right++;
                        }
                        else if (direction1[x]==2)
                        {
                            pro_loci_mid++;
                        }
                    }
                    
                    if (density[1]<=mean_density || density[2]<=mean_density || density[0]<=mean_density)
                    {
                        if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                        {
                            int order_loci[3]={1,2,3};
                            int n=3;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            int order_loci[2]={1,3};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                        {
                            int order_loci[2]={1,2};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                        {
                            int order_loci[2]={2,3};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            order=1;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                        {
                            order=2;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            order=3;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            shuffle(direction1,direction1+new_loci,RNG);
                            order=direction1[0];
                        }
                    }
                    else
                    {
                        int new_direction_number=0;
                        for (int dl=0;dl<8;dl++)
                        {
                            if (density[dl]<=mean_density)
                            {
                                new_direction_number=new_direction_number+1;
                            }
                        }
                        if (new_direction_number>0)
                        {
                            int *new_direction_for_migration=new int[new_direction_number];
                            int new_direction_number_for_migration=0;
                            for (int dl=0;dl<8;dl++)
                            {
                                if (density[dl]<=mean_density)
                                {
                                    new_direction_for_migration[new_direction_number_for_migration]=direction[dl];
                                    new_direction_number_for_migration=new_direction_number_for_migration+1;
                                }
                            }
                            int new_loci_for_migratio_number=0;
                            for (int density_nonzero=0; density_nonzero<new_direction_number_for_migration; density_nonzero++)
                            {
                                if (new_direction_for_migration[density_nonzero]>0)
                                {
                                    new_loci_for_migratio_number=new_loci_for_migratio_number+1;
                                }
                            }
                            int *new_loci_for_migration=new int[new_loci_for_migratio_number];
                            int locinumber=0;
                            for (int loci_mig=0;loci_mig<new_direction_number_for_migration;loci_mig++)
                            {
                                if (new_direction_for_migration[loci_mig]>0)
                                {
                                    new_loci_for_migration[locinumber]=new_direction_for_migration[loci_mig];
                                    locinumber=locinumber+1;
                                }
                            }
                            if (locinumber>0)
                            {
                                shuffle(new_loci_for_migration,new_loci_for_migration+locinumber,RNG);
                                order=new_loci_for_migration[0];
                            }
                            delete[] new_direction_for_migration;
                            new_direction_for_migration = NULL;
                            delete[] new_loci_for_migration;
                            new_loci_for_migration = NULL;
                        }
                    }
                }
                else if (cell_array(i,23)==3)
                {
                    int pro_loci_left=0;
                    int pro_loci_right=0;
                    int pro_loci_mid=0;
                    for (int x=0; x<new_loci; x++)
                    {
                        if (direction1[x]==2)
                        {
                            pro_loci_left++;
                        }
                        else if (direction1[x]==4)
                        {
                            pro_loci_right++;
                        }
                        else if (direction1[x]==3)
                        {
                            pro_loci_mid++;
                        }
                    }
                    if (density[2]<=mean_density || density[1]<=mean_density || density[3]<=mean_density)
                    {
                        if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                        {
                            int order_loci[3]={2,3,4};
                            int n=3;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            int order_loci[2]={2,4};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                        {
                            int order_loci[2]={2,3};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                        {
                            int order_loci[2]={3,4};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            order=2;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                        {
                            order=3;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            order=4;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            shuffle(direction1,direction1+new_loci,RNG);
                            order=direction1[0];
                        }
                    }
                    else
                    {
                        int new_direction_number=0;
                        for (int dl=0;dl<8;dl++)
                        {
                            if (density[dl]<=mean_density)
                            {
                                new_direction_number=new_direction_number+1;
                            }
                        }
                        if (new_direction_number>0)
                        {
                            int *new_direction_for_migration=new int[new_direction_number];
                            int new_direction_number_for_migration=0;
                            for (int dl=0;dl<8;dl++)
                            {
                                if (density[dl]<=mean_density)
                                {
                                    new_direction_for_migration[new_direction_number_for_migration]=direction[dl];
                                    new_direction_number_for_migration=new_direction_number_for_migration+1;
                                }
                            }
                            int new_loci_for_migratio_number=0;
                            for (int density_nonzero=0; density_nonzero<new_direction_number_for_migration; density_nonzero++)
                            {
                                if (new_direction_for_migration[density_nonzero]>0)
                                {
                                    new_loci_for_migratio_number=new_loci_for_migratio_number+1;
                                }
                            }
                            int *new_loci_for_migration=new int[new_loci_for_migratio_number];
                            int locinumber=0;
                            for (int loci_mig=0;loci_mig<new_direction_number_for_migration;loci_mig++)
                            {
                                if (new_direction_for_migration[loci_mig]>0)
                                {
                                    new_loci_for_migration[locinumber]=new_direction_for_migration[loci_mig];
                                    locinumber=locinumber+1;
                                }
                            }
                            if (locinumber>0)
                            {
                                shuffle(new_loci_for_migration,new_loci_for_migration+locinumber,RNG);
                                order=new_loci_for_migration[0];
                            }
                            delete[] new_direction_for_migration;
                            new_direction_for_migration = NULL;
                            delete[] new_loci_for_migration;
                            new_loci_for_migration = NULL;
                        }
                    }
                }
                else if (cell_array(i,23)==4)
                {
                    int pro_loci_left=0;
                    int pro_loci_right=0;
                    int pro_loci_mid=0;
                    for (int x=0; x<new_loci; x++)
                    {
                        if (direction1[x]==3)
                        {
                            pro_loci_left++;
                        }
                        else if (direction1[x]==5)
                        {
                            pro_loci_right++;
                        }
                        else if (direction1[x]==4)
                        {
                            pro_loci_mid++;
                        }
                    }
                    
                    if (density[3]<=mean_density || density[2]<=mean_density || density[4]<=mean_density)
                    {
                        if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                        {
                            int order_loci[3]={3,4,5};
                            int n=3;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            int order_loci[2]={3,5};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                        {
                            int order_loci[2]={3,4};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                        {
                            int order_loci[2]={4,5};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            order=3;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                        {
                            order=4;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            order=5;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            shuffle(direction1,direction1+new_loci,RNG);
                            order=direction1[0];
                        }
                    }
                    else
                    {
                        int new_direction_number=0;
                        for (int dl=0;dl<8;dl++)
                        {
                            if (density[dl]<=mean_density)
                            {
                                new_direction_number=new_direction_number+1;
                            }
                        }
                        if (new_direction_number>0)
                        {
                            int *new_direction_for_migration=new int[new_direction_number];
                            int new_direction_number_for_migration=0;
                            for (int dl=0;dl<8;dl++)
                            {
                                if (density[dl]<=mean_density)
                                {
                                    new_direction_for_migration[new_direction_number_for_migration]=direction[dl];
                                    new_direction_number_for_migration=new_direction_number_for_migration+1;
                                }
                            }
                            int new_loci_for_migratio_number=0;
                            for (int density_nonzero=0; density_nonzero<new_direction_number_for_migration; density_nonzero++)
                            {
                                if (new_direction_for_migration[density_nonzero]>0)
                                {
                                    new_loci_for_migratio_number=new_loci_for_migratio_number+1;
                                }
                            }
                            int *new_loci_for_migration=new int[new_loci_for_migratio_number];
                            int locinumber=0;
                            for (int loci_mig=0;loci_mig<new_direction_number_for_migration;loci_mig++)
                            {
                                if (new_direction_for_migration[loci_mig]>0)
                                {
                                    new_loci_for_migration[locinumber]=new_direction_for_migration[loci_mig];
                                    locinumber=locinumber+1;
                                }
                            }
                            if (locinumber>0)
                            {
                                shuffle(new_loci_for_migration,new_loci_for_migration+locinumber,RNG);
                                order=new_loci_for_migration[0];
                            }
                            delete[] new_direction_for_migration;
                            new_direction_for_migration = NULL;
                            delete[] new_loci_for_migration;
                            new_loci_for_migration = NULL;
                        }
                    }
                }
                else if (cell_array(i,23)==5)
                {
                    int pro_loci_left=0;
                    int pro_loci_right=0;
                    int pro_loci_mid=0;
                    for (int x=0; x<new_loci; x++)
                    {
                        if (direction1[x]==4)
                        {
                            pro_loci_left++;
                        }
                        else if (direction1[x]==6)
                        {
                            pro_loci_right++;
                        }
                        else if (direction1[x]==5)
                        {
                            pro_loci_mid++;
                        }
                    }
                    if (density[4]<=mean_density  || density[3]<=mean_density || density[5]<=mean_density)
                    {
                        if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                        {
                            int order_loci[3]={4,5,6};
                            int n=3;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            int order_loci[2]={4,6};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                        {
                            int order_loci[2]={4,5};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                        {
                            int order_loci[2]={5,6};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            order=4;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                        {
                            order=5;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            order=6;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            shuffle(direction1,direction1+new_loci,RNG);
                            order=direction1[0];
                        }
                    }
                    else
                    {
                        int new_direction_number=0;
                        for (int dl=0;dl<8;dl++)
                        {
                            if (density[dl]<=mean_density)
                            {
                                new_direction_number=new_direction_number+1;
                            }
                        }
                        if (new_direction_number>0)
                        {
                            int *new_direction_for_migration=new int[new_direction_number];
                            int new_direction_number_for_migration=0;
                            for (int dl=0;dl<8;dl++)
                            {
                                if (density[dl]<=mean_density)
                                {
                                    new_direction_for_migration[new_direction_number_for_migration]=direction[dl];
                                    new_direction_number_for_migration=new_direction_number_for_migration+1;
                                }
                            }
                            int new_loci_for_migratio_number=0;
                            for (int density_nonzero=0; density_nonzero<new_direction_number_for_migration; density_nonzero++)
                            {
                                if (new_direction_for_migration[density_nonzero]>0)
                                {
                                    new_loci_for_migratio_number=new_loci_for_migratio_number+1;
                                }
                            }
                            int *new_loci_for_migration=new int[new_loci_for_migratio_number];
                            int locinumber=0;
                            for (int loci_mig=0;loci_mig<new_direction_number_for_migration;loci_mig++)
                            {
                                if (new_direction_for_migration[loci_mig]>0)
                                {
                                    new_loci_for_migration[locinumber]=new_direction_for_migration[loci_mig];
                                    locinumber=locinumber+1;
                                }
                            }
                            if (locinumber>0)
                            {
                                shuffle(new_loci_for_migration,new_loci_for_migration+locinumber,RNG);
                                order=new_loci_for_migration[0];
                            }
                            delete[] new_direction_for_migration;
                            new_direction_for_migration = NULL;
                            delete[] new_loci_for_migration;
                            new_loci_for_migration = NULL;
                        }
                    }
                }
                else if (cell_array(i,23)==6)
                {
                    int pro_loci_left=0;
                    int pro_loci_right=0;
                    int pro_loci_mid=0;
                    for (int x=0; x<new_loci; x++)
                    {
                        if (direction1[x]==5)
                        {
                            pro_loci_left++;
                        }
                        else if (direction1[x]==7)
                        {
                            pro_loci_right++;
                        }
                        else if (direction1[x]==6)
                        {
                            pro_loci_mid++;
                        }
                    }
                    
                    if (density[5]<=mean_density || density[4]<=mean_density || density[6]<=mean_density)
                    {
                        if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                        {
                            int order_loci[3]={5,6,7};
                            int n=3;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            int order_loci[2]={5,7};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                        {
                            int order_loci[2]={5,6};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                        {
                            int order_loci[2]={6,7};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            order=5;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                        {
                            order=6;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            order=7;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            shuffle(direction1,direction1+new_loci,RNG);
                            order=direction1[0];
                        }
                    }
                    else
                    {
                        int new_direction_number=0;
                        for (int dl=0;dl<8;dl++)
                        {
                            if (density[dl]<=mean_density)
                            {
                                new_direction_number=new_direction_number+1;
                            }
                        }
                        if (new_direction_number>0)
                        {
                            int *new_direction_for_migration=new int[new_direction_number];
                            int new_direction_number_for_migration=0;
                            for (int dl=0;dl<8;dl++)
                            {
                                if (density[dl]<=mean_density)
                                {
                                    new_direction_for_migration[new_direction_number_for_migration]=direction[dl];
                                    new_direction_number_for_migration=new_direction_number_for_migration+1;
                                }
                            }
                            int new_loci_for_migratio_number=0;
                            for (int density_nonzero=0; density_nonzero<new_direction_number_for_migration; density_nonzero++)
                            {
                                if (new_direction_for_migration[density_nonzero]>0)
                                {
                                    new_loci_for_migratio_number=new_loci_for_migratio_number+1;
                                }
                            }
                            int *new_loci_for_migration=new int[new_loci_for_migratio_number];
                            int locinumber=0;
                            for (int loci_mig=0;loci_mig<new_direction_number_for_migration;loci_mig++)
                            {
                                if (new_direction_for_migration[loci_mig]>0)
                                {
                                    new_loci_for_migration[locinumber]=new_direction_for_migration[loci_mig];
                                    locinumber=locinumber+1;
                                }
                            }
                            if (locinumber>0)
                            {
                                shuffle(new_loci_for_migration,new_loci_for_migration+locinumber,RNG);
                                order=new_loci_for_migration[0];
                            }
                            delete[] new_direction_for_migration;
                            new_direction_for_migration = NULL;
                            delete[] new_loci_for_migration;
                            new_loci_for_migration = NULL;
                        }
                    }
                }
                else if (cell_array(i,23)==7)
                {
                    int pro_loci_left=0;
                    int pro_loci_right=0;
                    int pro_loci_mid=0;
                    for (int x=0; x<new_loci; x++)
                    {
                        if (direction1[x]==6)
                        {
                            pro_loci_left++;
                        }
                        else if (direction1[x]==8)
                        {
                            pro_loci_right++;
                        }
                        else if (direction1[x]==7)
                        {
                            pro_loci_mid++;
                        }
                    }
                    if (density[6]<=mean_density || density[5]<=mean_density || density[7]<=mean_density)
                    {
                        if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==1)
                        {
                            int order_loci[3]={6,7,8};
                            int n=3;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            int order_loci[2]={6,8};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_mid==1 && pro_loci_right==0)
                        {
                            int order_loci[2]={6,7};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_right==1 && pro_loci_mid==1 && pro_loci_left==0)
                        {
                            int order_loci[2]={7,8};
                            int n=2;
                            shuffle(order_loci,order_loci+n,RNG);
                            order=order_loci[0];
                        }
                        else if (pro_loci_left==1 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            order=6;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==1)
                        {
                            order=7;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==1 && pro_loci_mid==0)
                        {
                            order=8;
                        }
                        else if (pro_loci_left==0 && pro_loci_right==0 && pro_loci_mid==0)
                        {
                            shuffle(direction1,direction1+new_loci,RNG);
                            order=direction1[0];
                        }
                    }
                    else
                    {
                        int new_direction_number=0;
                        for (int dl=0;dl<8;dl++)
                        {
                            if (density[dl]<=mean_density)
                            {
                                new_direction_number=new_direction_number+1;
                            }
                        }
                        if (new_direction_number>0)
                        {
                            int *new_direction_for_migration=new int[new_direction_number];
                            int new_direction_number_for_migration=0;
                            for (int dl=0;dl<8;dl++)
                            {
                                if (density[dl]<=mean_density)
                                {
                                    new_direction_for_migration[new_direction_number_for_migration]=direction[dl];
                                    new_direction_number_for_migration=new_direction_number_for_migration+1;
                                }
                            }
                            int new_loci_for_migratio_number=0;
                            for (int density_nonzero=0; density_nonzero<new_direction_number_for_migration; density_nonzero++)
                            {
                                if (new_direction_for_migration[density_nonzero]>0)
                                {
                                    new_loci_for_migratio_number=new_loci_for_migratio_number+1;
                                }
                            }
                            int *new_loci_for_migration=new int[new_loci_for_migratio_number];
                            int locinumber=0;
                            for (int loci_mig=0;loci_mig<new_direction_number_for_migration;loci_mig++)
                            {
                                if (new_direction_for_migration[loci_mig]>0)
                                {
                                    new_loci_for_migration[locinumber]=new_direction_for_migration[loci_mig];
                                    locinumber=locinumber+1;
                                }
                            }
                            if (locinumber>0)
                            {
                                shuffle(new_loci_for_migration,new_loci_for_migration+locinumber,RNG);
                                order=new_loci_for_migration[0];
                            }
                            delete[] new_direction_for_migration;
                            new_direction_for_migration = NULL;
                            delete[] new_loci_for_migration;
                            new_loci_for_migration = NULL;
                        }
                    }
                }
                /////////////////////////////////////////////////////migration///////////////////////////////////////
                if (order == 1)
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
                }
                else if (order == 2)
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
                }
                else if (order == 3)
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
                }
                else if (order == 4)
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
                }
                else if (order == 5)
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
                }
                else if (order == 6)
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
                }
                else if (order == 7)
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
                }
                else if (order == 8)
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
                }
                delete[] direction1;
                direction1 = NULL;
                cell_array(i,20)=0;
            }
        }
    }
    else if (cell_array(i,9)==2) //K cells
    {
        if (cell_array(i,14)==0) // big shape
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
            ////////////////////////////////////////////////migration direction dudgement////////////////////////////////////////////
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
                    shuffle(direction1, direction1+new_loci,RNG);
                    order=direction1[0];
                }
                ////////////////////////////////////////////////////////////////////migration//////////////////////////////////////////////////////
                if (order==1)
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
                }
                else if (order==2)
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
                }
                else if (order==3)
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
                }
                else if (order==4)
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
                }
                else if (order==5)
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
                }
                else if (order==6)
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
                }
                else if (order==7)
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
                }
                else if (order==8)
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
                }
                delete[] direction1;
                direction1 = NULL;
                cell_array(i,20)=0;
                cell_array(i,24)=1;
            }
        }
        else
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
                    shuffle(direction1, direction1+new_loci,RNG);
                    order=direction1[0];
                }
                if (order==1)
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
                }
                else if (order==2)
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
                }
                else if (order==3)
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
                }
                else if (order==4)
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
                }
                else if (order==5)
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
                }
                else if (order==6)
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
                }
                else if (order==7)
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
                }
                else if (order==8)
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
                }
                delete[] direction1;
                direction1 = NULL;
                cell_array(i,20)=0;
                cell_array(i,24)=1;
            }
        }
    }
    migration_judgement=migration_judgement+1;
    gsl_rng_free(r6);
}
#endif /* migration_hpp */
