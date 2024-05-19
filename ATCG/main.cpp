//
//  main.cpp
// ATCG
//
//  Created by Tao Lee on 5/21/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//
//
// delta:R1-R0
// ar=70
// density=0.35
// delta=ar-sqrt((((8*density)-1)*ar^2)/3)
//
// ar=70
// delta=12
// density=(ar^2+3*(ar-delta)*(ar-delta))/(8*ar^2)


#include <iostream>
#include <omp.h>
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
#include "density_dependent_growth.hpp"
#include "Low_density_initial_growth.hpp"
#include "free_living_growth.hpp"
#define BZ_THREADSAFE
#include <omp.h>

using namespace std;
using namespace blitz;
#pragma pack(8)

static const char *short_options = "x:y:R:r:m:a:b:d:c:M:S:K:q:t:p:g:D:T:u:z:L:F:k:O";
static const struct option long_options[] = {
    {"Visual_range_x", required_argument, NULL, 'x'},
    {"Visual_range_y", required_argument, NULL, 'y'},
    {"R0", optional_argument, NULL, 'R'},
    {"R1", optional_argument, NULL, 'r'},
    {"mix_ratio_initial", optional_argument, NULL, 'm'},
    {"alpha", optional_argument, NULL, 'a'},
    {"beta", optional_argument, NULL, 'b'},
    {"DDM", optional_argument, NULL, 'd'},
    {"chemotaxis", optional_argument, NULL, 'c'},
    {"bunderD", optional_argument, NULL, 'B'},
    {"migration_rate_r_mean", optional_argument, NULL, 'M'},
    {"migration_rate_r_mean_quia", optional_argument, NULL, 'S'},
    {"migration_rate_K_mean", optional_argument, NULL, 'K'},
    {"beta_distribution_alpha", optional_argument, NULL, 'q'},
    {"beta_distribution_expected", optional_argument, NULL, 't'},
    {"beta_distribution_alpha_mig_time", optional_argument, NULL, 'p'},
    {"beta_distribution_expected_mig_time", optional_argument, NULL, 'g'},
    {"deathjudge", optional_argument, NULL, 'D'},
    {"time_interval", optional_argument, NULL, 'T'},
    {"utralsmall", optional_argument, NULL, 'u'},
    {"allpng", optional_argument, NULL, 'z'},
    {"Low_density_initial", optional_argument, NULL, 'L'},
    {"Single_cell", optional_argument, NULL, 'F'},
    {"K_formation_rate", optional_argument, NULL, 'k'},
    {"threads", optional_argument, NULL, 'O'},
    {"DynamicThreads", optional_argument, NULL, 'P'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char *argv[])
{
    int Visual_range_x=3000;
    int Visual_range_y=3000;
    double R0=60;
    double R1=50;
    double mix_ratio_initial=0.5;
    double alpha=2.2;
    double beta=0;
    int DDM=1;
    int chemotaxis=1;
    double bunderD=0.9;
    double migration_rate_r_mean=200;
    double migration_rate_r_mean_quia=0.5;
    double migration_rate_K_mean=0.25;
    double beta_distribution_alpha=0.01;////
    double beta_distribution_expected=0.15;//////0.2
    double beta_distribution_alpha_mig_time=0.005;
    double beta_distribution_expected_mig_time=0.3;
    double deathjudge=0.005;
    double time_interval=1008;
    int utralsmall=0; //1 yes, 0 no
    int allpng=0;
    int Low_density_initial=1;
    int Single_cell=1;
    int opt = 0;
    double K_formation_rate=0.05;
    int threads=16;
    int DynamicThreads = 1;
    while( (opt = getopt_long(argc, argv, short_options, long_options, NULL)) != -1){
        switch (opt){
            case '?':
                fprintf(stdout, "Usage: %s --Visual_range_x=<int> --Visual_range_y=<int> --R0=<double> --R1=<double> --mix_ratio_initial=<double> --alpha=<double> --beta=<double> --DDM=<int> --chemotaxis=<int> --bunderD=<double> --migration_rate_r_mean=<double> --migration_rate_r_mean_quia=<double> --migration_rate_K_mean=<double> --beta_distribution_alpha=<double> --beta_distribution_expected=<double> --beta_distribution_alpha_mig_time=<double> --beta_distribution_expected_mig_time=<double> --deathjudge=<double> --time_interval=<double> --utralsmall=<int> --allpng=<int> --Low_density_initial=<int>, --Single_cell=<int> --K_formation_rate=<double> --threads=<int> --DynamicThreads=<int>", argv[0]);
                return 0;
            case 'x':
                Visual_range_x = atoi(optarg);
                break;
            case 'y':
                Visual_range_y = atoi(optarg);
                break;
            case 'R':
                if(optarg == NULL){
                    R0 = 60;
                }
                else {
                    R0 = atof(optarg);
                }
                break;
            case 'r':
                if(optarg == NULL){
                    R1 = 50;
                }
                else {
                    R1 = atof(optarg);
                }
                break;
            case 'm':
                if(optarg == NULL){
                    mix_ratio_initial = 0.5;
                }
                else {
                    mix_ratio_initial = atof(optarg);
                }
                break;
            case 'a':
                if(optarg == NULL){
                    alpha = 2.2;
                }
                else {
                    alpha = atof(optarg);
                }
                break;
            case 'b':
                if(optarg == NULL){
                    beta = 0;
                }
                else {
                    beta = atof(optarg);
                }
                break;
            case 'd':
                if(optarg == NULL){
                    DDM = 0;
                }
                else {
                    DDM = atoi(optarg);
                }
                break;
            case 'c':
                if(optarg == NULL){
                    chemotaxis = 0;
                }
                else {
                    chemotaxis = atoi(optarg);
                }
                break;
            case 'B':
                if(optarg == NULL){
                    bunderD = 0.9;
                }
                else {
                    bunderD = atof(optarg);
                }
                break;
            case 'M':
                if(optarg == NULL){
                    migration_rate_r_mean = 5;
                }
                else {
                    migration_rate_r_mean = atof(optarg);
                }
                break;
            case 'S':
                if(optarg == NULL){
                    migration_rate_r_mean_quia = 0.5;
                }
                else {
                    migration_rate_r_mean_quia = atof(optarg);
                }
                break;
            case 'K':
                if(optarg == NULL){
                    migration_rate_K_mean = 0.25;
                }
                else {
                    migration_rate_K_mean = atof(optarg);
                }
                break;

            case 'q':
                if(optarg == NULL){
                    beta_distribution_alpha = 0.03;
                }
                else {
                    beta_distribution_alpha = atof(optarg);
                }
                break;
            case 't':
                if(optarg == NULL){
                    beta_distribution_expected = 0.3;
                }
                else {
                    beta_distribution_expected = atof(optarg);
                }
                break;
            case 'p':
                if(optarg == NULL){
                    beta_distribution_alpha_mig_time = 0.005;
                }
                else {
                    beta_distribution_alpha_mig_time = atof(optarg);
                }
                break;
            case 'g':
                if(optarg == NULL){
                    beta_distribution_expected_mig_time = 0.2;
                }
                else {
                    beta_distribution_expected_mig_time = atof(optarg);
                }
                break;
            case 'D':
                if(optarg == NULL){
                    deathjudge = 0;
                }
                else {
                    deathjudge = atof(optarg);
                }
                break;
            case 'T':
                if(optarg == NULL){
                    time_interval = 501;
                }
                else {
                    time_interval = atof(optarg);
                }
                break;
            case 'u':
                if(optarg == NULL){
                    utralsmall = 0;
                }
                else {
                    utralsmall = atoi(optarg);
                }
                break;
            case 'z':
                if(optarg == NULL){
                    allpng = 1;
                }
                else {
                    allpng = atoi(optarg);
                }
                break;
            case 'L':
                if(optarg == NULL){
                    Low_density_initial = 0;
                }
                else {
                    Low_density_initial = atoi(optarg);
                }
                break;
            case 'F':
                if(optarg == NULL){
                    Single_cell = 0;
                }
                else {
                    Single_cell = atoi(optarg);
                }
                break;
            case 'k':
                if(optarg == NULL){
                    K_formation_rate = 0;
                }
                else {
                    K_formation_rate = atof(optarg);
                }
                break;
            case 'O':
                if(optarg == NULL){
                    threads = 1;
                }
                else {
                    threads = atoi(optarg);
                }
                break;
            case 'P':
                if(optarg == NULL){
                    DynamicThreads = 1;
                }
                else {
                    DynamicThreads = atoi(optarg);
                }
                break;
        }
    }
    if(Single_cell==1)
    {
        mix_ratio_initial=1;
        Low_density_initial=1;
        free_living_growth(Visual_range_x, Visual_range_y, R0, R1, mix_ratio_initial, alpha, beta, DDM, chemotaxis, migration_rate_r_mean, migration_rate_r_mean_quia, migration_rate_K_mean, deathjudge, time_interval, utralsmall, allpng,Single_cell, K_formation_rate,bunderD,beta_distribution_alpha, beta_distribution_expected, beta_distribution_alpha_mig_time, beta_distribution_expected_mig_time,threads, DynamicThreads);
        
    }
    else if(Single_cell==0)
    {
        if (Low_density_initial==0)
        {
            R1=1;
            Low_density_initial_growth(Visual_range_x, Visual_range_y, R0, R1, mix_ratio_initial, alpha, beta, DDM, chemotaxis, migration_rate_r_mean, migration_rate_r_mean_quia, migration_rate_K_mean, deathjudge, time_interval, utralsmall, allpng,bunderD,beta_distribution_alpha, beta_distribution_expected, beta_distribution_alpha_mig_time, beta_distribution_expected_mig_time,threads);
        }
        else
        {
            density_dependent_growth(Visual_range_x, Visual_range_y, R0, R1, mix_ratio_initial, alpha, beta, DDM, chemotaxis, migration_rate_r_mean, migration_rate_r_mean_quia, migration_rate_K_mean, deathjudge, time_interval, utralsmall, allpng,bunderD,beta_distribution_alpha, beta_distribution_expected, beta_distribution_alpha_mig_time, beta_distribution_expected_mig_time,threads, DynamicThreads);
        }
    }
    
    return 0;
}
