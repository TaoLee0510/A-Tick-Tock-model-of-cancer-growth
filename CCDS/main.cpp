//
//  main.cpp
//  CCDS
//
//  Created by Tao Lee on 5/21/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

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
#include "density_dependent_growth.hpp"
#include "Low_density_initial_growth.hpp"

using namespace std;
using namespace blitz;
#pragma pack(8)

static const char *short_options = "x:y:R:r:m:a:b:d:c:M:S:K:D:T:u:p:L:F";
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
    {"migration_rate_r_mean", optional_argument, NULL, 'M'},
    {"migration_rate_r_mean_quia", optional_argument, NULL, 'S'},
    {"migration_rate_K_mean", optional_argument, NULL, 'K'},
    {"deathjudge", optional_argument, NULL, 'D'},
    {"time_interval", optional_argument, NULL, 'T'},
    {"utralsmall", optional_argument, NULL, 'u'},
    {"allpng", optional_argument, NULL, 'p'},
    {"Low_density_initial", optional_argument, NULL, 'L'},
    {"free_living", optional_argument, NULL, 'F'},
    {NULL, 0, NULL, 0}
};

int main(int argc, char *argv[])
{
    int Visual_range_x=1000;
    int Visual_range_y=1000;
    double R0=60;
    double R1=1;
    double mix_ratio_initial=0.5;
    float alpha=2.2;
    float beta=0;
    int DDM=0;
    int chemotaxis=0;
    double migration_rate_r_mean=5;
    double migration_rate_r_mean_quia=0.5;
    double migration_rate_K_mean=0.25;
    double deathjudge=0;
    double time_interval=505;
    int utralsmall=0; //1 yes, 0 no
    int allpng=0;
    int Low_density_initial=0;
    int free_living=0;
    int opt = 0;
    while( (opt = getopt_long(argc, argv, short_options, long_options, NULL)) != -1){
        switch (opt){
            case '?':
                fprintf(stdout, "Usage: %s --Visual_range_x=<int> --Visual_range_y=<int> --R0=<double> --R1=<double> --mix_ratio_initial=<double> --alpha=<float> --beta=<float> --DDM=<int> --chemotaxis=<int> --migration_rate_r_mean=<double> --migration_rate_r_mean_quia=<double> --migration_rate_K_mean=<double> --deathjudge=<double> --time_interval=<double> --utralsmall=<int> --allpng=<int> --Low_density_initial=<int>, --free_living=<int>", argv[0]);
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
            case 'p':
                if(optarg == NULL){
                    allpng = 0;
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
                    free_living = 0;
                }
                else {
                    free_living = atoi(optarg);
                }
                break;
        }
    }
    if(free_living==1)
    {
        mix_ratio_initial=1;
    }
    if (Low_density_initial==0)
    {
        Low_density_initial_growth(Visual_range_x, Visual_range_y, R0, R1, mix_ratio_initial, alpha, beta, DDM, chemotaxis, migration_rate_r_mean, migration_rate_r_mean_quia, migration_rate_K_mean, deathjudge, time_interval, utralsmall, allpng,free_living);
    }
    else
    {
        density_dependent_growth(Visual_range_x, Visual_range_y, R0, R1, mix_ratio_initial, alpha, beta, DDM, chemotaxis, migration_rate_r_mean, migration_rate_r_mean_quia, migration_rate_K_mean, deathjudge, time_interval, utralsmall, allpng,free_living);
    }
    return 0;
}
