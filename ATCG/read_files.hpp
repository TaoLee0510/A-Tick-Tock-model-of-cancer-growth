//
//  read_files.hpp
//  ATCG
//
//  Created by Tao Lee on 5/24/24.
//  Copyright © 2024 Tao Lee. All rights reserved.
//

#ifndef read_files_hpp
#define read_files_hpp

#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include "CountLines.hpp"

using namespace std;
using namespace blitz;
void read_file(Array<double,2> &cell_array,Array<long,2> &cell_trace, Array<double,2> &Parameters_array, string Cell_arry_file,string Cell_trace_arry_file,string Parameters)
{
    
    
    
    ifstream file;
    int LINES;
    file.open(Cell_arry_file,ios::in);
    if(file.fail())
    {
        cout<<"File not exits."<<endl;
        file.close();
    }
    else//file exits
    {
        LINES=CountLines(Cell_arry_file);
        cell_array.resize(LINES,31);
        int cols=31;
        int num=0;
        while(!file.eof()) //read file
        {
            for(int i=0;i<cols;++i)
            {
                    file >> cell_array(num+1,i+1);
            }
            num++;
        }
    }
    file.close(); //Close File
    
    
    
    ifstream file1;
    int LINES1;
    file1.open(Cell_trace_arry_file,ios::in);
    if(file1.fail())
    {
        cout<<"File not exits."<<endl;
        file1.close();
    }
    else//file exits
    {
        LINES1=CountLines(Cell_trace_arry_file);
        cell_trace.resize(LINES,150);
        int cols=150;
        int num=0;
        while(!file1.eof()) //read file
        {
            for(int i=1;i<=cols;++i)
            {

                    file1 >> cell_trace(num+1,i+1);
            }
            num++;
        }
    }
    file1.close(); //Close File
    

    
    std::ifstream file2(Parameters);
    std::vector<std::vector<int>> array2D;
    std::string line;
     
        if (file2.is_open()) {
            while (std::getline(file2, line)) {
                std::istringstream iss(line);
                std::vector<int> row;
                std::string dummy;
                // abandon first two cols
                if (!(iss >> dummy >> dummy)) {
                    continue;
                }
                // read from the third line
                int data;
                while (iss >> data) {
                    row.push_back(data);
                }
                if (!row.empty()) {
                    array2D.push_back(row);
                }
            }
            file.close();
        } else {
            std::cerr << "Can not open the file" << std::endl;
        }
    for (int i=0;i<39;++i)
    {
        Parameters_array(i+1,1)=array2D[i][0];
    }
    
}

#endif /* read_files_hpp */
