//
//  read_files.hpp
//  ATCG
//
//  Created by Tao Lee on 5/24/24.
//  Copyright © 2024 Tao Lee. All rights reserved.
//

#include "read_files.hpp"

#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "CountLines.hpp"
#include "cell_store.hpp"
#include "cell_trace.hpp"
#include "recovery_parameters.hpp"

using namespace std;
void read_file(CellStore &cells,CellTraceStore &cell_trace, RecoveryParameters &parameters, string Cell_arry_file,string Cell_trace_arry_file,string Parameters)
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
        cells.resize(LINES);
        int cols=cells.column_count();
        for(int row=1;row<=LINES;++row)
        {
            for(int col=1;col<=cols;++col)
            {
                file >> cells.column(col)[row - 1];
            }
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
        cell_trace.resize(LINES1,150);

        int cols=150;
        for(int row=1;row<=LINES1;++row)
        {
            for(int i=0;i<cols;++i)
            {

                    file1 >> cell_trace(row,i+1);
            }
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
            file2.close();
        } else {
            std::cerr << "Can not open the file" << std::endl;
        }
    for (int i=0;i<39;++i)
    {
        parameters(i+1)=array2D[i][0];
    }

}
