//
//  CountLines.hpp
//  ATCG
//
//  Created by Tao Lee on 5/24/24.
//  Copyright © 2024 Tao Lee. All rights reserved.
//

#ifndef CountLines_hpp
#define CountLines_hpp

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int CountLines(string Cell_arry_file)
{
    ifstream ReadFile;
    int n=0;
    string tmp;
    ReadFile.open(Cell_arry_file,ios::in);
    if(ReadFile.fail())
    {
        return 0;
    }
    else
    {
        while(getline(ReadFile,tmp,'\n'))
        {
            n++;
        }
        ReadFile.close();
        return n;
    }
}

#endif /* CountLines_hpp */
