//
//  sortRow.hpp
//  CCDS
//
//  Created by Tao Lee on 11/16/22.
//  Copyright © 2022 Tao Lee. All rights reserved.
//

#ifndef sortRow_hpp
#define sortRow_hpp


#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <omp.h>
#define BZ_THREADSAFE
#define BZ_THREADSAFE_USE_OPENMP

struct DA
{
    double data;
    int index;
};
int cmp1( const void *a ,const void *B)
{
    
    return (*(DA *)a).data > (*(DA *)B).data ? 1 : -1;
}



using namespace blitz;
using namespace std;
void sortRow(Array<float, 2> &cell_array, Array<float, 2> cell_array1, int Col, int colnum_to_sort,int threads)
{
    switch (threads)
    {
        case 1:
        {
            Range all = Range::all();
            int C_16= cell_array.rows();
            
            struct array{
                double data;
                int index;
            };
            
            array *p=new array[C_16];
            
            
            for (int CN=0; CN<C_16; CN++)
            {
                p[CN].data = cell_array(CN+1,colnum_to_sort);
                p[CN].index = CN+1;
            }
            
            qsort(p,C_16,sizeof(p[0]),cmp1);
            
            cell_array1.resize(C_16, Col);
            
            for (int CNx=0; CNx<C_16; CNx++)
            {
                cell_array1(CNx+1,all)=cell_array(p[CNx].index,all);
            }
            
            cell_array(all,all)=0;
            cell_array(all,all)=cell_array1(all,all);
            delete[] p;
            p=NULL;
            break;
        }
            
        default:
        {
            omp_set_num_threads(threads);
            switch (Col)
            {
                case 31:
                {
                    Range all = Range::all();
                    int C_16= cell_array.rows();
                    //    double cell_array_16[C_16];
                    
                    struct array{
                        double data;
                        int index;
                    };
                    
                    array *p=new array[C_16];
                    
                    
                    for (int CN=0; CN<C_16; CN++)
                    {
                        //        cell_array_16[CN]=cell_array(CN+1,colnum_to_sort);
                        p[CN].data = cell_array(CN+1,colnum_to_sort);
                        p[CN].index = CN+1;
                    }
                    
                    //    qsort(p, sizeof(cell_array_16)/sizeof(cell_array_16[0]), sizeof(p[0]), cmp);
                    qsort(p,C_16,sizeof(p[0]),cmp1);
                    
                    cell_array1.resize(C_16, Col);
                    
#pragma omp parallel for schedule(static)
                    {
                        for (int CNx=0; CNx<C_16; CNx++)
                        {
                            cell_array1(CNx+1,1)=cell_array(p[CNx].index,1);
                            cell_array1(CNx+1,2)=cell_array(p[CNx].index,2);
                            cell_array1(CNx+1,3)=cell_array(p[CNx].index,3);
                            cell_array1(CNx+1,4)=cell_array(p[CNx].index,4);
                            cell_array1(CNx+1,5)=cell_array(p[CNx].index,5);
                            cell_array1(CNx+1,6)=cell_array(p[CNx].index,6);
                            cell_array1(CNx+1,7)=cell_array(p[CNx].index,7);
                            cell_array1(CNx+1,8)=cell_array(p[CNx].index,8);
                            cell_array1(CNx+1,9)=cell_array(p[CNx].index,9);
                            cell_array1(CNx+1,10)=cell_array(p[CNx].index,10);
                            cell_array1(CNx+1,11)=cell_array(p[CNx].index,11);
                            cell_array1(CNx+1,12)=cell_array(p[CNx].index,12);
                            cell_array1(CNx+1,13)=cell_array(p[CNx].index,13);
                            cell_array1(CNx+1,14)=cell_array(p[CNx].index,14);
                            cell_array1(CNx+1,15)=cell_array(p[CNx].index,15);
                            cell_array1(CNx+1,16)=cell_array(p[CNx].index,16);
                            cell_array1(CNx+1,17)=cell_array(p[CNx].index,17);
                            cell_array1(CNx+1,18)=cell_array(p[CNx].index,18);
                            cell_array1(CNx+1,19)=cell_array(p[CNx].index,19);
                            cell_array1(CNx+1,20)=cell_array(p[CNx].index,20);
                            cell_array1(CNx+1,21)=cell_array(p[CNx].index,21);
                            cell_array1(CNx+1,22)=cell_array(p[CNx].index,22);
                            cell_array1(CNx+1,23)=cell_array(p[CNx].index,23);
                            cell_array1(CNx+1,24)=cell_array(p[CNx].index,24);
                            cell_array1(CNx+1,25)=cell_array(p[CNx].index,25);
                            cell_array1(CNx+1,26)=cell_array(p[CNx].index,26);
                            cell_array1(CNx+1,27)=cell_array(p[CNx].index,27);
                            cell_array1(CNx+1,28)=cell_array(p[CNx].index,28);
                            cell_array1(CNx+1,29)=cell_array(p[CNx].index,29);
                            cell_array1(CNx+1,30)=cell_array(p[CNx].index,30);
                            cell_array1(CNx+1,31)=cell_array(p[CNx].index,31);
                            //                cell_array1(CNx+1,all)=cell_array(p[CNx].index,all);
                        }
                    }
                    
                    cell_array(all,all)=0;
                    cell_array(all,all)=cell_array1(all,all);
                    delete[] p;
                    p=NULL;
                    break;
                }
                    
                case 28:
                {
                    Range all = Range::all();
                    int C_16= cell_array.rows();
                    struct array{
                        double data;
                        int index;
                    };
                    
                    array *p=new array[C_16];
                    
                    
                    for (int CN=0; CN<C_16; CN++)
                    {
                        p[CN].data = cell_array(CN+1,colnum_to_sort);
                        p[CN].index = CN+1;
                    }
                    
                    qsort(p,C_16,sizeof(p[0]),cmp1);
                    
                    cell_array1.resize(C_16, Col);
                    
#pragma omp parallel for schedule(static)
                    {
                        for (int CNx=0; CNx<C_16; CNx++)
                        {
                            cell_array1(CNx+1,1)=cell_array(p[CNx].index,1);
                            cell_array1(CNx+1,2)=cell_array(p[CNx].index,2);
                            cell_array1(CNx+1,3)=cell_array(p[CNx].index,3);
                            cell_array1(CNx+1,4)=cell_array(p[CNx].index,4);
                            cell_array1(CNx+1,5)=cell_array(p[CNx].index,5);
                            cell_array1(CNx+1,6)=cell_array(p[CNx].index,6);
                            cell_array1(CNx+1,7)=cell_array(p[CNx].index,7);
                            cell_array1(CNx+1,8)=cell_array(p[CNx].index,8);
                            cell_array1(CNx+1,9)=cell_array(p[CNx].index,9);
                            cell_array1(CNx+1,10)=cell_array(p[CNx].index,10);
                            cell_array1(CNx+1,11)=cell_array(p[CNx].index,11);
                            cell_array1(CNx+1,12)=cell_array(p[CNx].index,12);
                            cell_array1(CNx+1,13)=cell_array(p[CNx].index,13);
                            cell_array1(CNx+1,14)=cell_array(p[CNx].index,14);
                            cell_array1(CNx+1,15)=cell_array(p[CNx].index,15);
                            cell_array1(CNx+1,16)=cell_array(p[CNx].index,16);
                            cell_array1(CNx+1,17)=cell_array(p[CNx].index,17);
                            cell_array1(CNx+1,18)=cell_array(p[CNx].index,18);
                            cell_array1(CNx+1,19)=cell_array(p[CNx].index,19);
                            cell_array1(CNx+1,20)=cell_array(p[CNx].index,20);
                            cell_array1(CNx+1,21)=cell_array(p[CNx].index,21);
                            cell_array1(CNx+1,22)=cell_array(p[CNx].index,22);
                            cell_array1(CNx+1,23)=cell_array(p[CNx].index,23);
                            cell_array1(CNx+1,24)=cell_array(p[CNx].index,24);
                            cell_array1(CNx+1,25)=cell_array(p[CNx].index,25);
                            cell_array1(CNx+1,26)=cell_array(p[CNx].index,26);
                            cell_array1(CNx+1,27)=cell_array(p[CNx].index,27);
                            cell_array1(CNx+1,28)=cell_array(p[CNx].index,28);
                        }
                    }
                    
                    cell_array(all,all)=0;
                    cell_array(all,all)=cell_array1(all,all);
                    delete[] p;
                    p=NULL;
                    break;
                }
                    
            }
            break;
        }
    }
}

#endif /* sortRow_hpp */
