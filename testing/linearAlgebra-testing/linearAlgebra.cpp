#include <iostream>
#include <vector> 

#include <matrix.hpp>

int main() 
{
    double data1[] = {1.0, 2.0, 3.0, 4.0};
    double data2[] = {2.0, 3.0, 1.0, 5.0, 3.0, 8.0, 1.0, 2.0, 3.0};
    
    Matrix<double> mat1(2,2,data1);
    Matrix<double> mat2(3,3,data2);
    
    mat2.print();
    
    std::cout<<"========="<<std::endl;
    
    mat2.appendCol(Matrix<double>(3, 1, {{1.1}, {1.2}, {1.3}}));
    
    mat2.removeCol(0);
    
    mat2.print();
}
