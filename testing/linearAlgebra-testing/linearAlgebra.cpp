#include <iostream>
#include <vector> 

#include <matrix.hpp>

int main() 
{
    double data1[] = {1.0, 2.0, 3.0, 4.0};
    double data2[] = {2.0, 3.0, 1.0, 5.0, 3.0, 8.0};
    
    Matrix2<double> mat1(2,2,data1);
    Matrix2<double> mat2(2,3,data2);
    
    Matrix2<double> mat3 = mat1;
    
    Matrix2<double> mat4;
    
    mat4 = mat2;
    
}
