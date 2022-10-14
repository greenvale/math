#include <iostream>
#include <vector> 

#include <matrix.hpp>

int main() 
{
    Matrix<double> mat1(3, 3, {{1.0, 2.0, 3.0}, {1.0, 2.0, 2.0}, {2.0, 4.0, 4.0}});
    
    mat1.appendCol(Matrix<double>(1, {1.0, 3.0, 2.0}));
    
    mat1.print();
    
    std::cout << "======" << std::endl;
    
    mat1.echelonForm();
    
    mat1.print();
    
    std::cout << "======" << std::endl;
    
    mat1.reducedEchelonForm();
    
    mat1.print();

    std::cout << mat1.getNumPivotCols() << std::endl;
}
