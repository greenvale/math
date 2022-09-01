#include <iostream>
#include <vector> 

#include <linearAlgebra.hpp>

int main() 
{
    std::vector<std::vector<double>> inputs = { {2.0, 2.0, 4.0, 18.0}, {3.0, -1.0, 2.0, 7.0}, {1.5, 0.5, -1.0, -0.5} };
    
    mathlib::LinearAlgebra<double>::printMatrix_stdvec(inputs);
    
    std::vector<double> result = mathlib::LinearAlgebra<double>::gaussElim_stdvec(3, inputs);
    
    mathlib::LinearAlgebra<double>::printMatrix_stdvec({ result });
}
