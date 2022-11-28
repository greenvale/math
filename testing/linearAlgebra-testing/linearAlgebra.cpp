#include <iostream>
#include <vector> 

#include "../../LinearAlgebra.hpp"

int main() 
{

    mathlib::Matrix A({2, 2}, {{2.0, 1.0}, {1.0, -1.0}});
    mathlib::Matrix b({2, 1}, {{5.0}, {1.0}});

    mathlib::Matrix x = mathlib::Matrix::solve_GJ(A, b);

    x.display();
}
