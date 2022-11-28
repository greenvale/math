#include <iostream>
#include <vector> 

#include "../../LinearAlgebra.hpp"

int main() 
{
    /*
    mathlib::Matrix A({2, 2}, {{1.0, -1.0}, {2.0, 1.0}});
    mathlib::Matrix b({2, 1}, {{3.0}, {2.0}});

    mathlib::Matrix sys = A;
    sys.insertCols(2, b);

    sys.display();

    std::cout << std::endl;

    sys.axpyRow(0, 1, -0.5);

    sys.display();
    */

    mathlib::Matrix mat1({5, 3}, {{0.0, 0.0, 0.0}, {1.0, 2.0, 0.0},  {1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {2.0, 0.0, 1.0}});

    mat1.display();

    std::cout << std::endl;

    mat1.echelonForm();

    mat1.display();

}
