#include <iostream>
#include <vector> 

#include "../../LinearAlgebra.hpp"

int main() 
{

    //mathlib::Matrix A({4, 5}, {{1.0, 0.0, 1.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 0.0, 1.0, 0.0}, {0.0, 1.0, 0.0, 0.0, 1.0}});
    mathlib::Matrix A({4, 5}, {{1.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 1.0}});

    std::tuple<mathlib::Matrix, std::vector<std::vector<unsigned int>>, unsigned int> tup = A.rowRedEchelonForm();

    mathlib::Matrix rr = std::get<0>(tup);
    
    rr.display();

    std::cout << std::endl;

    mathlib::Matrix ns = A.nullBasis();
    ns.display();
}