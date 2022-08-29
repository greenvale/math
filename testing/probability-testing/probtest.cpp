#include <iostream>

#include <probability.hpp>

int main() 
{
    double myRand = mathlib::Probability::randomNumber();
    
    
    std::cout << myRand << std::endl;
    
}
