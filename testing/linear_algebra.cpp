#include <iostream>
#include <string>
#include <tuple>

#include "../linear_algebra.hpp"

int main()
{
    gv::matrix mat1({2,3},
    {{1.0, 2.0, 3.0}, 
    {2.0, 3.0, 4.0}});

    gv::matrix mat2 = {{2,3}, 
    {{2.0,2.0,2.0},
    {3.0,3.0,3.0}}};
 
    gv::matrix mat3 = {{3,6},
    {{1,2,3,4,5,6},{7,8,9,1,2,3},{4,5,6,7,8,9}}};

    gv::matrix mat4 = {{6,3},
    {{1,2,3},{4,5,6},{7,8,9},{1,2,3},{4,5,6},{7,8,9}}};

    gv::matrix mat0 = {{3,3}, 0.0};

    //mat3.insert_rows(0, {{1,6},0.0});
    mat3.insert_rows(0, {{1,6},0.0});
    mat3.insert_rows(3, {{1,6},0.0});
    
    mat3.print();

    std::cout << "\n\n";

    //mat3.sink_zero_rows();

    mat3.print();

}