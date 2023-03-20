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
 
    gv::matrix mat3 = {{5,6},
    {{1,2,3,4,5,6},{7,8,9,1,2,3},{4,5,6,7,8,9}, {9,8,7,6,5,4}, {6,5,4,3,2,1}}};

    gv::matrix mat4 = {{6,3},
    {{1,2,3},{4,5,6},{7,8,9},{1,2,3},{4,5,6},{7,8,9}}};

    gv::matrix mat0 = {{3,3}, 0.0};

    gv::matrix zero({1,6},0.0);

    mat3.print();
    std::cout << "\n\n";

    std::tuple<gv::matrix, std::vector<std::pair<size_t,size_t>>, size_t> res = mat3.reduced_row_echelon_form();
    
    std::get<0>(res).print();
    std::cout << std::get<2>(res) << std::endl;
}