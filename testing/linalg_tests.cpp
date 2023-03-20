#include <iostream>
#include <string>
#include <tuple>
#include <chrono>

#include "math.h"
#include "../linear_algebra.hpp"

/*******************************************************************************************/
// TEST HELPER FUNCTIONS
/*******************************************************************************************/

std::pair<double,uint64_t> test_soln_GJ(const gv::matrix& A, const gv::matrix& x)
{
    gv::matrix b = A * x;

    gv::matrix aug = A;
    aug.insert_cols(A.size().second, b);

    //std::tuple<gv::matrix, std::vector<std::pair<size_t,size_t>>, size_t> tup = aug.reduced_row_echelon_form();
    //std::get<0>(tup).print();

    auto start = std::chrono::high_resolution_clock::now();

    gv::matrix x_soln = gv::linalg::solve_GJ(A, b);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    /*
    // PRINT SOLUTION INFO
    std::cout << "SYSTEM: \nA:\n";
    A.print();
    std::cout << "b:\n";
    b.print();
    std::cout << "SOLUTION: \n";
    if (x_soln.is_empty())
        std::cout << "\n==!! NO UNIQUE SOLN FOUND !!==\n";
    x_soln.print();
    std::cout << "SHOULD BE: \n";
    x.print();
    */

    return {(x_soln - x).mag(), duration.count()};
}

gv::matrix rand_matrix(std::pair<size_t,size_t> size, const double& min, const double& max)
{
    std::vector<double> vals(size.first * size.second);
    for (int i = 0; i < size.first * size.second; ++i)
        vals[i] = min + ((double)rand() / (double)RAND_MAX) * (max - min);
    return gv::matrix::from_vector(vals).resize(size.first, size.second);
}

/*******************************************************************************************/
// MAIN TEST FUNCTIONS
/*******************************************************************************************/

void LINALG_SOLN_GJ_TEST(int num_tests, double min, double max, const double& tol)
{
    std::cout << "\tRUNNING LINALG GJ SOLN TEST..." << std::endl;

    int passes = 0;
    std::vector<double> diffs(num_tests);
    std::vector<uint64_t> dur(num_tests);

    for (int i = 0; i < num_tests; ++i)
    {
        gv::matrix A = rand_matrix({i+2,i+2}, min, max);
        gv::matrix x = rand_matrix({i+2,1}, min, max);
        std::pair<double, uint64_t> pr = test_soln_GJ(A, x);
        diffs[i] = pr.first;
        dur[i] = pr.second;

        std::cout.precision(5);
        std::cout.width(20);
        std::cout << "TEST " << i << "\t";
        
        std::cout.precision(5);
        std::cout.width(20);
        std::cout << i+2 << "x" << i+2;

        std::cout.precision(5);
        std::cout.width(20);
        std::cout << "\t " << diffs[i];
        passes += (diffs[i] > tol);

        std::cout.precision(5);
        std::cout.width(20);
        std::cout << "\t " << dur[i] << " microsec";

        std::cout << "\n";
    }

    std::cout << "\n\n" << passes << "/" << num_tests << " LINALG GJ SOLN TESTS PASSED " << tol << " TOL\n\n";
}

/*******************************************************************************************/

int main()
{
    LINALG_SOLN_GJ_TEST(1000, -10000, 10000, 1e-8);
}