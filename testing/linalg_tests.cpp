#include <iostream>
#include <string>
#include <tuple>
#include <chrono>
#include <assert.h>

#include "math.h"
#include "../linear_algebra.hpp"

/*******************************************************************************************/
// TEST HELPER FUNCTIONS
/*******************************************************************************************/

std::pair<double,uint64_t> test_soln_GJ(const gv::matrix& A, const gv::matrix& x)
{
    gv::matrix b = gv::matrix::matmul(A, x);

    gv::matrix aug = A;
    aug.insert_cols(A.size()[1], b);

    auto start = std::chrono::high_resolution_clock::now();

    gv::matrix x_soln = gv::linalg::solve_GJ(A, b);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    return {(x_soln - x).mag(), duration.count()};
}

std::pair<double,uint64_t> test_invert_GJ(const gv::matrix& A)
{
    auto start = std::chrono::high_resolution_clock::now();

    gv::matrix inv = gv::linalg::invert_GJ(A);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    return {(gv::matrix::matmul(A,inv) - gv::matrix::identity(A.size()[0])).mag(), duration.count()};
}

std::pair<double,uint64_t> test_null_basis_GJ(const gv::matrix& A)
{
    auto tup = A.rref();
    auto pivots = std::get<1>(tup);

    auto start = std::chrono::high_resolution_clock::now();

    gv::matrix nb = gv::linalg::null_basis_GJ(A);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    // null basis should have correct number of vectors
    assert(nb.size()[1] == A.size()[1] - pivots.size());

    if (!nb.is_empty())
        return {(gv::matrix::matmul(A,nb)).mag(), duration.count()};
    else
        return {0, duration.count()};
}

gv::matrix generate_null_basis_mat(const size_t& nrows, const size_t& ncols)
{
    assert(nrows > 0 && ncols > 0);

    size_t k = rand() % nrows;
    std::cout << "k: " << k << std::endl;

    gv::matrix A = gv::matrix::random({nrows - k, ncols}, -100, 100);
    for (size_t i = 0; i < k; ++i)
    {
        gv::matrix sum({1,ncols}, 0);
        for (size_t j = 0; j < nrows - k; ++j)
        {
            //double s = rand() / (double)RAND_MAX);
            double s = 1.0;
            sum += s * A.get_region({j,0}, {1,ncols});
        }
        A.insert_rows((rand() % A.size()[0]), sum);
    }
    return A;
}

/*******************************************************************************************/
// MAIN TEST FUNCTIONS
/*******************************************************************************************/

void LINALG_SOLN_GJ_TEST(size_t num_tests, double min, double max)
{
    std::cout << "\tRUNNING LINALG GJ SOLN TEST..." << std::endl;

    std::vector<double> diffs(num_tests);
    std::vector<uint64_t> dur(num_tests);

    for (size_t i = 0; i < num_tests; ++i)
    {
        gv::matrix A = gv::matrix::random({i+2,i+2}, min, max);
        gv::matrix x = gv::matrix::random({i+2,1}, min, max);
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

        std::cout.precision(5);
        std::cout.width(20);
        std::cout << "\t " << dur[i] << " microsec";

        std::cout << "\n";
    }
}

void LINALG_INV_GJ_TEST(size_t num_tests, double min, double max)
{
    std::cout << "\tRUNNING LINALG GJ INVERT TEST..." << std::endl;

    std::vector<double> diffs(num_tests);
    std::vector<uint64_t> dur(num_tests);

    for (size_t i = 0; i < num_tests; ++i)
    {
        gv::matrix A = gv::matrix::random({i+2,i+2}, min, max);
        std::pair<double, uint64_t> pr = test_invert_GJ(A);
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

        std::cout.precision(5);
        std::cout.width(20);
        std::cout << "\t " << dur[i] << " microsec";

        std::cout << "\n";
    }
}

/*******************************************************************************************/

int main()
{
    
}