/*
Linear Algebra library
William Denny, 1st September 2022
*/

#include <iostream>
#include <vector>
#include <assert.h>
#include <array>

namespace mathlib
{
    
    template <class T>
    class LinearAlgebra
    {
    public:
        static std::vector<T> gaussElim_stdvec(const int& n, const std::vector<std::vector<T>>& matrix);
        
        static void printMatrix_stdvec(const std::vector<std::vector<T>>& matrix);
    };
    
    //====================================================================================================================
    /*
    Gauss's method for solving linear system of equations (not using linear algebra class - just vector class)
    */
    template <class T>
    std::vector<T> LinearAlgebra<T>::gaussElim_stdvec(const int& n, const std::vector<std::vector<T>>& matrix)
    {
        std::vector<std::vector<T>> matrix_cpy = matrix;
        
        // first check dimensions of matrix rows are all correct
        assert(matrix_cpy.size() == n);
        for (int i = 0; i < n; ++i)
        {
            assert(matrix_cpy[i].size() == n+1);
        }
        
        // eliminate lower triangular elements to get upper triangular matrix
        for (int i = 0; i < n; ++i)
        {
            // normalise row against ith entry
            T normEntry = matrix_cpy[i][i];
            assert(normEntry != 0.0);
            
            for (int j = 0; j < n + 1; ++j)
            {
                matrix_cpy[i][j] /= normEntry;
            }
            
            // eliminate this variable from row below by subtracting row multiplied by coefficient of same variable below
            if (i < n - 1)
            {
                for (int j = i + 1; j < n; ++j)
                {
                    T elimEntry = matrix_cpy[j][i];
                    
                    for (int k = 0; k < n + 1; ++k)
                    {
                        matrix_cpy[j][k] -= elimEntry * matrix_cpy[i][k];
                    }
                }
            }
        }
        
        // eliminate upper triangular elements to get diagonal matrix
        for (int i = n - 1; i > 0; --i)
        {
            for (int j = 0; j < i; ++j)
            {
                matrix_cpy[j][n] -= matrix_cpy[j][i] * matrix_cpy[i][n];
                matrix_cpy[j][i] = 0.0;
            }
        }
        
        // extract results
        std::vector<T> result = {};
        
        for (int i = 0; i < n; ++i)
        {
            result.push_back(matrix_cpy[i][n]);
        }
        
        return result;
    }
    
    //====================================================================================================================
    /* 
    Print matrix implemented using std::vector
    */
    template <class T>
    void LinearAlgebra<T>::printMatrix_stdvec(const std::vector<std::vector<T>>& matrix)
    {
        std::cout << "===================================" << std::endl;
        
        int colNum = matrix[0].size();
        for (int i = 0; i < matrix.size(); ++i)
        {
            assert(matrix[i].size() == colNum);
            
            for (int j = 0; j < matrix[i].size(); ++j)
            {
                std::cout << matrix[i][j] << "\t";
            }
            
            std::cout << std::endl;
        }
        
        std::cout << "===================================" << std::endl;
    }
}
