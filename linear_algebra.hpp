#pragma once
/* Linear Algebra
William Denny
*/

#include <vector>
#include <iostream>
#include <assert.h>
#include <functional>
#include <map>
#include <tuple>
#include <iomanip>
#include <string>
#include <cstring>
#include <algorithm>
#include "math.h"

namespace gv
{

/* ************************************************************************* */

// Matrix class

class matrix
{
private:
    std::pair<size_t,size_t> m_size;
    std::vector<double> m_data;

public:
    matrix();
    matrix(const std::pair<size_t,size_t>& size);
    matrix(const std::pair<size_t,size_t>& size, const double& val);
    matrix(const std::pair<size_t,size_t>& size, const std::vector<std::vector<double>>& data);

    // operator overloading
    friend bool     operator==(const matrix& lh, const matrix& rh);
    friend matrix   operator+ (const matrix& lh, const matrix& rh);
    friend matrix   operator+ (const double& lh, const matrix& rh);
    friend matrix   operator+ (const matrix& lh, const double& rh);
           void     operator+=                  (const matrix& rh);
           void     operator+=                  (const double& rh);
    friend matrix   operator- (const matrix& lh, const matrix& rh);
    friend matrix   operator- (const double& lh, const matrix& rh);
    friend matrix   operator- (const matrix& lh, const double& rh);
           void     operator-=                  (const matrix& rh);
           void     operator-=                  (const double& rh);
    friend matrix   operator* (const matrix& lh, const matrix& rh);
    friend matrix   operator* (const double& lh, const matrix& rh);
    friend matrix   operator* (const matrix& lh, const double& rh);
           void     operator*=                  (const double& rh);
    friend matrix   operator/ (const matrix& lh, const double& rh);
           void     operator/=                  (const double& rh);
           double&  operator[](const std::pair<size_t,size_t>& sub);
    
    // customised uniform operation for all elements given individual function
    void lambda(const std::function<double()>& func);
    void lambda(const std::function<double(double)>& func);
    void lambda(const std::function<double(size_t,size_t)>& func);

    // basic matrix manipulation
    void print() const;
    std::pair<size_t,size_t> size() const;
    size_t ind(const size_t& r, const size_t& c) const;
    double scal() const;
    matrix get_region(const std::pair<size_t,size_t>& sub, const std::pair<size_t,size_t>& size) const;
    void   set_region(const std::pair<size_t,size_t>& sub, const matrix& mat);
    void   insert_rows(const size_t& r, const matrix& mat);
    void   insert_cols(const size_t& c, const matrix& mat);
    void   remove_rows(const size_t& r, const size_t& n);
    void   remove_cols(const size_t& c, const size_t& n);
    matrix resize(const size_t& nrows, const size_t& ncols) const;
    matrix transpose() const;
    double mag() const;
    bool   is_empty() const;
    bool   is_square() const;
    bool   is_symmetric() const;
    std::vector<std::vector<double>> to_vector() const;

    // matrix instance types
    static matrix identity(const size_t& size);
    static matrix diag(const matrix& vals);
    static matrix empty();
    static matrix from_vector(const std::vector<double>& vec);
    static matrix random(const std::pair<size_t,size_t> size, const double& min, const double& max);

    // Gauss-Jordan elimination
    void swap_rows(const size_t& r0, const size_t& r1);
    void swap_cols(const size_t& c0, const size_t& c1);
    void scale_row(const size_t& r, const double& a);
    void axpy_row(const size_t& rx, const size_t& ry, const double& a);
    size_t leading_entry_col(const size_t& r) const;
    bool is_zero_row(const size_t& r) const;
    void sink_zero_rows();
    std::tuple<matrix, std::vector<std::pair<size_t,size_t>>, size_t> reduced_row_echelon_form() const;

    matrix null_basis() const;

    // determinants
    matrix minor(const size_t& r, const size_t& c) const;
    double determinant() const;

};

/* ************************************************************************* */

// linear algebra functions
namespace linalg
{
    matrix solve_GJ(const matrix& A, const matrix& b);
    matrix invert_GJ(const matrix& mat);
};

/* ************************************************************************* */

/* default ctor */
matrix::matrix() 
{
    m_size = {0, 0};
}

/* ctor */
matrix::matrix(const std::pair<size_t,size_t>& size)
{
    m_size = size;
    m_data = std::vector<double>(size.first * size.second);
}

/* ctor for matrix of given size with constant value */ 
matrix::matrix(const std::pair<size_t,size_t>& size, const double& val)
{
    m_size = size;
    m_data = std::vector<double>(size.first * size.second, val);
}

/* ctor with nested vector matrix form */ 
matrix::matrix(const std::pair<size_t,size_t>& size, const std::vector<std::vector<double>>& data)
{
    assert(data.size() == size.first);
    m_size = size;

    for (size_t i = 0; i < m_size.first; ++i)
    {
        assert(data[i].size() == size.second);
        m_data.insert(m_data.end(), data[i].begin(), data[i].end());
    }
}

/* ************************************************************************* 
Operator overloading */

/* comparison operator */
bool operator==(const matrix& lh, const matrix& rh)
{
    if (lh.m_size != rh.m_size)
        return false;
    for (size_t i = 0; i < lh.m_data.size(); ++i)
        if (lh.m_data[i] != rh.m_data[i])
            return false;
    return true;
}

/* matrix + matrix */
matrix operator+(const matrix& lh, const matrix& rh)
{
    assert(lh.m_size == rh.m_size);
    matrix result(lh.m_size);
    for (size_t i = 0; i < result.m_data.size(); ++i)
        result.m_data[i] = lh.m_data[i] + rh.m_data[i];    
    return result;
}

/* double + matrix */
matrix operator+(const double& lh, const matrix& rh)
{
    matrix result(rh.m_size);
    for (size_t i = 0; i < result.m_data.size(); ++i)
        result.m_data[i] = lh + rh.m_data[i];    
    return result;
}

/* matrix + double */
matrix operator+(const matrix& lh, const double& rh)
{
    matrix result(lh.m_size);
    for (size_t i = 0; i < result.m_data.size(); ++i)
        result.m_data[i] = lh.m_data[i] + rh;    
    return result;
}

/* */
void matrix::operator+=(const matrix& rh)
{
    assert(m_size == rh.m_size); 
    for (size_t i = 0; i < m_data.size(); ++i)
        m_data[i] += rh.m_data[i];
}

/* */
void matrix::operator+=(const double& rh)
{ 
    for (size_t i = 0; i < m_data.size(); ++i)
        m_data[i] += rh;
}

/* matrix - matrix */
matrix operator-(const matrix& lh, const matrix& rh)
{
    assert(lh.m_size == rh.m_size);
    matrix result(lh.m_size);
    for (size_t i = 0; i < lh.m_data.size(); ++i)
        result.m_data[i] = lh.m_data[i] - rh.m_data[i];    
    return result;
}

/* double - matrix */
matrix operator-(const double& lh, const matrix& rh)
{
    matrix result(rh.m_size);
    for (size_t i = 0; i < result.m_data.size(); ++i)
        result.m_data[i] = lh - rh.m_data[i];
    return result;
}

/* matrix - double */
matrix operator-(const matrix& lh, const double& rh)
{
    matrix result(lh.m_size);
    for (size_t i = 0; i < result.m_data.size(); ++i)
        result.m_data[i] = lh.m_data[i] - rh;    
    return result;
}

/* */
void matrix::operator-=(const matrix& rh)
{
    assert(m_size == rh.m_size); 
    for (size_t i = 0; i < m_data.size(); ++i)
        m_data[i] -= rh.m_data[i];
}

/* */
void matrix::operator-=(const double& rh)
{
    for (size_t i = 0; i < m_data.size(); ++i)
        m_data[i] -= rh;
}

/* matrix * matrix */
matrix operator*(const matrix& lh, const matrix& rh)
{
    assert(lh.m_size.second == rh.m_size.first); // num cols for lh must equate num rows for rh
    matrix result({lh.m_size.first, rh.m_size.second}, 0.0); // result size is (num rows for lh, num cols for rh)
    
    for (size_t i = 0; i < lh.m_size.first; ++i) // loop through rows in lh
        for (size_t j = 0; j < rh.m_size.second; ++j) // loop through cols in rh
            for (size_t k = 0; k < lh.m_size.second; ++k) // for element (i, j) take dot product of row i in lh and col j in rh
                result.m_data[result.ind(i, j)] += lh.m_data[lh.ind(i, k)] * rh.m_data[rh.ind(k, j)];

    return result;
}

/* scalar * matrix */
matrix operator*(const double& lh, const matrix& rh)
{
    matrix result = rh;
    for (size_t i = 0; i < rh.m_data.size(); ++i)
        result.m_data[i] *= lh;
    return result;
}

/* matrix * scalar */
matrix operator*(const matrix& lh, const double& rh)
{
    matrix result = lh;
    for (size_t i = 0; i < lh.m_data.size(); ++i)
        result.m_data[i] *= rh;
    return result;
}

/* */
void matrix::operator*=(const double& rh)
{
    for (size_t i = 0; i < m_data.size(); ++i)
        m_data[i] *= rh;
}

/* matrix / scalar */
matrix operator/(const matrix& lh, const double& rh)
{
    matrix result = lh;
    for (size_t i = 0; i < lh.m_data.size(); ++i)
        result.m_data[i] /= rh;
    return result;
}

/* */
void matrix::operator/=(const double& rh)
{
    for (size_t i = 0; i < m_data.size(); ++i)
        m_data[i] /= rh;
}

/* */
double& matrix::operator[](const std::pair<size_t,size_t>& sub)
{
    return m_data[ind(sub.first, sub.second)];
}

/* lambda operation for all values in matrix taking no params */
void matrix::lambda(const std::function<double()>& func)
{
    for (size_t i = 0; i < m_data.size(); ++i)
        m_data[i] = func();
}

/* lambda operation for all values in matrix taking value as param */
void matrix::lambda(const std::function<double(double)>& func)
{
    for (size_t i = 0; i < m_data.size(); ++i)
        m_data[i] = func(m_data[i]);
}

/* lambda operation for all values in matrix taking sub as param */
void matrix::lambda(const std::function<double(size_t, size_t)>& func)
{
    for (size_t i = 0; i < m_size.first; ++i)
        for (size_t j = 0; j < m_size.second; ++j)
            m_data[ind(i, j)] = func(i, j);
}

/* ************************************************************************* 
Matrix manipulation */

/* prints matrix to output stream */
void matrix::print() const
{
    for (size_t i = 0; i < m_size.first; ++i)
    {
        for (size_t j = 0; j < m_size.second; ++j)
        {
            if (j > 0)
                std::cout << "\t";
            std::cout.width(10);
            std::cout.precision(3);
            std::cout << m_data[ind(i, j)];
        }
        std::cout << "\n";
    }
}

/* returns size of matrix */
std::pair<size_t,size_t> matrix::size() const
{
    return m_size;
}

/* returns index in valArr of element at position (r, c) - in row-ordered matrix storage */
size_t matrix::ind(const size_t& r, const size_t& c) const
{
    assert(r < m_size.first && c < m_size.second);
    return r * m_size.second + c;
}

/* returns single value in matrix if the matrix is a scalar (i.e. 1x1 matrix) */
double matrix::scal() const
{
    assert(m_size.first == 1 && m_size.second == 1);
    return m_data[0];
}

/* returns matrix of region starting at sub with given size */
matrix matrix::get_region(const std::pair<size_t,size_t>& sub, const std::pair<size_t,size_t>& size) const
{
    assert(sub.first < m_size.first && sub.second < m_size.second);
    assert(sub.first + size.first <= m_size.first && sub.second + size.second <= m_size.second);
    
    matrix result(size);
    for (size_t i = 0; i < size.first; ++i)
        for (size_t j = 0; j < size.second; ++j)
            result.m_data[result.ind(i, j)] = m_data[ind(sub.first + i, sub.second + j)];
    return result;
}

/* sets region starting at sub with given size and given matrix */
void matrix::set_region(const std::pair<size_t,size_t>& sub, const matrix& mat) 
{
    assert(sub.first < m_size.first && sub.second < m_size.second);
    assert(sub.first + mat.m_size.first <= m_size.first && sub.second + mat.m_size.second <= m_size.second);
    
    for (size_t i = 0; i < mat.m_size.first; ++i)
        for (size_t j = 0; j < mat.m_size.second; ++j)
            m_data[ind(sub.first + i, sub.second + j)] = mat.m_data[mat.ind(i, j)];
}

/* insert rows */
void matrix::insert_rows(const size_t& r, const matrix& mat)
{
    assert(r <= m_size.first);
    assert(mat.m_size.second == m_size.second);
    
    m_size.first += mat.m_size.first;
    m_data.insert(m_data.begin() + r*m_size.second, mat.m_data.begin(), mat.m_data.end());    
}

/* insert cols */
void matrix::insert_cols(const size_t& c, const matrix& mat)
{
    assert(c <= m_size.second);
    assert(mat.m_size.first == m_size.first);

    m_size.second += mat.m_size.second;
    for (size_t i = 0; i < m_size.first; ++i)
        m_data.insert(m_data.begin() + i*m_size.second + c, 
            mat.m_data.begin() + i*mat.m_size.second, 
            mat.m_data.begin() + (i+1)*mat.m_size.second);
}

/* remove rows */
void matrix::remove_rows(const size_t& r, const size_t& n)
{
    assert(r + n <= m_size.first);
    
    m_size.first -= n;
    m_data.erase(m_data.begin() + r*m_size.second, m_data.begin() + (r+n)*m_size.second);    
}

/* remove cols */
void matrix::remove_cols(const size_t& c, const size_t& n)
{
    assert(c + n <= m_size.second);

    m_size.second -= n;
    for (size_t i = 0; i < m_size.first; ++i)
        m_data.erase(m_data.begin() + i*m_size.second + c, 
            m_data.begin() + i*m_size.second + c + n);
}

/* returns resized matrix */
matrix matrix::resize(const size_t& nrows, const size_t& ncols) const
{
    assert(nrows * ncols == m_data.size());

    matrix result = *this;
    result.m_size.first = nrows;
    result.m_size.second = ncols;
    return result;
}

/* returns transpose of matrix */
matrix matrix::transpose() const
{
    matrix result({m_size.second, m_size.first});
    
    for (size_t i = 0; i < m_size.first; ++i)
        for (size_t j = 0; j < m_size.second; ++j)
            result.m_data[result.ind(j, i)] = m_data[ind(i, j)];
    
    return result;
}

/* returns the magnitude of matrix by taking L2 norm of all values in vector style */
double matrix::mag() const
{
    double sum = 0.0;
    for (auto& val : m_data)
        sum += val*val;
    return (double)sqrt(sum);
}

/* checks if matrix is empty */
bool matrix::is_empty() const
{
    return (m_size.first == 0 && m_size.second == 0 && m_data.size() == 0);
}

/* checks if matrix is square */
bool matrix::is_square() const
{
    return (m_size.first == m_size.second);
}

/* checks if matrix is symmetric */
bool matrix::is_symmetric() const
{
    if (m_size.first != m_size.second)
        return false;
    
    for (size_t i = 0; i < m_size.first; ++i)
        for (size_t j = 0; j < m_size.first; ++j)
            if (m_data[ind(i, j)] != m_data[ind(j, i)])
                return false;
    
    return true;
}

/* returns vector form of matrix */
std::vector<std::vector<double>> matrix::to_vector() const
{
    std::vector<std::vector<double>> result(m_size.first);
    for (size_t i = 0; i < m_size.first; ++i)
        result[i] = std::vector<double>(m_data.begin() + i*m_size.second, 
            m_data.begin() + (i+1)*m_size.second);
    
    return result;
}

/* ************************************************************************* 
Common matrix types */

/* returns identity matrix with given size */
matrix matrix::identity(const size_t& size)
{
    matrix ident({size, size}, 0.0);
    for (size_t i = 0; i < size; ++i)
        ident.m_data[ident.ind(i, i)] = 1.0;
    
    return ident;
}

/* creates diagonal matrix given col vector of values */
matrix matrix::diag(const matrix& vals)
{
    assert(vals.m_size.first > 0);
    assert(vals.m_size.second == 1);

    matrix result({vals.m_size.first, vals.m_size.first}, 0.0);
    for (size_t i = 0; i < vals.m_size.first; ++i)
        result.m_data[result.ind(i, i)] = vals.m_data[vals.ind(i, 0)];
    
    return result;
}

/* returns empty matrix */
matrix matrix::empty()
{
    matrix result;
    return result;
}

/* returns col matrix from vector */
matrix matrix::from_vector(const std::vector<double>& vec)
{
    matrix result;
    result.m_size = {vec.size(), 1};
    result.m_data = vec;
    return result;
}

/* returns random matrix of given size over given range */
matrix matrix::random(std::pair<size_t,size_t> size, const double& min, const double& max)
{
    matrix result;
    result.m_size = size;
    result.m_data = std::vector<double>(size.first * size.second);
    for (size_t i = 0; i < size.first * size.second; ++i)
        result.m_data[i] = min + ((double)rand() / (double)RAND_MAX)*(max - min);
    return result;
}

/* ************************************************************************* 
Gauss-Jordan elimination */

/* swap two rows */
void matrix::swap_rows(const size_t& r0, const size_t& r1) 
{
    for (size_t i = 0; i < m_size.second; ++i)
        std::swap(m_data[ind(r0, i)], m_data[ind(r1, i)]);
}

/* swap two cols */
void matrix::swap_cols(const size_t& c0, const size_t& c1) 
{
    for (size_t i = 0; i < m_size.second; ++i)
        std::swap(m_data[ind(i, c0)], m_data[ind(i, c1)]);
}

/* scale a row */
void matrix::scale_row(const size_t& r, const double& a)
{
    for (size_t i = 0; i < m_size.second; ++i)
        m_data[ind(r, i)] *= a;
}

/* row axpy */
void matrix::axpy_row(const size_t& ry, const size_t& rx, const double& a)
{
    for (size_t i = 0; i < m_size.second; ++i)
        m_data[ind(ry, i)] += a * m_data[ind(rx, i)];
}

/* get leading entry col for a given row index */
size_t matrix::leading_entry_col(const size_t& r) const
{
    for (size_t i = 0; i < m_size.second; ++i)
        if (m_data[ind(r, i)] != 0.0)
            return i;
    return m_size.second;
}

/* returns true if row is filled with zeros */
bool matrix::is_zero_row(const size_t& r) const
{
    for (size_t i = 0; i < m_size.second; ++i)
        if (m_data[ind(r, i)] != 0.0)
            return false;
    return true;
}

/* send all zero rows to bottom of matrix */
void matrix::sink_zero_rows()
{
    int nzrows = 0;
    int i = 0;
    while (i < m_size.first - nzrows)
    {
        if (is_zero_row(i))
        {
            swap_rows(0, i);
            std::rotate(m_data.begin(), m_data.begin() + m_size.second, m_data.end());
            ++nzrows;
        }
        else
            ++i;
    }
}

/* gauss-jordan elimination to get row-reduced echelon form
    - soln refers to whether a square matrix A, b vector system of equations would have a soln
        given reduced row echelon form on their augmented matrix. The values of soln are enumerated by:
            0 = no solns exist (inconsistent)
            1 = unique soln exists
            2 = infinite solns exist
*/
std::tuple<matrix, std::vector<std::pair<size_t,size_t>>, size_t> matrix::reduced_row_echelon_form() const
{
    // copy matrix to get row-reduced matrix
    matrix rrm(*this);

    // send all zero rows to bottom of the copied matrix
    rrm.sink_zero_rows();

    // locate pivots
    // target_row and target_col are the targets for the next pivot, these are pushed downwards and rightwards as the pivots are found
    // to find a pivot, the rows starting from the target_row are scanned on each column starting at the target column
    // once a pivot is found its row is swapped to the target row and the row is scaled so the pivot element = 1
    // the type of soln present is recorded, this is relevant if A,b sys of equations is given in augmented matrix 
    size_t target_row, target_col;
    size_t soln = 1;
    std::vector<std::pair<size_t,size_t>> pivots;

    target_row = 0;

    // iterate across columns
    for (target_col = 0; target_col < m_size.second; ++target_col)
    {
        // if target_row is out of range then stop searching for pivots
        if (target_row >= m_size.first)
            break;
        
        // iterate down column to find first row with non-zero elem
        for (size_t i = target_row; i < m_size.first; ++i)
        {
            if (rrm.m_data[ind(i, target_col)] != 0.0)
            {
                // swap this row with target_row
                rrm.swap_rows(target_row, i);
                
                // scale row s.t. pivot elem = 1
                rrm.scale_row(target_row, 1.0 / rrm.m_data[ind(target_row, target_col)]);

                // record pivot
                pivots.push_back({target_row, target_col});

                // make sure the pivot is the only non-zero elem in this column
                for (size_t j = 0; j < m_size.first; ++j)
                    if (j != target_row && rrm.m_data[ind(j, target_col)] != 0.0)
                        rrm.axpy_row(j, target_row, -rrm.m_data[ind(j, target_col)] / rrm.m_data[ind(target_row, target_col)]);   

                // for the case of A,b sys of equations augmented matrix
                // if the pivot is in the b cols then this system is inconsistent
                // record this in the soln parameter
                if (target_col >= rrm.m_size.first)
                    soln = 0;

                // increment target_row as this row now has a pivot
                target_row++;

                // stop iterating down column
                break;
            }
        }
    }

    // for the case of A,b sys of equations augmented matrix
    // if there are less pivots than number of rows, the sys has infinite solns
    if (soln != 0 && pivots.size() < rrm.m_size.first)
        soln = 2;

    return {rrm, pivots, soln};
}

/* solve system of equations for unique soln use gauss-jordan elimination
- requires A is square 
- requires b has 1 column
- only provides output if unique solution exists, otherwise returns zero matrix
*/
matrix linalg::solve_GJ(const matrix& A, const matrix& b)
{
    assert(A.is_square() && A.size().first == b.size().first && b.size().second == 1);

    // augment A matrix with b vector
    matrix aug(A);
    aug.insert_cols(A.size().second, b);

    // get reduced row echelon form of augmented matrix
    std::tuple<matrix, std::vector<std::pair<size_t,size_t>>, size_t> tup = aug.reduced_row_echelon_form();

    // if soln exists then augmented matrix, in the case of dim=n, has form [ In  s ]
    // where In is identity matrix of size n and s is the solution vector
    if (std::get<2>(tup) == 1)
        return std::get<0>(tup).get_region({0, aug.size().second - 1}, {aug.size().first, 1});
    else
        return matrix::empty();
}

// inverts a matrix using gauss-jordan elimination
matrix linalg::invert_GJ(const matrix& mat)
{
    assert(mat.is_square());

    // augment matrix with identity matrix
    matrix aug(mat);
    aug.insert_cols(mat.size().second, matrix::identity(mat.size().first));

    // get row reduced echelon form
    std::tuple<matrix, std::vector<std::pair<size_t,size_t>>, size_t> tup = aug.reduced_row_echelon_form();

    // if soln == 1 then matrix invertible so return augmented part as inverted matrix
    if (std::get<2>(tup) == 1)
        return std::get<0>(tup).get_region({0, mat.size().first}, {mat.size().first, mat.size().first});
    else
        return matrix::empty();
}

} // namespace gv