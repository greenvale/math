#pragma once

/*  LINEAR ALGEBRA
    ^^^^^^^^^^^^^^

    William Denny (greenvale)

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
#include <array>

namespace gv
{

/* ************************************************************************* */

// Matrix class

class matrix
{
private:
    std::array<size_t,2> m_size;
    std::vector<double> m_data;

public:
    matrix();
    matrix(const std::array<size_t,2>& size);
    matrix(const std::array<size_t,2>& size, const double& val);
    matrix(const std::array<size_t,2>& size, const std::vector<std::vector<double>>& data);

    // operator overloading for element-wise operations
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
    friend matrix   operator/ (const matrix& lh, const matrix& rh);
    friend matrix   operator/ (const matrix& lh, const double& rh);
           void     operator/=                  (const double& rh);
           void     operator/=                  (const matrix& rh);
           double&  operator[](const std::array<size_t,2>& sub);
           

    // matrix multiplication
           matrix   mul                         (const matrix& rh) const;
    static matrix   matmul    (const matrix& lh, const matrix& rh);
    
    // customised uniform operation for all elements given individual function
    void lambda(const std::function<double()>& func);
    void lambda(const std::function<double(double)>& func);
    void lambda(const std::function<double(size_t,size_t)>& func);

    // basic matrix manipulation
    void print() const;
    std::array<size_t,2> size() const;
    size_t ind(const size_t& r, const size_t& c) const;
    double scal() const;
    matrix get_region(const std::array<size_t,2>& sub, const std::array<size_t,2>& size) const;
    void   set_region(const std::array<size_t,2>& sub, const matrix& mat);
    void   insert_rows(const size_t& r, const matrix& mat);
    void   insert_cols(const size_t& c, const matrix& mat);
    void   remove_rows(const size_t& r, const size_t& n);
    void   remove_cols(const size_t& c, const size_t& n);
    void   duplicate_rows(const size_t& r, const size_t& nrows, const size_t& n);
    void   duplicate_cols(const size_t& c, const size_t& ncols, const size_t& n);
    void   move_rows(const size_t& r0, const size_t& r1, const size_t& n);
    void   move_cols(const size_t&c0, const size_t& c1, const size_t& n);
    matrix resize(const size_t& nrows, const size_t& ncols) const;
    matrix transpose() const;
    matrix sum(const int& dim) const;
    double mag() const;
    bool   is_empty() const;
    bool   is_square() const;
    bool   is_symmetric() const;
    std::vector<std::vector<double>> to_vector() const;

    // matrix instance types
    static matrix identity(const size_t& size);
    static matrix diag(const matrix& vals);
    static matrix empty();
    static matrix vector(const std::vector<double>& vec);
    static matrix random(const std::array<size_t,2> size, const double& min, const double& max);

    // elementary row operations
    void swap_rows(const size_t& r0, const size_t& r1);
    void swap_cols(const size_t& c0, const size_t& c1);
    void scale_row(const size_t& r, const double& a);
    void axpy_row(const size_t& rx, const size_t& ry, const double& a);
    size_t leading_entry_col(const size_t& r) const;
    bool is_zero_row(const size_t& r) const;
    void sink_zero_rows();
    std::tuple<matrix, std::vector<size_t>, size_t> rref() const;

    // determinant
    matrix minor(const size_t& r, const size_t& c) const;
    double det() const;

};

/* ************************************************************************* */

// linear algebra functions
namespace linalg
{
    // Gauss-Jordan elimination
    matrix solve_GJ(const matrix& A, const matrix& b);
    matrix invert_GJ(const matrix& mat);
    matrix null_basis_GJ(const matrix& mat);
};

/* ************************************************************************* */

/* default ctor */
matrix::matrix() 
{
    m_size = {0, 0};
}

/* ctor */
matrix::matrix(const std::array<size_t,2>& size)
{
    m_size = size;
    m_data = std::vector<double>(size[0] * size[1]);
}

/* ctor for matrix of given size with constant value */ 
matrix::matrix(const std::array<size_t,2>& size, const double& val)
{
    m_size = size;
    m_data = std::vector<double>(size[0] * size[1], val);
}

/* ctor with nested vector matrix form */ 
matrix::matrix(const std::array<size_t,2>& size, const std::vector<std::vector<double>>& data)
{
    assert(data.size() == size[0]);
    m_size = size;

    for (size_t i = 0; i < m_size[0]; ++i)
    {
        assert(data[i].size() == size[1]);
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

/* matrix * matrix - multiply element by element */
matrix operator*(const matrix& lh, const matrix& rh)
{
    assert(lh.m_size == rh.m_size);
    matrix result = lh;
    for (size_t i = 0; i < lh.m_data.size(); ++i)
        result.m_data[i] *= rh.m_data[i];
    return result;
}

/* matrix multiplication (member function) */
matrix matrix::mul(const matrix& rh) const
{
    assert(m_size[1] == rh.m_size[0]); // num cols for lh must equate num rows for rh
    matrix result({m_size[0], rh.m_size[1]}, 0.0); // result size is (num rows for lh, num cols for rh)
    
    for (size_t i = 0; i < m_size[0]; ++i) // loop through rows in lh
        for (size_t j = 0; j < rh.m_size[1]; ++j) // loop through cols in rh
            for (size_t k = 0; k < m_size[1]; ++k) // for element (i, j) take dot product of row i in lh and col j in rh
                result.m_data[result.ind(i, j)] += m_data[ind(i, k)] * rh.m_data[rh.ind(k, j)];

    return result;
}

/* matrix multiplication (static function) */
matrix matrix::matmul(const matrix& lh, const matrix& rh)
{
    return lh.mul(rh);
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

/* matrix / matrix - divide element by element */
matrix operator/(const matrix &lh, const matrix& rh)
{
    assert(lh.m_size == rh.m_size);
    matrix result = lh;
    for (size_t i = 0; i < lh.m_data.size(); ++i)
        result.m_data[i] /= rh.m_data[i];
    return result;
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

/*  divided one matrix by other element by element */
void matrix::operator/=(const matrix& rh)
{
    assert(m_size == rh.m_size);
    for (size_t i = 0; i < m_data.size(); ++i)
        m_data[i] /= rh.m_data[i];
}

/* */
double& matrix::operator[](const std::array<size_t,2>& sub)
{
    return m_data[ind(sub[0], sub[1])];
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
    for (size_t i = 0; i < m_size[0]; ++i)
        for (size_t j = 0; j < m_size[1]; ++j)
            m_data[ind(i, j)] = func(i, j);
}

/* ************************************************************************* 
Matrix manipulation */

/* prints matrix to output stream */
void matrix::print() const
{
    for (size_t i = 0; i < m_size[0]; ++i)
    {
        for (size_t j = 0; j < m_size[1]; ++j)
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
std::array<size_t,2> matrix::size() const
{
    return m_size;
}

/* returns index in valArr of element at position (r, c) - in row-ordered matrix storage */
size_t matrix::ind(const size_t& r, const size_t& c) const
{
    assert(r < m_size[0] && c < m_size[1]);
    return r * m_size[1] + c;
}

/* returns single value in matrix if the matrix is a scalar (i.e. 1x1 matrix) */
double matrix::scal() const
{
    assert(m_size[0] == 1 && m_size[1] == 1);
    return m_data[0];
}

/* returns matrix of region starting at sub with given size */
matrix matrix::get_region(const std::array<size_t,2>& sub, const std::array<size_t,2>& size) const
{
    assert(sub[0] < m_size[0] && sub[1] < m_size[1]);
    assert(sub[0] + size[0] <= m_size[0] && sub[1] + size[1] <= m_size[1]);
    
    matrix result(size);
    for (size_t i = 0; i < size[0]; ++i)
        for (size_t j = 0; j < size[1]; ++j)
            result.m_data[result.ind(i, j)] = m_data[ind(sub[0] + i, sub[1] + j)];
    return result;
}

/* sets region starting at sub with given size and given matrix */
void matrix::set_region(const std::array<size_t,2>& sub, const matrix& mat) 
{
    assert(sub[0] < m_size[0] && sub[1] < m_size[1]);
    assert(sub[0] + mat.m_size[0] <= m_size[0] && sub[1] + mat.m_size[1] <= m_size[1]);
    
    for (size_t i = 0; i < mat.m_size[0]; ++i)
        for (size_t j = 0; j < mat.m_size[1]; ++j)
            m_data[ind(sub[0] + i, sub[1] + j)] = mat.m_data[mat.ind(i, j)];
}

/* insert rows */
void matrix::insert_rows(const size_t& r, const matrix& mat)
{
    assert(r <= m_size[0]);
    assert(mat.m_size[1] == m_size[1]);
    
    m_size[0] += mat.m_size[0];
    m_data.insert(m_data.begin() + r*m_size[1], mat.m_data.begin(), mat.m_data.end());    
}

/* insert cols */
void matrix::insert_cols(const size_t& c, const matrix& mat)
{
    assert(c <= m_size[1]);
    assert(mat.m_size[0] == m_size[0]);

    m_size[1] += mat.m_size[1];
    for (size_t i = 0; i < m_size[0]; ++i)
        m_data.insert(m_data.begin() + i*m_size[1] + c, 
            mat.m_data.begin() + i*mat.m_size[1], 
            mat.m_data.begin() + (i+1)*mat.m_size[1]);
}

/* remove rows */
void matrix::remove_rows(const size_t& r, const size_t& n)
{
    assert(r + n <= m_size[0]);
    
    m_size[0] -= n;
    m_data.erase(m_data.begin() + r*m_size[1], m_data.begin() + (r+n)*m_size[1]);    
}

/* remove cols */
void matrix::remove_cols(const size_t& c, const size_t& n)
{
    assert(c + n <= m_size[1]);

    m_size[1] -= n;
    for (size_t i = 0; i < m_size[0]; ++i)
        m_data.erase(m_data.begin() + i*m_size[1] + c, 
            m_data.begin() + i*m_size[1] + c + n);
}

/* duplicate rows at position n times */
void matrix::duplicate_rows(const size_t& r, const size_t& nrows, const size_t& n)
{
    assert(r < m_size[0] && r + n <= m_size[0]);

    std::vector<double> tmp(m_data.begin() + m_size[1]*r, m_data.begin() + m_size[1]*(r + nrows));

    m_size[0] += nrows*n;
    for (size_t i = 0; i < n; ++i)
        m_data.insert(m_data.begin() + r*m_size[1], tmp.begin(), tmp.end());
}

/* duplicate cols at position n times */
void matrix::duplicate_cols(const size_t& c, const size_t& ncols, const size_t& n)
{
    assert(c < m_size[1]);

    m_size[1] += ncols*n;
    for (size_t i = 0; i < m_size[0]; ++i)
    {
        std::vector<double> tmp(m_data.begin() + m_size[1]*i + c, m_data.begin() + m_size[1]*i + c+ncols);
        for (size_t j = 0; j < n; ++j)
            m_data.insert(m_data.begin() + i*m_size[1] + c, tmp.begin(), tmp.end());
    }
}

/* moves n rows starting at index r0 to the current r1 position */
void matrix::move_rows(const size_t& r0, const size_t& r1, const size_t& n)
{
    assert(r0 < m_size[0] && r1 <= m_size[0] && r0 + n <= m_size[0]);

    size_t k = 0;
    if (r1 < r0)
        k = n;
    std::vector<double> tmp(m_data.begin() + m_size[1]*r0, 
                            m_data.begin() + m_size[1]*(r0+n));
    m_data.insert(m_data.begin() + m_size[1]*r1,            tmp.begin(), tmp.end());
    m_data.erase( m_data.begin() + m_size[1]*(r0+k), 
                  m_data.begin() + m_size[1]*(r0+k+n));
}

/* moves n cols starting at index c0 to the current c1 position */
void matrix::move_cols(const size_t& c0, const size_t& c1, const size_t& n)
{
    assert(c0 < m_size[1] && c1 <= m_size[1] && c0 + n <= m_size[1]);

    size_t k = 0;
    if (c1 < c0)
        k = n;

    for (size_t i = 0; i < m_size[0]; ++i)
    {
        std::vector<double> tmp(m_data.begin() + m_size[1]*i + c0, 
                                m_data.begin() + m_size[1]*i + c0+n);
        m_data.insert(m_data.begin() + m_size[1]*i + c1,    tmp.begin(), tmp.end());
        m_data.erase( m_data.begin() + m_size[1]*i + c0+k, 
                      m_data.begin() + m_size[1]*i + c0+k+n);
        
    }
}

/* returns resized matrix */
matrix matrix::resize(const size_t& nrows, const size_t& ncols) const
{
    assert(nrows * ncols == m_data.size());

    matrix result = *this;
    result.m_size[0] = nrows;
    result.m_size[1] = ncols;
    return result;
}

/* returns transpose of matrix */
matrix matrix::transpose() const
{
    matrix result({m_size[1], m_size[0]});
    
    for (size_t i = 0; i < m_size[0]; ++i)
        for (size_t j = 0; j < m_size[1]; ++j)
            result.m_data[result.ind(j, i)] = m_data[ind(i, j)];
    
    return result;
}

/* returns sum of all elements in matrix */
matrix matrix::sum(const int& dim) const
{
    assert(dim == -1 || dim == 0 || dim == 1);
    if (dim == -1) // sum all values
    {
        matrix result({1,1});
        result.m_data[0] = 0.0;
        for (auto& d : m_data)
            result.m_data[0] += d;
        return result;        
    }
    else if (dim == 0) // sum along each row
    {
        matrix result({m_size[0], 1});
        for (size_t i = 0; i < m_size[0]; ++i)
        {
            result.m_data[i] = 0.0;
            for (size_t j = 0; j < m_size[1]; ++j)
                result.m_data[i] += m_data[ind(i, j)];
        }
        return result;
    }
    else // dim == 1 : sum along each col
    {
        matrix result({1, m_size[1]});
        for (size_t i = 0; i < m_size[1]; ++i)
        {
            result.m_data[i] = 0.0;
            for (size_t j = 0; j < m_size[0]; ++j)
                result.m_data[i] += m_data[ind(j, i)];
        }
        return result;
    }
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
    return (m_size[0] == 0 && m_size[1] == 0 && m_data.size() == 0);
}

/* checks if matrix is square */
bool matrix::is_square() const
{
    return (m_size[0] == m_size[1]);
}

/* checks if matrix is symmetric */
bool matrix::is_symmetric() const
{
    if (m_size[0] != m_size[1])
        return false;
    
    for (size_t i = 0; i < m_size[0]; ++i)
        for (size_t j = 0; j < m_size[0]; ++j)
            if (m_data[ind(i, j)] != m_data[ind(j, i)])
                return false;
    
    return true;
}

/* returns vector form of matrix */
std::vector<std::vector<double>> matrix::to_vector() const
{
    std::vector<std::vector<double>> result(m_size[0]);
    for (size_t i = 0; i < m_size[0]; ++i)
        result[i] = std::vector<double>(m_data.begin() + i*m_size[1], 
            m_data.begin() + (i+1)*m_size[1]);
    
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
    assert(vals.m_size[0] > 0);
    assert(vals.m_size[1] == 1);

    matrix result({vals.m_size[0], vals.m_size[0]}, 0.0);
    for (size_t i = 0; i < vals.m_size[0]; ++i)
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
matrix matrix::vector(const std::vector<double>& vec)
{
    matrix result;
    result.m_size = {vec.size(), 1};
    result.m_data = vec;
    return result;
}

/* returns random matrix of given size over given range */
matrix matrix::random(std::array<size_t,2> size, const double& min, const double& max)
{
    matrix result;
    result.m_size = size;
    result.m_data = std::vector<double>(size[0] * size[1]);
    for (size_t i = 0; i < size[0] * size[1]; ++i)
        result.m_data[i] = min + ((double)rand() / (double)RAND_MAX)*(max - min);
    return result;
}

/* ************************************************************************* 
Gauss-Jordan elimination */

/* swap two rows */
void matrix::swap_rows(const size_t& r0, const size_t& r1) 
{
    for (size_t i = 0; i < m_size[1]; ++i)
        std::swap(m_data[ind(r0, i)], m_data[ind(r1, i)]);
}

/* swap two cols */
void matrix::swap_cols(const size_t& c0, const size_t& c1) 
{
    for (size_t i = 0; i < m_size[1]; ++i)
        std::swap(m_data[ind(i, c0)], m_data[ind(i, c1)]);
}

/* scale a row */
void matrix::scale_row(const size_t& r, const double& a)
{
    for (size_t i = 0; i < m_size[1]; ++i)
        m_data[ind(r, i)] *= a;
}

/* row axpy */
void matrix::axpy_row(const size_t& ry, const size_t& rx, const double& a)
{
    for (size_t i = 0; i < m_size[1]; ++i)
    {
        m_data[ind(ry, i)] += a * m_data[ind(rx, i)];
    }
}

/* get leading entry col for a given row index */
size_t matrix::leading_entry_col(const size_t& r) const
{
    for (size_t i = 0; i < m_size[1]; ++i)
        if (m_data[ind(r, i)] != 0.0)
            return i;
    return m_size[1];
}

/* returns true if row is filled with zeros */
bool matrix::is_zero_row(const size_t& r) const
{
    for (size_t i = 0; i < m_size[1]; ++i)
        if (m_data[ind(r, i)] != 0.0)
            return false;
    return true;
}

/* send all zero rows to bottom of matrix */
void matrix::sink_zero_rows()
{
    int nzrows = 0;
    int i = 0;
    while (i < m_size[0] - nzrows)
    {
        if (is_zero_row(i))
        {
            swap_rows(0, i);
            std::rotate(m_data.begin(), m_data.begin() + m_size[1], m_data.end());
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
std::tuple<matrix, std::vector<size_t>, size_t> matrix::rref() const
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
    std::vector<size_t> pivots;

    target_row = 0;

    // iterate across columns
    for (target_col = 0; target_col < m_size[1]; ++target_col)
    {
        // if target_row is out of range then stop searching for pivots
        if (target_row >= m_size[0])
            break;
        
        // iterate down column to find first row with non-zero elem
        for (size_t i = target_row; i < m_size[0]; ++i)
        {
            if (rrm.m_data[ind(i, target_col)] != 0.0)
            {
                // swap this row with target_row
                rrm.swap_rows(target_row, i);
                
                // scale row s.t. pivot elem = 1
                rrm.scale_row(target_row, 1.0 / rrm.m_data[ind(target_row, target_col)]);

                // record pivot
                pivots.push_back(target_col);

                // make sure the pivot is the only non-zero elem in this column
                for (size_t j = 0; j < m_size[0]; ++j)
                    if (j != target_row && rrm.m_data[ind(j, target_col)] != 0.0)
                    {
                        rrm.axpy_row(j, target_row, -rrm.m_data[ind(j, target_col)] / rrm.m_data[ind(target_row, target_col)]);
                        rrm.m_data[ind(j, target_col)] = 0.0;
                    }   

                // for the case of A,b sys of equations augmented matrix
                // if the pivot is in the b cols then this system is inconsistent
                // record this in the soln parameter
                if (target_col >= rrm.m_size[0])
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
    if (soln != 0 && pivots.size() < rrm.m_size[0])
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
    assert(A.is_square() && A.size()[0] == b.size()[0] && b.size()[1] == 1);

    // augment A matrix with b vector
    matrix aug(A);
    aug.insert_cols(A.size()[1], b);

    // get reduced row echelon form of augmented matrix
    auto tup = aug.rref();

    // if soln exists then augmented matrix, in the case of dim=n, has form [ In  s ]
    // where In is identity matrix of size n and s is the solution vector
    if (std::get<2>(tup) == 1)
        return std::get<0>(tup).get_region({0, aug.size()[1] - 1}, {aug.size()[0], 1});
    else
        return matrix::empty();
}

// inverts a matrix using gauss-jordan elimination
matrix linalg::invert_GJ(const matrix& mat)
{
    assert(mat.is_square());

    // augment matrix with identity matrix
    matrix aug(mat);
    aug.insert_cols(mat.size()[1], matrix::identity(mat.size()[0]));

    // get row reduced echelon form
    auto tup = aug.rref();

    // if soln == 1 then matrix invertible so return augmented part as inverted matrix
    if (std::get<2>(tup) == 1)
        return std::get<0>(tup).get_region({0, mat.size()[0]}, {mat.size()[0], mat.size()[0]});
    else
        return matrix::empty();
}

// returns nullspace basis by solving equation Ax = 0 for x (x = matrix of nullspace vectors) using Gauss-Jordan elimination
/*  Summary of method:
    - obtains the row reduced form
    - then rearranges columns (and thus variables) to get it in the form
                     r                  k-r
              /                                    \
         r   |    I_[r x r]        C_[r x (k-r)]    |
             |                                      |
        n-r  |  0_[(n-r) x r]    0_[(n-r) x (k-r)]  |
              \                                    /
        
        with dimension [n x k]
        
        where:
            -> r       =  the rank of this matrix
            -> k       =  the number of cols (num of variables)
            -> k-r     =  the nullity (num of null basis vectors)
            -> n       =  the number of rows (num of equations)
            -> I       =  identity
            -> C       =  some matrix
            -> 0       =  zero matrix
            -> [r x c] =  the dimension of each matrix with r = rows, c = cols

    - the matrix of null-space basis vectors then has the form

      /                   \
     |   - C_[r x (k-r)]   |
     |                     |
     |  I_[(k-r) x (k-r)]  |
      \                   /

        with dimension [k x (k-r)]

    this means that multiplying this matrix with the row-reduced rearranged matrix, 
    you get a zero matrix of dimension [n * (k - r)] as expected for k-r null-space basis vectors

    Note on implementation:
    To prevent having to copy original matrix and then rearrange columns to calculate C,
    C is calculated by first taking the 'diminished' part of C which is already on the RHS and then
    adding the non-pivot columns on the RHS in order. Then after adding the identity below, the rows are 
    rearranged using an inverse of the permutation map used to position the non-pivot cols.

*/
matrix linalg::null_basis_GJ(const matrix& mat)
{
    // get row reduced form of matrix
    auto tup = mat.rref();

    matrix mat_rref = std::get<0>(tup);
    auto pivots = std::get<1>(tup);

    // if rank = num_cols then null space = { 0 }
    if (pivots.size() == mat.size()[1])
        return matrix::empty();

    // all cols that don't contain a pivot need to be on LHS
    // iterating in col direction any col without a pivot is removed
    // and placed at the end and the col index is recorded
    std::vector<size_t> nonpiv_cols;
    size_t c = 0;
    size_t i = 0;
    while (i < pivots.size())
    {
        if (pivots[i] != c)
        {
            nonpiv_cols.push_back(c);
            ++c;
        }
        else
        {
            ++i; ++c;
        }
    }

    size_t Cdim_width = mat.size()[1] - *(pivots.end() - 1) - 1;

    // initialise null basis matrix with correct size
    matrix nb = gv::matrix({mat.size()[1], mat.size()[1] - pivots.size()}, 0.0);

    // all non-pivot cols are sent to the RHS of matrix so LHS is diminished C matrix part
    nb.set_region({0,0}, -1.0 * mat_rref.get_region({0, *(pivots.end() - 1) + 1}, {pivots.size(), Cdim_width}));

    for (size_t i = 0; i < nonpiv_cols.size(); ++i)
    {
        nb.set_region({0,Cdim_width+i}, -1.0 * mat_rref.get_region({0, nonpiv_cols[i]}, {mat.size()[0], 1}));
    }

    // add identity part
    nb.set_region({pivots.size(), 0}, gv::matrix::identity(mat.size()[1] - pivots.size()));

    // rearrange rows to match mapping of non-pivot cols
    for (size_t i = 0; i < nonpiv_cols.size(); ++i)
        nb.move_rows(mat.size()[1] - nonpiv_cols.size() + i, nonpiv_cols[i], 1);

    return nb;
}

// returns (r,c) minor of matrix
matrix minor(const size_t& r, const size_t& c) const
{
    matrix result(*this);
    result.remove_rows(r, 1);
    result.remove_cols(c, 1);
    return result;
}

// returns determinant using recursive formula - purely pedagogical as not efficient
double matrix::det() const
{
    // matrix must be square
    assert(m_size[0] == m_size[1]);
    
    if (m_size[0] == 1 && m_size[1] == 1)
        return m_data[0];
    else
    {
        double d = 0;
        double sign = 1;
        for (size_t i = 0; i < m_size[1]; ++i)
        {
            matrix min = minor(0, i);
            d += sign * m_data[ind(0, i)] * min.det();
            sign *= -1.0;
        }
        return det;
    }
}

} // namespace gv