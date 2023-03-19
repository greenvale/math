/* Linear Algebra
William Denny
*/
#pragma once

#include <vector>
#include <iostream>
#include <assert.h>
#include <functional>
#include <map>
#include <tuple>
#include <iomanip>
#include <string>
#include <cstring>

namespace gv
{

/* **************************************************************************************************
    MATRIX
************************************************************************************************** */

class matrix
{
private:
    std::pair<size_t,size_t> m_size;
    size_t m_num_elems;
    double* m_data;

public:
    matrix();
    matrix(const std::pair<size_t,size_t>& size);
    matrix(const std::pair<size_t,size_t>& size, const double& val);
    matrix(const std::pair<size_t,size_t>& size, const std::vector<std::vector<double>>& data);
    matrix(const matrix& rh); // copy ctor
    matrix(matrix&& rh); // move ctor
    ~matrix(); // dtor

    // operator overloading
           matrix&  operator=(const matrix& rh); // copy assignment operator
           matrix&  operator=(matrix&& rh); // move assignment operator
    friend bool     operator==(const matrix& lh, const matrix& rh);
    friend matrix   operator+(const matrix& lh, const matrix& rh);
           void     operator+=(const matrix& rh);
    friend matrix   operator-(const matrix& lh, const matrix& rh);
           void     operator-=(const matrix& rh);
    friend matrix   operator*(const matrix& lh, const matrix& rh);
    friend matrix   operator*(const double& lh, const matrix& rh);
    friend matrix   operator*(const matrix& lh, const double& rh);
           void     operator*=(const double& rh);
    friend matrix   operator/(const matrix& lh, const double& rh);
           void     operator/=(const double& rh);
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
    void set_region(const std::pair<size_t,size_t>& sub, const matrix& mat);
    void insert_rows(const size_t& r, const matrix& mat);
    void insert_cols(const size_t& c, const matrix& mat);
    void remove_rows(const size_t& r, const size_t& n);
    void remove_cols(const size_t& c, const size_t& n);
    matrix resize(const size_t& nrows, const size_t& ncols) const;
    matrix transpose() const;
    bool is_empty() const;
    bool is_square() const;
    bool is_symmetric() const;
    std::vector<std::vector<double>> to_vector() const;

    // matrix instance types
    static matrix identity(const size_t& size);
    static matrix diag(const matrix& vals);
    static matrix empty();
    static matrix from_vector(const std::vector<double>& vec);

    // Gauss-Jordan elimination
    void swap_rows(const size_t& r0, const size_t& r1);
    void swap_cols(const size_t& c0, const size_t& c1);
    void scale_row(const size_t& r, const double& a);
    void axpy_row(const size_t& rx, const size_t& ry, const double& a);
    size_t leading_entry_col(const size_t& r) const;
    void sink_zero_rows();
    std::tuple<matrix, std::vector<std::pair<size_t,size_t>>, size_t> row_reduced_echelon_form() const;
    static matrix solve_GJ(const matrix& A, const matrix& b);
    matrix inverse_GJ() const;
    matrix null_basis() const;

    // determinants
    matrix minor(const size_t& r, const size_t& c) const;
    double determinant() const;

};

/* ************************************************************************* */

/* default ctor */
matrix::matrix() 
{
    m_size = {0, 0};
    m_num_elems = 0;
    m_data = nullptr;
}

/* ctor */
matrix::matrix(const std::pair<size_t,size_t>& size)
{
    m_size = size;
    m_num_elems = size.first * size.second;
    m_data = new double[m_num_elems];
}

/* ctor for matrix of given size with constant value */ 
matrix::matrix(const std::pair<size_t,size_t>& size, const double& val)
{
    m_size = size;
    m_num_elems = size.first * size.second;
    m_data = new double[m_num_elems];

    for (size_t i = 0; i < m_num_elems; ++i)
        m_data[i] = val;
}

/* ctor with nested vector matrix form */ 
matrix::matrix(const std::pair<size_t,size_t>& size, const std::vector<std::vector<double>>& data)
{
    assert(data.size() == size.first);
    m_size = size;
    m_num_elems = m_size.first * m_size.second;
    m_data = new double[m_num_elems];

    for (size_t i = 0; i < m_size.first; ++i)
    {
        assert(data[i].size() == size.second);
        for (size_t j = 0; j < m_size.second; ++j)
            m_data[ind(i, j)] = data[i][j];
    }
}

/* copy ctor */
matrix::matrix(const matrix& rh)
{
    m_size = rh.m_size;
    m_num_elems = rh.m_num_elems;
    m_data = new double[m_num_elems];
    std::memcpy(m_data, rh.m_data, m_num_elems*sizeof(double));
}

/* move ctor (takes temporary rvalue) */
matrix::matrix(matrix&& rh)
{
    m_size = rh.m_size;
    m_num_elems = rh.m_num_elems;
    m_data = rh.m_data;
    rh.m_size = {};
    rh.m_num_elems = 0;
    rh.m_data = nullptr;
}

/* dtor */
matrix::~matrix()
{
    delete[] m_data;
}

/* ************************************************************************* 
Operator overloading */

/* copy assignment operator for lvalue */
matrix& matrix::operator=(const matrix& rh)
{
    // prevent self-assignment
    if (this == &rh)
        return *this;

    // if size is not the same, change the size and array size
    if (m_size != rh.m_size)
    {
        if (m_data != nullptr)
            delete[] m_data;

        m_size = rh.m_size;
        m_num_elems = rh.m_num_elems;
        m_data = new double[m_num_elems];
    }

    // copy values
    std::memcpy(m_data, rh.m_data, m_num_elems*sizeof(double));
    
    return *this;
}

/* move assignment operator for rvalue */
matrix& matrix::operator=(matrix&& rh)
{
    // prevent self-assignment
    if (this == &rh)
        return *this;

    m_size = rh.m_size;
    m_num_elems = rh.m_num_elems;
    m_data = rh.m_data;
    rh.m_size = {};
    rh.m_num_elems = 0;
    rh.m_data = nullptr;

    return *this;
}

/* comparison operator */
bool operator==(const matrix& lh, const matrix& rh)
{
    if (lh.m_size != rh.m_size)
        return false;
    for (size_t i = 0; i < lh.m_num_elems; ++i)
        if (lh.m_data[i] != rh.m_data[i])
            return false;
    return true;
}

/* matrix + matrix */
matrix operator+(const matrix& lh, const matrix& rh)
{
    assert(lh.m_size == rh.m_size);
    matrix result(lh.m_size);
    for (size_t i = 0; i < result.m_num_elems; ++i)
        result.m_data[i] = lh.m_data[i] + rh.m_data[i];    
    return result;
}

/* */
void matrix::operator+=(const matrix& rh)
{
    assert(m_size == rh.m_size); 
    for (size_t i = 0; i < m_size.first; ++i)
        for (size_t j = 0; j < m_size.second; ++j)
            m_data[ind(i, j)] += rh.m_data[ind(i, j)];
}

/* matrix - matrix */
matrix operator-(const matrix& lh, const matrix& rh)
{
    assert(lh.m_size == rh.m_size);
    matrix result(lh.m_size);
    for (size_t i = 0; i < lh.m_num_elems; ++i)
        result.m_data[i] = lh.m_data[i] - rh.m_data[i];    
    return result;
}

/* */
void matrix::operator-=(const matrix& rh)
{
    assert(m_size == rh.m_size); 
    for (size_t i = 0; i < m_size.first; ++i)
        for (size_t j = 0; j < m_size.second; ++j)
            m_data[ind(i, j)] -= rh.m_data[ind(i, j)];
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
    for (size_t i = 0; i < rh.m_num_elems; ++i)
        result.m_data[i] *= lh;
    return result;
}

/* matrix * scalar */
matrix operator*(const matrix& lh, const double& rh)
{
    matrix result = lh;
    for (size_t i = 0; i < lh.m_num_elems; ++i)
        result.m_data[i] *= rh;
    return result;
}

/* */
void matrix::operator*=(const double& rh)
{
    for (size_t i = 0; i < m_num_elems; ++i)
        m_data[i] *= rh;
}

/* matrix / scalar */
matrix operator/(const matrix& lh, const double& rh)
{
    matrix result = lh;
    for (size_t i = 0; i < lh.m_num_elems; ++i)
        result.m_data[i] /= rh;
    return result;
}

/* */
void matrix::operator/=(const double& rh)
{
    for (size_t i = 0; i < m_num_elems; ++i)
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
    for (size_t i = 0; i < m_num_elems; ++i)
        m_data[i] = func();
}

/* lambda operation for all values in matrix taking value as param */
void matrix::lambda(const std::function<double(double)>& func)
{
    for (size_t i = 0; i < m_num_elems; ++i)
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

void matrix::print() const
{
    for (size_t i = 0; i < m_size.first; ++i)
    {
        for (size_t j = 0; j < m_size.second; ++j)
        {
            std::cout.width(10);
            std::cout.precision(3);
            std::cout << m_data[ind(i, j)] << "\t";
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
    assert(mat.m_size.first > 0);
    
    // construct new matrix in tmp 
    double* tmp = new double[m_num_elems + mat.m_num_elems];
    
    // copy rows before insertion
    if (r > 0)
        std::memcpy(tmp, m_data, r*m_size.second*sizeof(double));

    // copy inserted rows
    std::memcpy(tmp + r*m_size.second, mat.m_data,
        mat.m_num_elems*sizeof(double));

    // copy rows after insertion
    if (m_size.first - r > 0)
    {
        std::memcpy(tmp + (r + mat.m_size.first)*m_size.second, m_data + r*m_size.second,
            (m_size.first - r)*m_size.second*sizeof(double));
    }
    // set new size
    m_size.first += mat.m_size.first;

    // point array data to tmp
    delete[] m_data;
    m_data = tmp;
}

/* insert cols */
void matrix::insert_cols(const size_t& c, const matrix& mat)
{
    assert(c <= m_size.second);
    assert(mat.m_size.first == m_size.first);
    assert(mat.m_size.second > 0);

    // construct new matrix in tmp
    double* tmp = new double[m_num_elems + mat.m_num_elems];

    // copy current col num
    size_t ncols = m_size.second;

    // set new size
    m_size.second += mat.m_size.second;

    // allocate values in tmp matrix
    for (size_t i = 0; i < m_size.first; ++i)
        for (size_t j = 0; j < m_size.second; ++j)
            if (j < c)
                tmp[ind(i,j)] = m_data[ncols*i + j];
            else if (j >= c && j < mat.m_size.second + c)
                tmp[ind(i,j)] = mat.m_data[mat.ind(i, j - c)];
            else
                tmp[ind(i,j)] = m_data[ncols*i + (j - mat.m_size.second)];
    
    // point array data to tmp
    delete[] m_data;
    m_data = tmp;
}

/* remove rows */
void matrix::remove_rows(const size_t& r, const size_t& n)
{
    assert(r + n <= m_size.first);
    assert(n > 0);

    // construct new matrix in tmp
    double* tmp = new double[m_num_elems - n*m_size.second];

    // copy rows before removal
    if (r > 0)
        std::memcpy(tmp, m_data, r*m_size.second*sizeof(double));

    // copy rows after removal
    if (m_size.first - r - n > 0)
        std::memcpy(tmp + r*m_size.second, m_data + (r+n)*m_size.second,
            (m_size.first - r - n)*m_size.second*sizeof(double));

    // set new size
    m_size.first -= n;

    // point array data to tmp
    delete[] m_data;
    m_data = tmp;
}

/* remove cols */
void matrix::remove_cols(const size_t& c, const size_t& n)
{
    assert(c + n <= m_size.second);
    assert(n > 0);

    // construct new matrix in tmp
    double* tmp = new double[m_num_elems - n*m_size.first];

    // copy current col num
    size_t ncols = m_size.second;

    // set new size
    m_size.second -= n;

    // allocate values in tmp matrix
    for (size_t i = 0; i < m_size.first; ++i)
        for (size_t j = 0; j < m_size.second; ++j)
            if (j < c)
                tmp[ind(i,j)] = m_data[ncols*i + j];
            else
                tmp[ind(i,j)] = m_data[ncols*i + j + n];

    // point array data to tmp
    delete[] m_data;
    m_data = tmp;
}

/* returns resized matrix */
matrix matrix::resize(const size_t& nrows, const size_t& ncols) const
{
    assert(nrows * ncols == m_num_elems);

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

/* checks if matrix is empty */
bool matrix::is_empty() const
{
    return (m_size.first == 0 && m_size.second == 0 && m_num_elems == 0);
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
    result.m_num_elems = vec.size();
    result.m_data = new double[result.m_num_elems];
    for (size_t i = 0; i < result.m_num_elems; ++i)
        result.m_data[i] = vec[i];
    return result;
}

/* ************************************************************************* 
Gauss-Jordan elimination */

/* swap two rows */
void matrix::swap_rows(const size_t& r0, const size_t& r1) 
{
    double tmp;
    for (size_t i = 0; i < m_size.second; ++i)
    {
        tmp = m_data[ind(r0, i)];
        m_data[ind(r0, i)] = m_data[ind(r1, i)];
        m_data[ind(r1, i)] = tmp;
    }
}

/* swap two cols */
void matrix::swap_cols(const size_t& c0, const size_t& c1) 
{
    double tmp;
    for (size_t i = 0; i < m_size.first; ++i)
    {
        tmp = m_data[ind(i, c0)];
        m_data[ind(i, c0)] = m_data[ind(i, c1)];
        m_data[ind(i, c1)] = tmp;
    }
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

/* send all zero rows to bottom of matrix */
void matrix::sink_zero_rows()
{
    size_t zero_targ = m_size.first - 1;
    for (size_t i = 0; i < m_size.first; ++i)
    {
        if (i >= zero_targ)
            break;
        
        bool zero1 = true;
        for (size_t j = 0; j < m_size.second; ++j)
            if (m_data[ind(i, j)] != 0.0)
            {
                zero1 = false;
                break;
            }
        if (zero1 == true && zero_targ > i)
        {
            // adjust zero targ if the row to be swapped with is also a zero row
            while (zero_targ > i)
            {
                bool zero2 = true;
                for (size_t j = 0; j < m_size.second; ++j)
                    if (m_data[ind(zero_targ, j)] != 0.0)
                    {
                        zero2 = false;
                        break;
                    }
                if (zero2 == true)
                    zero_targ--;
                else
                    break;
            }
            swap_rows(i, zero_targ);
            zero_targ--;
        }
    }
}

/* gauss-jordan elimination to get row-reduced echelon form 
    - returns:  tuple of reduced matrix, pivot subs and type of solution available
    
    - summary of function:
        1. moves rows with all zeros to bottom of matrix
        2. identify first column that isn't all zeros and bring first leading row in this column to top
        3. perform echelon algorithm
        4. identify type of solution available assuming if matrix is in augmented form

    - type of solution available is returned as enum:
        0 = inconsistent system
        1 = unique solution available
        2 = infinite solutions
*/
std::tuple<matrix, std::vector<std::pair<size_t,size_t>>, size_t> matrix::row_reduced_echelon_form() const
{
    size_t solutionType = 1;
    matrix mat(*this);

    // send all rows to bottom of matrix
    sink_zero_rows();

    // locate pivots
    // loop down each row for each col to locate pivots and place them at the top most 'target row', which increments each time a pivot is found
    // then loop down each row in the pivot's column and axpy each row to make sure there are zeros below each pivot in the pivot's column
    // then increment the target col to find the next pivot
    size_t targetCol = 0;
    size_t targetRow = 0;
    std::vector<std::pair<size_t,size_t>> pivots = {};

    while ((targetRow < mat.m_size.first) && (targetCol < mat.m_size.second))
    {
        bool found = false;
        for (size_t i = targetRow; i < mat.m_size.first; ++i)
        {
            if (mat.m_data[mat.ind(i, targetCol)] != 0.0)
            {
                found = true;
                mat.swap_rows(targetRow, i); // once finding the row with the non zero term in the target column, swap this to be at the target row
                mat.scale_row(targetRow, 1.0 / mat.m_data[mat.ind(targetRow, targetCol)]); // then scale the target row s.t. the target element is 1
                pivots.push_back(std::pair<size_t,size_t>({targetRow, targetCol})); // this is a pivot so add it to the pivot vector
                if (targetCol >= mat.m_size.first) // this condition might need to be corrected... assumes that any col beyond the row size (assuming matrix is square) shouldn't have a pivot for consistency
                {
                    solutionType = 0; // pivot located outside on LH square region of size_terest
                }
                break;
            }
        }

        // check if a pivot was found in this column
        if (found == true)
        {
            for (size_t i = targetRow + 1; i < mat.m_size.first; ++i)
            {
                if (mat.m_data[mat.ind(i, targetCol)] != 0.0)
                {
                    mat.axpy_row(i, targetRow, -mat.m_data[mat.ind(i, targetCol)] / mat.m_data[mat.ind(targetRow, targetCol)]); // make sure all elements below the target row in target col are zero
                }
            }
            targetRow++;
        }
        targetCol++;
    }

    // ensure all zeros above each pivot
    for (size_t i = pivots.size() - 1; i >= 0; --i)
    {
        for (size_t j = 0; j < pivots[i].first; ++j)
        {
            if (mat.m_data[mat.ind(j, pivots[i].second)] != 0.0)
            {
                mat.axpy_row(j, pivots[i].first, -mat.m_data[mat.ind(j, pivots[i].second)] / mat.m_data[mat.ind(pivots[i].first, pivots[i].second)]);
            }
        }
    }

    if ((solutionType != 0) && (pivots.size() < mat.m_size.first)) // check that sufficient pivots are found (= number of rows) otherwise there are infinite solutions 
    {
        solutionType = 2;
    }

    return std::tuple<matrix, std::vector<std::pair<size_t,size_t>>, size_t>(mat, pivots, solutionType);
}

/* solve system of equations 
    requires that A is square
*/
matrix matrix::solve_GJ(const matrix& A, const matrix& b)
{
    matrix augMat(A);
    augMat.insert_cols(A.m_size.second, b);

    std::tuple<matrix, std::vector<std::pair<size_t,size_t>>, size_t> tup = augMat.row_reduced_echelon_form();
    matrix reducedAugMat = std::get<0>(tup);
    size_t solutionType = std::get<2>(tup);

    if (solutionType == 1)
        // if there is a single solution then return this
        return reducedAugMat.get_region({0, reducedAugMat.m_size.second - 1}, {reducedAugMat.m_size.first, 1});
    else
        // otherwise return empty matrix
        return matrix::empty();
}

/* inverts a matrix by augmenting on right with identity and reducing to row echelon form thereby obtaining identity on the left and inverse on the right */
matrix matrix::inverse_GJ() const
{ 
    if (m_size.first != m_size.second)
        return matrix::empty(); // if not a square matrix then return empty matrix as obviously not invertible
    
    matrix augMat(*this);
    augMat.insert_cols(m_size.second, matrix::identity(m_size.first));

    std::tuple<matrix, std::vector<std::pair<size_t,size_t>>, size_t> tup = augMat.row_reduced_echelon_form();
    matrix reducedAugMat = std::get<0>(tup);
    size_t solutionType = std::get<2>(tup);

    if (solutionType == 1) // correct number of pivots found, all are located in the LH square region of size_terest
    {
        return reducedAugMat.get_region({0, m_size.first}, {m_size.first, m_size.first});
    }
    else
        return matrix::empty();
}

/* finds the null space of a matrix and returns its basis vectors 
    - obtains the row reduced form
    - then rearranges columns (and thus variables) to get it in the form
    { { I_[r * r] , C_[r * (k-r)] }, { 0_[(n-r) * r] , 0_[(n-r) * (k-r)] } }
        with dimension [n * k]
        where:
            -> r is the rank of this matrix
            -> k is the number of cols (num of variables)
            -> n is the number of rows (num of equations)
            -> I is identity
            -> C is arbitrary matrix
            -> 0 is zero matrix
            -> [r * c] is the dimension of each matrix with r = rows, c = cols

    - the matrix of null-space basis vectors then has the form
    { { -C_[r * (k-r)] } , { I_[(k-r) * (k-r)] } }
        with dimension [k * (k-r)]
    this means that multiplying this matrix with the row-reduced rearranged matrix, you get a zero matrix of dimension [n * (k - r)] as expected for k-r null-space basis vectors
    then the rows of the null space basis vectors are rearranged according to the rearrangement map to get the identity on the top LH of the row reduced matrix
    this ensures that the variables are in the original order provided by the cols of the matrix
*/
matrix matrix::null_basis() const
{
    // get the row echelon form of this matrix
    std::tuple<matrix, std::vector<std::pair<size_t,size_t>>, size_t> tup = row_reduced_echelon_form();
    matrix mat = std::get<0>(tup);
    std::vector<std::pair<size_t,size_t>> pivots = std::get<1>(tup);
    size_t solutionType = std::get<2>(tup);
    size_t rank = pivots.size();

    if (rank == m_size.second)
        return matrix::empty(); // there are same number of pivots as variables therefore nullspace is { 0.0 }

    // create a mapping for variable rearrangement to ensure pivots are contained in a identity matrix on top LH corner of matrix
    std::vector<size_t> rearrangeMap = {};
    for (size_t i = 0; i < m_size.second; ++i)
        rearrangeMap.push_back(i);
    
    // go through pivot positions to determine if swaps need to be made
    size_t targetCol = 0;
    for (size_t i = 0; i < rank; ++i)
    {
        if (pivots[i].second != targetCol)
        {
            // need to swap columns, record this in the rearragement mapping
            mat.swap_cols(pivots[i].second, targetCol);
            rearrangeMap[pivots[i].second] = targetCol;
            rearrangeMap[targetCol] = pivots[i].second;
        }
        targetCol++;
    }
    
    // construct null basis matrix
    matrix nb({m_size.second, m_size.second - rank}); // create matrix for null basis vectors
    nb.set_region({0, 0}, -1.0 * mat.get_region({0, rank}, {rank, m_size.second - rank})); // adds the -C arbitrary matrix
    nb.set_region({rank, 0}, matrix::identity(m_size.second - rank)); // adds the identity 

    // rearrange rows to get original ordering of variables
    for (size_t i = 0; i < m_size.second; ++i)
    {
        if (rearrangeMap[i] != i)
        {
            nb.swap_rows(rearrangeMap[i], i);
            rearrangeMap[rearrangeMap[i]] = rearrangeMap[i];
            rearrangeMap[i] = i;
        }
    }

    return nb;
}

/* ************************************************************************* 
Determinants */

/* returns the minor of the current matrix */
matrix matrix::minor(const size_t& r, const size_t& c) const
{
    matrix result(*this);
    result.remove_rows(r, 1);
    result.remove_cols(c, 1);
    return result;
}

/* calculates the determinant using recursive formula - note this is not computationally efficient, just a proof of concept */
double matrix::determinant() const
{
    if ((m_size.first == 1) && (m_size.second == 1))
    {
        return m_data[0];
    }
    else
    {
        double det = 0.0;
        double sign = 1.0;
        for (size_t i = 0; i < m_size.second; ++i)
        {
            matrix min = minor(0, i);
            det += sign * m_data[ind(0, i)] * min.determinant();
            sign *= -1.0;  
        }
        return det;
    }
}

} // namespace gv
