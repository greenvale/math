/* Linear Algebra library 
- Matrix class is higher-performance than v1
*/
#pragma once

#include <vector>
#include <iostream>
#include <assert.h>

/* **************************************************************************************************
    MATRIX
************************************************************************************************** */
namespace mathlib
{

class Matrix
{
private:
    std::vector<unsigned int> m_size;
    std::vector<double> m_elemArr;
public:
    Matrix();
    Matrix(const std::vector<unsigned int>& size);
    Matrix(const std::vector<unsigned int>& size, const double& val);
    //Matrix(const std::vector<unsigned int>& size, const std::vector<double>& elemArr); - removed to prevent ambiguous error!
    Matrix(const std::vector<unsigned int>& size, const std::vector<std::vector<double>>& mat); // accepts matrix in vector-vector form
    Matrix(const Matrix& rh); // copy ctor

    Matrix& operator=(const Matrix& rh); // copy assignment operator
    friend bool operator==(const Matrix& lh, const Matrix& rh);
    friend Matrix operator+(const Matrix& lh, const Matrix& rh);
    void operator+=(const Matrix& rh);
    friend Matrix operator-(const Matrix& lh, const Matrix& rh);
    void operator-=(const Matrix& rh);
    friend Matrix operator*(const Matrix& lh, const Matrix& rh);
    friend Matrix operator*(const double& lh, const Matrix& rh);
    friend Matrix operator*(const Matrix& lh, const double& rh);
    void operator*=(const Matrix& rh);
    void operator*=(const double& rh);
    friend Matrix operator/(const Matrix& lh, const double& rh);
    void operator/=(const double& rh);

    void display() const;
    std::vector<unsigned int> size() const;
    unsigned int ind(const unsigned int& r, const unsigned int& c) const;
    double get(const std::vector<unsigned int>& sub) const;
    void set(const std::vector<unsigned int>& sub, const double& val);
    Matrix getRegion(const std::vector<unsigned int>& sub, const std::vector<unsigned int>& size) const;
    void setRegion(const std::vector<unsigned int>& sub, const Matrix& mat);
};

/* default ctor */
Matrix::Matrix() 
{
    //std::cout << "Default ctor" << std::endl;
    this->m_size = {0, 0};
    this->m_elemArr = {};
}

/* */ 
Matrix::Matrix(const std::vector<unsigned int>& size)
{
    assert(size.size() == 2); // ensure size has correct form
    this->m_size = size;
    this->m_elemArr = std::vector<double>(size[0] * size[1]);
}

/* */ 
Matrix::Matrix(const std::vector<unsigned int>& size, const double& val)
{
    assert(size.size() == 2); // ensure size has correct form
    this->m_size = size;
    this->m_elemArr = std::vector<double>(size[0] * size[1], val);
}

/* */
/*
Matrix::Matrix(const std::vector<unsigned int>& size, const std::vector<double>& elemArr)
{
    assert(size.size() == 2); // ensure size has correct form
    assert(elemArr.size() == size[0] * size[1]); // ensure elem array has correct form
    this->m_size = size;
    this->m_elemArr = elemArr;
}
*/

/* ctor accepting matrix argument in vector-vector form */ 
Matrix::Matrix(const std::vector<unsigned int>& size, const std::vector<std::vector<double>>& mat)
{
    assert(size.size() == 2); // ensure size has correct form
    assert(size[0] == mat.size()); // ensure mat num rows is consistent with size
    this->m_size = size;
    this->m_elemArr = std::vector<double>(size[0] * size[1]);
    for (unsigned int i = 0; i < size[0]; ++i) // loop through rows
    {
        assert(mat[i].size() == size[1]); // ensure num cols is consistent for each row
        for (unsigned int j = 0; j < size[1]; ++j) // loop through cols
        {
            this->m_elemArr[this->ind(i, j)] = mat[i][j];
        }
    }
}

/* copy ctor */
Matrix::Matrix(const Matrix& rh)
{
    //std::cout << "Copy ctor" << std::endl;
    this->m_size = rh.m_size;
    this->m_elemArr = rh.m_elemArr;
}

/* ************************************************************************* 
Operator overloading */

/* copy assignment operator */
Matrix& Matrix::operator=(const Matrix& rh)
{
    //std::cout << "= operator" << std::endl;
    this->m_size = rh.m_size;
    this->m_elemArr = rh.m_elemArr;
}

/* comparison operator */
bool operator==(const Matrix& lh, const Matrix& rh)
{
    return ((lh.m_size == rh.m_size) && (lh.m_elemArr == rh.m_elemArr));
}

/* matrix + matrix */
Matrix operator+(const Matrix& lh, const Matrix& rh)
{
    //std::cout << "+ operator" << std::endl;
    assert(lh.m_size == rh.m_size);
    Matrix result(lh.m_size);
    for (unsigned int i = 0; i < lh.m_size[0] * lh.m_size[1]; ++i) // loop through elements
    {
        result.m_elemArr[i] = lh.m_elemArr[i] + rh.m_elemArr[i];    
    }
    return result;
}

/* */
void Matrix::operator+=(const Matrix& rh)
{
    assert(this->m_size == rh.m_size); 
    for (unsigned int i = 0; i < this->m_size[0]; ++i)
    {
        for (unsigned int j = 0; j < this->m_size[1]; ++j)
        {
            this->m_elemArr[this->ind(i, j)] += rh.m_elemArr[rh.ind(i, j)];
        }
    }
}

/* matrix - matrix */
Matrix operator-(const Matrix& lh, const Matrix& rh)
{
    assert(lh.m_size == rh.m_size);
    Matrix result(lh.m_size);
    for (unsigned int i = 0; i < lh.m_size[0] * lh.m_size[1]; ++i) // loop through elements
    {
        result.m_elemArr[i] = lh.m_elemArr[i] - rh.m_elemArr[i];    
    }
    return result;
}

/* */
void Matrix::operator-=(const Matrix& rh)
{
    assert(this->m_size == rh.m_size); 
    for (unsigned int i = 0; i < this->m_size[0]; ++i)
    {
        for (unsigned int j = 0; j < this->m_size[1]; ++j)
        {
            this->m_elemArr[this->ind(i, j)] -= rh.m_elemArr[rh.ind(i, j)];
        }
    }
}

/* matrix * matrix */
Matrix operator*(const Matrix& lh, const Matrix& rh)
{
    assert(lh.m_size[1] == rh.m_size[0]); // num cols for lh must equate num rows for rh
    Matrix result({lh.m_size[0], rh.m_size[1]}, 0.0); // result size is (num rows for lh, num cols for rh)
    for (unsigned int i = 0; i < lh.m_size[0]; ++i) // loop through rows in lh
    {
        for (unsigned int j = 0; j < rh.m_size[1]; ++j) // loop through cols in rh
        {
            for (unsigned int k = 0; k < lh.m_size[1]; ++k) // for element (i, j) take dot product of row i in lh and col j in rh
            {
                result.m_elemArr[result.ind(i, j)] += lh.m_elemArr[lh.ind(i, k)] * rh.m_elemArr[rh.ind(k, j)];
            }
        }
    }
    return result;
}

/* scalar * matrix */
Matrix operator*(const double& lh, const Matrix& rh)
{
    Matrix result = rh;
    for (unsigned int i = 0; i < rh.m_size[0] * rh.m_size[1]; ++i)
    {
        result.m_elemArr[i] *= lh;
    }
    return result;
}

/* matrix * scalar */
Matrix operator*(const Matrix& lh, const double& rh)
{
    Matrix result = lh;
    for (unsigned int i = 0; i < lh.m_size[0] * lh.m_size[1]; ++i)
    {
        result.m_elemArr[i] *= rh;
    }
    return result;
}

/* */
void Matrix::operator*=(const Matrix& rh)
{
    Matrix result = (*this) * rh;
    (*this) = result;
}

/* */
void Matrix::operator*=(const double& rh)
{
    for (unsigned int i = 0; i < this->m_size[0]; ++i)
    {
        for (unsigned int j = 0; j < this->m_size[1]; ++j)
        {
            this->m_elemArr[this->ind(i, j)] *= rh;
        }
    }
}

/* matrix / scalar */
Matrix operator/(const Matrix& lh, const double& rh)
{
    Matrix result = lh;
    for (unsigned int i = 0; i < lh.m_size[0] * lh.m_size[1]; ++i)
    {
        result.m_elemArr[i] /= rh;
    }
    return result;
}

/* */
void Matrix::operator/=(const double& rh)
{
    for (unsigned int i = 0; i < this->m_size[0]; ++i)
    {
        for (unsigned int j = 0; j < this->m_size[1]; ++j)
        {
            this->m_elemArr[this->ind(i, j)] /= rh;
        }
    }
}

/* ************************************************************************* 
Matrix manipulation */

void Matrix::display() const
{
    for (unsigned int i = 0; i < this->m_size[0]; ++i) // loop through rows
    {
        for (unsigned int j = 0; j < this->m_size[1]; ++j) // loop through cols
        {
            std::cout << this->m_elemArr[this->ind(i, j)] << "\t";
        }
        std::cout << "\n";
    }
}

/* returns index in elemArr of element at position (r, c) - in row-ordered matrix storage */
unsigned int Matrix::ind(const unsigned int& r, const unsigned int& c) const
{
    return r * this->m_size[1] + c;
}

/* returns size of matrix */
std::vector<unsigned int> Matrix::size() const
{
    return this->m_size;
}

/* */
double Matrix::get(const std::vector<unsigned int>& sub) const
{
    //assert((sub.size() == 2) && (sub[0] < this->m_size[0]) && (sub[1] < this->m_size[1])); // ensure sub has correct form
    return this->m_elemArr[this->ind(sub[0], sub[1])];
}

/* */
void Matrix::set(const std::vector<unsigned int>& sub, const double& val)
{
    //assert((sub.size() == 2) && (sub[0] < this->m_size[0]) && (sub[1] < this->m_size[1])); // ensure sub has correct form
    this->m_elemArr[this->ind(sub[0], sub[1])] = val;
}

/* returns matrix of region starting at sub with given size */
Matrix Matrix::getRegion(const std::vector<unsigned int>& sub, const std::vector<unsigned int>& size) const
{
    Matrix mat(size);
    for (unsigned int i = 0; i < size[0]; ++i)
    {
        for (unsigned int j = 0; j < size[1]; ++j)
        {
            mat.m_elemArr[mat.ind(i, j)] = this->m_elemArr[this->ind(sub[0] + i, sub[1] + j)];
        }
    }
    return mat;
}

/* sets region starting at sub with given size and given matrix */
void Matrix::setRegion(const std::vector<unsigned int>& sub, const Matrix& mat) 
{
    for (unsigned int i = 0; i < mat.m_size[0]; ++i)
    {
        for (unsigned int j = 0; j < mat.m_size[1]; ++j)
        {
            this->m_elemArr[mat.ind(sub[0] + i, sub[1] + j)] = mat.m_elemArr[mat.ind(i, j)];
        }
    }
}

} // namespace mathlib
