/* Linear Algebra library 
- Matrix class is higher-performance than v1
*/
#pragma once

#include <vector>
#include <iostream>
#include <assert.h>
#include <functional>

/* **************************************************************************************************
    MATRIX
************************************************************************************************** */
namespace mathlib
{

class Matrix
{
private:
    std::vector<unsigned int> m_size;
    unsigned int m_numElems;
    double* m_valArr;
public:
    Matrix();
    Matrix(const std::vector<unsigned int>& size);
    Matrix(const std::vector<unsigned int>& size, const double& val);
    Matrix(const std::vector<unsigned int>& size, const std::vector<std::vector<double>>& mat); // accepts input matrix in vector-vector form
    Matrix(const Matrix& rh); // copy ctor
    Matrix(Matrix&& rh); // move ctor
    ~Matrix(); // dtor

    // operator overloading
    Matrix& operator=(const Matrix& rh); // copy assignment operator
    Matrix& operator=(Matrix&& rh); // move assignment operator
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
    
    // customised operations with lambda expressions
    Matrix operation(const std::function<double(double)>& func);

    // matrix manipulation
    void display() const;
    std::vector<unsigned int> size() const;
    unsigned int ind(const unsigned int& r, const unsigned int& c) const;
    double get(const std::vector<unsigned int>& sub) const;
    void set(const std::vector<unsigned int>& sub, const double& val);
    Matrix getRegion(const std::vector<unsigned int>& sub, const std::vector<unsigned int>& size) const;
    void setRegion(const std::vector<unsigned int>& sub, const Matrix& mat);
    void insertRows(const unsigned int& rowSub, const Matrix& mat);
    void insertCols(const unsigned int& colSub, const Matrix& mat);
    void removeRows(const unsigned int& rowSub, const unsigned int& numRows);
    void removeCols(const unsigned int& colSub, const unsigned int& numCols);
    Matrix transpose();

    // Gauss-Jordan elimination


};

/* default ctor */
Matrix::Matrix() 
{
    //std::cout << "Default ctor" << std::endl;
    this->m_size = {0, 0};
    this->m_numElems = 0;
    this->m_valArr = nullptr;
}

/* ctor - uninitialised matrix with given size */ 
Matrix::Matrix(const std::vector<unsigned int>& size)
{
    //std::cout << "Ctor 1" << std::endl;
    assert(size.size() == 2); // ensure size has correct form
    this->m_size = size;
    this->m_numElems = size[0] * size[1];
    this->m_valArr = new double[this->m_numElems];//std::vector<double>(size[0] * size[1]);
}

/* */ 
Matrix::Matrix(const std::vector<unsigned int>& size, const double& val)
{
    //std::cout << "Ctor 2" << std::endl;
    assert(size.size() == 2); // ensure size has correct form
    this->m_size = size;
    this->m_numElems = size[0] * size[1];
    this->m_valArr = new double[this->m_numElems]; // std::vector<double>(size[0] * size[1], val);
    for (int i = 0; i < this->m_numElems; ++i)
    {
        this->m_valArr[i] = val;
    }
}

/* ctor accepting matrix argument in vector-vector form */ 
Matrix::Matrix(const std::vector<unsigned int>& size, const std::vector<std::vector<double>>& mat)
{
    //std::cout << "Ctor 3" << std::endl;
    assert(size.size() == 2); // ensure size has correct form
    assert(size[0] == mat.size()); // ensure mat num rows is consistent with size
    this->m_size = size;
    this->m_numElems = size[0] * size[1];
    this->m_valArr = new double[this->m_numElems]; //std::vector<double>(size[0] * size[1]);
    for (unsigned int i = 0; i < this->m_size[0]; ++i) // loop through rows
    {
        assert(mat[i].size() == size[1]); // ensure num cols is consistent for each row
        for (unsigned int j = 0; j < this->m_size[1]; ++j) // loop through cols
        {
            this->m_valArr[this->ind(i, j)] = mat[i][j];
        }
    }
}

/* copy ctor */
Matrix::Matrix(const Matrix& rh)
{
    //std::cout << "Copy ctor" << std::endl;
    this->m_size = rh.m_size;
    this->m_numElems = rh.m_numElems;
    this->m_valArr = new double[this->m_numElems];
    for (int i = 0; i < this->m_numElems; ++i)
    {
        this->m_valArr[i] = rh.m_valArr[i];
    }
}

/* move ctor (takes temporary rvalue) */
Matrix::Matrix(Matrix&& rh)
{
    //std::cout << "Move ctor" << std::endl;
    this->m_size = rh.m_size;
    this->m_numElems = rh.m_numElems;
    this->m_valArr = rh.m_valArr;
    rh.m_size = {};
    rh.m_numElems = 0;
    rh.m_valArr = nullptr;
}

/* dtor */
Matrix::~Matrix()
{
    //std::cout << "Dtor" << std::endl;
    delete[] this->m_valArr;
}

/* ************************************************************************* 
Operator overloading */

/* copy assignment operator for lvalue */
Matrix& Matrix::operator=(const Matrix& rh)
{
    //std::cout << "copy assignment operator" << std::endl;

    // prevent self-assignment
    if (this == &rh)
    {
        return *this;
    }

    // if size is not the same, change the size and array size
    if (this->m_size != rh.m_size)
    {
        if (this->m_valArr != nullptr)
        {
            delete[] this->m_valArr;
        }
        this->m_size = {};
        this->m_numElems = 0;
        this->m_valArr = new double[this->m_numElems];
        this->m_size = rh.m_size;
        this->m_numElems = rh.m_numElems;
    }

    // copy values
    for (int i = 0; i < this->m_numElems; ++i)
    {
        this->m_valArr[i] = rh.m_valArr[i];
    }
    
    return *this;
}

/* move assignment operator for rvalue */
Matrix& Matrix::operator=(Matrix&& rh)
{
    //std::cout << "move assignment operator" << std::endl;

    // prevent self-assignment
    if (this == &rh)
    {
        return *this;
    }

    this->m_size = rh.m_size;
    this->m_numElems = rh.m_numElems;
    this->m_valArr = rh.m_valArr;
    rh.m_size = {};
    rh.m_numElems = 0;
    rh.m_valArr = nullptr;

    return *this;
}

/* comparison operator */
bool operator==(const Matrix& lh, const Matrix& rh)
{
    return ((lh.m_size == rh.m_size) && (lh.m_valArr == rh.m_valArr));
}

/* matrix + matrix */
Matrix operator+(const Matrix& lh, const Matrix& rh)
{
    //std::cout << "+ operator" << std::endl;
    assert(lh.m_size == rh.m_size);
    Matrix result(lh.m_size);
    for (unsigned int i = 0; i < result.m_numElems; ++i) // loop through elements
    {
        result.m_valArr[i] = lh.m_valArr[i] + rh.m_valArr[i];    
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
            this->m_valArr[this->ind(i, j)] += rh.m_valArr[rh.ind(i, j)];
        }
    }
}

/* matrix - matrix */
Matrix operator-(const Matrix& lh, const Matrix& rh)
{
    assert(lh.m_size == rh.m_size);
    Matrix result(lh.m_size);
    for (unsigned int i = 0; i < lh.m_numElems; ++i) // loop through elements
    {
        result.m_valArr[i] = lh.m_valArr[i] - rh.m_valArr[i];    
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
            this->m_valArr[this->ind(i, j)] -= rh.m_valArr[rh.ind(i, j)];
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
                result.m_valArr[result.ind(i, j)] += lh.m_valArr[lh.ind(i, k)] * rh.m_valArr[rh.ind(k, j)];
            }
        }
    }
    return result;
}

/* scalar * matrix */
Matrix operator*(const double& lh, const Matrix& rh)
{
    Matrix result = rh;
    for (unsigned int i = 0; i < rh.m_numElems; ++i)
    {
        result.m_valArr[i] *= lh;
    }
    return result;
}

/* matrix * scalar */
Matrix operator*(const Matrix& lh, const double& rh)
{
    Matrix result = lh;
    for (unsigned int i = 0; i < lh.m_numElems; ++i)
    {
        result.m_valArr[i] *= rh;
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
            this->m_valArr[this->ind(i, j)] *= rh;
        }
    }
}

/* matrix / scalar */
Matrix operator/(const Matrix& lh, const double& rh)
{
    Matrix result = lh;
    for (unsigned int i = 0; i < lh.m_numElems; ++i)
    {
        result.m_valArr[i] /= rh;
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
            this->m_valArr[this->ind(i, j)] /= rh;
        }
    }
}

/* Customised operation with lambda expression of form double func(double) */
Matrix Matrix::operation(const std::function<double(double)>& func)
{
    Matrix result(*this);
    for (int i = 0; i < this->m_numElems; ++i)
    {
        result.m_valArr[i] = func(this->m_valArr[i]);
    }
    return result;
}

/* ************************************************************************* 
Matrix manipulation */

void Matrix::display() const
{
    for (unsigned int i = 0; i < this->m_size[0]; ++i) // loop through rows
    {
        for (unsigned int j = 0; j < this->m_size[1]; ++j) // loop through cols
        {
            std::cout << this->m_valArr[this->ind(i, j)] << "\t";
        }
        std::cout << "\n";
    }
}

/* returns index in valArr of element at position (r, c) - in row-ordered matrix storage */
unsigned int Matrix::ind(const unsigned int& r, const unsigned int& c) const
{
    return r * this->m_size[1] + c;
}

/* returns size of matrix */
std::vector<unsigned int> Matrix::size() const
{
    return this->m_size;
}

/* returns element at sub */
double Matrix::get(const std::vector<unsigned int>& sub) const
{
    //assert((sub.size() == 2) && (sub[0] < this->m_size[0]) && (sub[1] < this->m_size[1])); // ensure sub has correct form
    return this->m_valArr[this->ind(sub[0], sub[1])];
}

/* sets element at sub */
void Matrix::set(const std::vector<unsigned int>& sub, const double& val)
{
    //assert((sub.size() == 2) && (sub[0] < this->m_size[0]) && (sub[1] < this->m_size[1])); // ensure sub has correct form
    this->m_valArr[this->ind(sub[0], sub[1])] = val;
}

/* returns matrix of region starting at sub with given size */
Matrix Matrix::getRegion(const std::vector<unsigned int>& sub, const std::vector<unsigned int>& size) const
{
    Matrix result(size);
    for (unsigned int i = 0; i < size[0]; ++i)
    {
        for (unsigned int j = 0; j < size[1]; ++j)
        {
            result.m_valArr[result.ind(i, j)] = this->m_valArr[this->ind(sub[0] + i, sub[1] + j)];
        }
    }
    return result;
}

/* sets region starting at sub with given size and given matrix */
void Matrix::setRegion(const std::vector<unsigned int>& sub, const Matrix& mat) 
{
    for (unsigned int i = 0; i < mat.m_size[0]; ++i)
    {
        for (unsigned int j = 0; j < mat.m_size[1]; ++j)
        {
            this->m_valArr[mat.ind(sub[0] + i, sub[1] + j)] = mat.m_valArr[mat.ind(i, j)];
        }
    }
}

/* insert rows */
void Matrix::insertRows(const unsigned int& rowSub, const Matrix& mat)
{
    assert(rowSub <= this->m_size[0]);
    assert(mat.m_size[1] == this->m_size[1]);

    // copy current array
    double* tmp = new double[this->m_numElems];
    for (unsigned int i = 0; i < this->m_numElems; ++i)
    {
        tmp[i] = this->m_valArr[i];
    }

    // increase array size
    if (this->m_valArr != nullptr)
    {
        delete[] this->m_valArr;
    }
    std::vector<unsigned int> tmpSize = this->m_size; // copy prev size
    this->m_size[0] += mat.m_size[0];
    this->m_numElems = this->m_size[0] * this->m_size[1];
    this->m_valArr = new double[this->m_numElems];

    // fill new array with corresponding values
    for (unsigned int i = 0; i < this->m_size[0]; ++i) // loop through rows
    {
        for (unsigned int j = 0; j < this->m_size[1]; ++j) // loop through cols
        {
            if (i < rowSub)
            {
                this->m_valArr[this->ind(i, j)] = tmp[i * (tmpSize[1]) + j];
            }
            if ((i >= rowSub) && (i < rowSub + mat.m_size[0]))
            {
                this->m_valArr[this->ind(i, j)] = mat.m_valArr[mat.ind(i - rowSub, j)]; // insert row from mat
            }
            else if (i >= rowSub + mat.m_size[0])
            {
                this->m_valArr[this->ind(i, j)] = tmp[(i - mat.m_size[0])*(tmpSize[1]) + j];
            }
        }
    }
    delete[] tmp;
}

/* insert cols */
void Matrix::insertCols(const unsigned int& colSub, const Matrix& mat)
{
    assert(colSub <= this->m_size[1]);
    assert(mat.m_size[0] == this->m_size[0]);

    // copy current array
    double* tmp = new double[this->m_numElems];
    for (unsigned int i = 0; i < this->m_numElems; ++i)
    {
        tmp[i] = this->m_valArr[i];
    }

    // increase array size
    if (this->m_valArr != nullptr)
    {
        delete[] this->m_valArr;
    }
    std::vector<unsigned int> tmpSize = this->m_size; // copy prev size
    this->m_size[1] += mat.m_size[1];
    this->m_numElems = this->m_size[0] * this->m_size[1];
    this->m_valArr = new double[this->m_numElems];

    // fill new array with corresponding values
    for (unsigned int i = 0; i < this->m_size[0]; ++i) // loop through rows
    {
        for (unsigned int j = 0; j < this->m_size[1]; ++j) // loop through cols
        {
            if (j < colSub)
            {
                this->m_valArr[this->ind(i, j)] = tmp[i * (tmpSize[1]) + j];
            }
            if ((j >= colSub) && (j < colSub + mat.m_size[1]))
            {
                this->m_valArr[this->ind(i, j)] = mat.m_valArr[mat.ind(i, j - colSub)]; // insert col from mat
            }
            else if (j >= colSub + mat.m_size[1])
            {
                this->m_valArr[this->ind(i, j)] = tmp[i * (tmpSize[1]) + j - mat.m_size[1]];
            }
        }
    }
    delete[] tmp;
}

/* remove rows */
void Matrix::removeRows(const unsigned int& rowSub, const unsigned int& numRows)
{
    assert(rowSub + numRows <= this->m_size[0]);
    assert(numRows > 0);

    // copy current array
    double* tmp = new double[this->m_numElems];
    for (unsigned int i = 0; i < this->m_numElems; ++i)
    {
        tmp[i] = this->m_valArr[i];
    }

    // decrease array size
    if (this->m_valArr != nullptr)
    {
        delete[] this->m_valArr;
    }
    std::vector<unsigned int> tmpSize = this->m_size; // copy prev size
    this->m_size[0] -= numRows;
    this->m_numElems = this->m_size[0] * this->m_size[1];
    this->m_valArr = new double[this->m_numElems];

    // fill new array with corresponding values
    for (unsigned int i = 0; i < this->m_size[0]; ++i) // loop through rows
    {
        for (unsigned int j = 0; j < this->m_size[1]; ++j) // loop through cols
        {
            if (i < rowSub)
            {
                this->m_valArr[this->ind(i, j)] = tmp[i * (tmpSize[1]) + j];
            }
            else if (i >= rowSub)
            {
                this->m_valArr[this->ind(i, j)] = tmp[(i + numRows)*(tmpSize[1]) + j];
            }
        }
    }
    delete[] tmp;
}

/* remove cols */
void Matrix::removeCols(const unsigned int& colSub, const unsigned int& numCols)
{
    assert(colSub + numCols <= this->m_size[1]);
    assert(numCols > 0);

    // copy current array
    double* tmp = new double[this->m_numElems];
    for (unsigned int i = 0; i < this->m_numElems; ++i)
    {
        tmp[i] = this->m_valArr[i];
    }

    // decrease array size
    if (this->m_valArr != nullptr)
    {
        delete[] this->m_valArr;
    }
    std::vector<unsigned int> tmpSize = this->m_size; // copy prev size
    this->m_size[1] -= numCols;
    this->m_numElems = this->m_size[0] * this->m_size[1];
    this->m_valArr = new double[this->m_numElems];

    // fill new array with corresponding values
    for (unsigned int i = 0; i < this->m_size[0]; ++i) // loop through rows
    {
        for (unsigned int j = 0; j < this->m_size[1]; ++j) // loop through cols
        {
            if (j < colSub)
            {
                this->m_valArr[this->ind(i, j)] = tmp[i * (tmpSize[1]) + j];
            }
            else if (j >= colSub)
            {
                this->m_valArr[this->ind(i, j)] = tmp[i * (tmpSize[1]) + j + numCols];
            }
        }
    }
    delete[] tmp;
}

/* returns transpose of matrix */
Matrix Matrix::transpose()
{
    Matrix result({this->m_size[1], this->m_size[0]});
    for (unsigned int i = 0; i < this->m_size[0]; ++i)
    {
        for (unsigned int j = 0; j < this->m_size[1]; ++j)
        {
            result.m_valArr[result.ind(j, i)] = this->m_valArr[this->ind(i, j)];
        }
    }
    return result;
}

/* ************************************************************************* 
Gauss-Jordan elimination */



} // namespace mathlib