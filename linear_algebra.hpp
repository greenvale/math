/* Linear Algebra library by William Denny (greenvale)
- Matrix class is higher-performance than v1
*/
#pragma once

#include <vector>
#include <iostream>
#include <assert.h>
#include <functional>
#include <map>
#include <tuple>

namespace mathlib
{

/* **************************************************************************************************
    MATRIX
************************************************************************************************** */

class Matrix
{
private:
    std::vector<unsigned int> m_size;
    unsigned int m_num_elems;
    double* m_data;
public:
    Matrix();
    Matrix(const std::vector<unsigned int>& size);
    Matrix(const std::vector<unsigned int>& size, const double& val);
    Matrix(const std::vector<unsigned int>& size, const std::vector<std::vector<double>>& mat); // accepts input matrix in vector-vector form
    Matrix(const Matrix& rh); // copy ctor
    Matrix(Matrix&& rh); // move ctor
    ~Matrix(); // dtor

    // operator overloading
    Matrix&         operator=(const Matrix& rh); // copy assignment operator
    Matrix&         operator=(Matrix&& rh); // move assignment operator
    friend bool     operator==(const Matrix& lh, const Matrix& rh);
    friend Matrix   operator+(const Matrix& lh, const Matrix& rh);
    void            operator+=(const Matrix& rh);
    friend Matrix   operator-(const Matrix& lh, const Matrix& rh);
    void            operator-=(const Matrix& rh);
    friend Matrix   operator*(const Matrix& lh, const Matrix& rh);
    friend Matrix   operator*(const double& lh, const Matrix& rh);
    friend Matrix   operator*(const Matrix& lh, const double& rh);
    void            operator*=(const Matrix& rh);
    void            operator*=(const double& rh);
    friend Matrix   operator/(const Matrix& lh, const double& rh);
    void            operator/=(const double& rh);
    double&         operator[](const std::vector<unsigned int>& sub);
    
    // customised uniform operation for all elements given individual function
    void operation(const std::function<double()>& func);
    void operation(const std::function<double(double)>& func);

    // basic matrix manipulation
    void display() const;
    std::vector<unsigned int> size() const;
    unsigned int ind(const unsigned int& r, const unsigned int& c) const;
    double get(const std::vector<unsigned int>& sub) const;
    void set(const std::vector<unsigned int>& sub, const double& val);
    Matrix get_region(const std::vector<unsigned int>& sub, const std::vector<unsigned int>& size) const;
    void set_region(const std::vector<unsigned int>& sub, const Matrix& mat);
    void insert_rows(const unsigned int& rowSub, const Matrix& mat);
    void insert_cols(const unsigned int& colSub, const Matrix& mat);
    void remove_rows(const unsigned int& rowSub, const unsigned int& numRows);
    void remove_cols(const unsigned int& colSub, const unsigned int& numCols);
    Matrix resize(const unsigned int& numRows, const unsigned int& numCols);
    Matrix transpose();
    bool is_empty();
    bool is_square();
    bool is_symmetric();

    // matrix instance types
    static Matrix identity(const unsigned int& size);
    static Matrix diag(const Matrix& vals);
    static Matrix empty();

    // Gauss-Jordan elimination
    void swap_rows(const unsigned int& r0, const unsigned int& r1);
    void swap_cols(const unsigned int& c0, const unsigned int& c1);
    void scale_row(const unsigned int& r, const double& a);
    void axpy_row(const unsigned int& rx, const unsigned int& ry, const double& a);
    unsigned int leading_entry_col(const unsigned int& r);
    std::tuple<Matrix, std::vector<std::vector<unsigned int>>, unsigned int> row_reduced_echelon_form();
    static Matrix solve_GJ(const Matrix& A, const Matrix& b);
    Matrix inverse_GJ();
    Matrix nullBasis();

    // determinants
    Matrix minor(const unsigned int& r, const unsigned int& c);
    double determinant();

};

/* default ctor */
Matrix::Matrix() 
{
    //std::cout << "Default ctor" << std::endl;
    this->m_size = {0, 0};
    this->m_num_elems = 0;
    this->m_data = nullptr;
}

/* ctor - uninitialised matrix with given size */ 
Matrix::Matrix(const std::vector<unsigned int>& size)
{
    //std::cout << "Ctor 1" << std::endl;
    assert(size.size() == 2); // ensure size has correct form
    this->m_size = size;
    this->m_num_elems = size[0] * size[1];
    this->m_data = new double[this->m_num_elems];//std::vector<double>(size[0] * size[1]);
}

/* */ 
Matrix::Matrix(const std::vector<unsigned int>& size, const double& val)
{
    //std::cout << "Ctor 2" << std::endl;
    assert(size.size() == 2); // ensure size has correct form
    this->m_size = size;
    this->m_num_elems = size[0] * size[1];
    this->m_data = new double[this->m_num_elems]; // std::vector<double>(size[0] * size[1], val);
    for (int i = 0; i < this->m_num_elems; ++i)
    {
        this->m_data[i] = val;
    }
}

/* ctor accepting matrix argument in vector-vector form */ 
Matrix::Matrix(const std::vector<unsigned int>& size, const std::vector<std::vector<double>>& mat)
{
    //std::cout << "Ctor 3" << std::endl;
    assert(size.size() == 2); // ensure size has correct form
    assert(size[0] == mat.size()); // ensure mat num rows is consistent with size
    this->m_size = size;
    this->m_num_elems = size[0] * size[1];
    this->m_data = new double[this->m_num_elems]; //std::vector<double>(size[0] * size[1]);
    for (unsigned int i = 0; i < this->m_size[0]; ++i) // loop through rows
    {
        assert(mat[i].size() == size[1]); // ensure num cols is consistent for each row
        for (unsigned int j = 0; j < this->m_size[1]; ++j) // loop through cols
        {
            this->m_data[this->ind(i, j)] = mat[i][j];
        }
    }
}

/* copy ctor */
Matrix::Matrix(const Matrix& rh)
{
    //std::cout << "Copy ctor" << std::endl;
    this->m_size = rh.m_size;
    this->m_num_elems = rh.m_num_elems;
    this->m_data = new double[this->m_num_elems];
    for (int i = 0; i < this->m_num_elems; ++i)
    {
        this->m_data[i] = rh.m_data[i];
    }
}

/* move ctor (takes temporary rvalue) */
Matrix::Matrix(Matrix&& rh)
{
    //std::cout << "Move ctor" << std::endl;
    this->m_size = rh.m_size;
    this->m_num_elems = rh.m_num_elems;
    this->m_data = rh.m_data;
    rh.m_size = {};
    rh.m_num_elems = 0;
    rh.m_data = nullptr;
}

/* dtor */
Matrix::~Matrix()
{
    //std::cout << "Dtor" << std::endl;
    delete[] this->m_data;
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
        if (this->m_data != nullptr)
        {
            delete[] this->m_data;
        }
        this->m_size = {};
        this->m_num_elems = 0;
        this->m_data = new double[this->m_num_elems];
        this->m_size = rh.m_size;
        this->m_num_elems = rh.m_num_elems;
    }

    // copy values
    for (int i = 0; i < this->m_num_elems; ++i)
    {
        this->m_data[i] = rh.m_data[i];
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
    this->m_num_elems = rh.m_num_elems;
    this->m_data = rh.m_data;
    rh.m_size = {};
    rh.m_num_elems = 0;
    rh.m_data = nullptr;

    return *this;
}

/* comparison operator */
bool operator==(const Matrix& lh, const Matrix& rh)
{
    return ((lh.m_size == rh.m_size) && (lh.m_data == rh.m_data));
}

/* matrix + matrix */
Matrix operator+(const Matrix& lh, const Matrix& rh)
{
    //std::cout << "+ operator" << std::endl;
    assert(lh.m_size == rh.m_size);
    Matrix result(lh.m_size);
    for (unsigned int i = 0; i < result.m_num_elems; ++i) // loop through elements
    {
        result.m_data[i] = lh.m_data[i] + rh.m_data[i];    
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
            this->m_data[this->ind(i, j)] += rh.m_data[rh.ind(i, j)];
        }
    }
}

/* matrix - matrix */
Matrix operator-(const Matrix& lh, const Matrix& rh)
{
    assert(lh.m_size == rh.m_size);
    Matrix result(lh.m_size);
    for (unsigned int i = 0; i < lh.m_num_elems; ++i) // loop through elements
    {
        result.m_data[i] = lh.m_data[i] - rh.m_data[i];    
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
            this->m_data[this->ind(i, j)] -= rh.m_data[rh.ind(i, j)];
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
                result.m_data[result.ind(i, j)] += lh.m_data[lh.ind(i, k)] * rh.m_data[rh.ind(k, j)];
            }
        }
    }
    return result;
}

/* scalar * matrix */
Matrix operator*(const double& lh, const Matrix& rh)
{
    Matrix result = rh;
    for (unsigned int i = 0; i < rh.m_num_elems; ++i)
    {
        result.m_data[i] *= lh;
    }
    return result;
}

/* matrix * scalar */
Matrix operator*(const Matrix& lh, const double& rh)
{
    Matrix result = lh;
    for (unsigned int i = 0; i < lh.m_num_elems; ++i)
    {
        result.m_data[i] *= rh;
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
            this->m_data[this->ind(i, j)] *= rh;
        }
    }
}

/* matrix / scalar */
Matrix operator/(const Matrix& lh, const double& rh)
{
    Matrix result = lh;
    for (unsigned int i = 0; i < lh.m_num_elems; ++i)
    {
        result.m_data[i] /= rh;
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
            this->m_data[this->ind(i, j)] /= rh;
        }
    }
}

/* */
double& Matrix::operator[](const std::vector<unsigned int>& sub)
{
    assert(sub.size() == 2);
    return this->m_data[this->ind(sub[0], sub[1])];
}

/* customised uniform operation for all elements given individual function */
void Matrix::operation(const std::function<double()>& func)
{
    for (int i = 0; i < this->m_num_elems; ++i)
    {
        this->m_data[i] = func();
    }
}

/* customised uniform operation for all elements given individual function with 1 parameter (normally element value) */
void Matrix::operation(const std::function<double(double)>& func)
{
    for (int i = 0; i < this->m_num_elems; ++i)
    {
        this->m_data[i] = func(this->m_data[i]);
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
            std::cout << this->m_data[this->ind(i, j)] << "\t";
        }
        std::cout << "\n";
    }
}

/* returns index in valArr of element at position (r, c) - in row-ordered matrix storage */
unsigned int Matrix::ind(const unsigned int& r, const unsigned int& c) const
{
    assert((r < this->m_size[0]) && (c < this->m_size[1]));
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
    return this->m_data[this->ind(sub[0], sub[1])];
}

/* sets element at sub */
void Matrix::set(const std::vector<unsigned int>& sub, const double& val)
{
    //assert((sub.size() == 2) && (sub[0] < this->m_size[0]) && (sub[1] < this->m_size[1])); // ensure sub has correct form
    this->m_data[this->ind(sub[0], sub[1])] = val;
}

/* returns matrix of region starting at sub with given size */
Matrix Matrix::get_region(const std::vector<unsigned int>& sub, const std::vector<unsigned int>& size) const
{
    Matrix result(size);
    for (unsigned int i = 0; i < size[0]; ++i)
    {
        for (unsigned int j = 0; j < size[1]; ++j)
        {
            result.m_data[result.ind(i, j)] = this->m_data[this->ind(sub[0] + i, sub[1] + j)];
        }
    }
    return result;
}

/* sets region starting at sub with given size and given matrix */
void Matrix::set_region(const std::vector<unsigned int>& sub, const Matrix& mat) 
{
    for (unsigned int i = 0; i < mat.m_size[0]; ++i)
    {
        for (unsigned int j = 0; j < mat.m_size[1]; ++j)
        {
            this->m_data[this->ind(sub[0] + i, sub[1] + j)] = mat.m_data[mat.ind(i, j)];
        }
    }
}

/* insert rows */
void Matrix::insert_rows(const unsigned int& rowSub, const Matrix& mat)
{
    assert(rowSub <= this->m_size[0]);
    assert(mat.m_size[1] == this->m_size[1]);

    // copy current array
    double* tmp = new double[this->m_num_elems];
    for (unsigned int i = 0; i < this->m_num_elems; ++i)
    {
        tmp[i] = this->m_data[i];
    }

    // increase array size
    if (this->m_data != nullptr)
    {
        delete[] this->m_data;
    }
    std::vector<unsigned int> tmpSize = this->m_size; // copy prev size
    this->m_size[0] += mat.m_size[0];
    this->m_num_elems = this->m_size[0] * this->m_size[1];
    this->m_data = new double[this->m_num_elems];

    // fill new array with corresponding values
    for (unsigned int i = 0; i < this->m_size[0]; ++i) // loop through rows
    {
        for (unsigned int j = 0; j < this->m_size[1]; ++j) // loop through cols
        {
            if (i < rowSub)
            {
                this->m_data[this->ind(i, j)] = tmp[i * (tmpSize[1]) + j];
            }
            if ((i >= rowSub) && (i < rowSub + mat.m_size[0]))
            {
                this->m_data[this->ind(i, j)] = mat.m_data[mat.ind(i - rowSub, j)]; // insert row from mat
            }
            else if (i >= rowSub + mat.m_size[0])
            {
                this->m_data[this->ind(i, j)] = tmp[(i - mat.m_size[0])*(tmpSize[1]) + j];
            }
        }
    }
    delete[] tmp;
}

/* insert cols */
void Matrix::insert_cols(const unsigned int& colSub, const Matrix& mat)
{
    assert(colSub <= this->m_size[1]);
    assert(mat.m_size[0] == this->m_size[0]);

    // copy current array
    double* tmp = new double[this->m_num_elems];
    for (unsigned int i = 0; i < this->m_num_elems; ++i)
    {
        tmp[i] = this->m_data[i];
    }

    // increase array size
    if (this->m_data != nullptr)
    {
        delete[] this->m_data;
    }
    std::vector<unsigned int> tmpSize = this->m_size; // copy prev size
    this->m_size[1] += mat.m_size[1];
    this->m_num_elems = this->m_size[0] * this->m_size[1];
    this->m_data = new double[this->m_num_elems];

    // fill new array with corresponding values
    for (unsigned int i = 0; i < this->m_size[0]; ++i) // loop through rows
    {
        for (unsigned int j = 0; j < this->m_size[1]; ++j) // loop through cols
        {
            if (j < colSub)
            {
                this->m_data[this->ind(i, j)] = tmp[i * (tmpSize[1]) + j];
            }
            if ((j >= colSub) && (j < colSub + mat.m_size[1]))
            {
                this->m_data[this->ind(i, j)] = mat.m_data[mat.ind(i, j - colSub)]; // insert col from mat
            }
            else if (j >= colSub + mat.m_size[1])
            {
                this->m_data[this->ind(i, j)] = tmp[i * (tmpSize[1]) + j - mat.m_size[1]];
            }
        }
    }
    delete[] tmp;
}

/* remove rows */
void Matrix::remove_rows(const unsigned int& rowSub, const unsigned int& numRows)
{
    assert(rowSub + numRows <= this->m_size[0]);
    assert(numRows > 0);

    // copy current array
    double* tmp = new double[this->m_num_elems];
    for (unsigned int i = 0; i < this->m_num_elems; ++i)
    {
        tmp[i] = this->m_data[i];
    }

    // decrease array size
    if (this->m_data != nullptr)
    {
        delete[] this->m_data;
    }
    std::vector<unsigned int> tmpSize = this->m_size; // copy prev size
    this->m_size[0] -= numRows;
    this->m_num_elems = this->m_size[0] * this->m_size[1];
    this->m_data = new double[this->m_num_elems];

    // fill new array with corresponding values
    for (unsigned int i = 0; i < this->m_size[0]; ++i) // loop through rows
    {
        for (unsigned int j = 0; j < this->m_size[1]; ++j) // loop through cols
        {
            if (i < rowSub)
            {
                this->m_data[this->ind(i, j)] = tmp[i * (tmpSize[1]) + j];
            }
            else if (i >= rowSub)
            {
                this->m_data[this->ind(i, j)] = tmp[(i + numRows)*(tmpSize[1]) + j];
            }
        }
    }
    delete[] tmp;
}

/* remove cols */
void Matrix::remove_cols(const unsigned int& colSub, const unsigned int& numCols)
{
    assert(colSub + numCols <= this->m_size[1]);
    assert(numCols > 0);

    // copy current array
    double* tmp = new double[this->m_num_elems];
    for (unsigned int i = 0; i < this->m_num_elems; ++i)
    {
        tmp[i] = this->m_data[i];
    }

    // decrease array size
    if (this->m_data != nullptr)
    {
        delete[] this->m_data;
    }
    std::vector<unsigned int> tmpSize = this->m_size; // copy prev size
    this->m_size[1] -= numCols;
    this->m_num_elems = this->m_size[0] * this->m_size[1];
    this->m_data = new double[this->m_num_elems];

    // fill new array with corresponding values
    for (unsigned int i = 0; i < this->m_size[0]; ++i) // loop through rows
    {
        for (unsigned int j = 0; j < this->m_size[1]; ++j) // loop through cols
        {
            if (j < colSub)
            {
                this->m_data[this->ind(i, j)] = tmp[i * (tmpSize[1]) + j];
            }
            else if (j >= colSub)
            {
                this->m_data[this->ind(i, j)] = tmp[i * (tmpSize[1]) + j + numCols];
            }
        }
    }
    delete[] tmp;
}

/* returns resized matrix */
Matrix Matrix::resize(const unsigned int& numRows, const unsigned int& numCols)
{
    assert(numRows * numCols == this->m_num_elems);
    Matrix result = *this;
    result.m_size[0] = numRows;
    result.m_size[1] = numCols;
    return result;
}

/* returns transpose of matrix */
Matrix Matrix::transpose()
{
    Matrix result({this->m_size[1], this->m_size[0]});
    for (unsigned int i = 0; i < this->m_size[0]; ++i)
    {
        for (unsigned int j = 0; j < this->m_size[1]; ++j)
        {
            result.m_data[result.ind(j, i)] = this->m_data[this->ind(i, j)];
        }
    }
    return result;
}

/* checks if matrix is empty */
bool Matrix::is_empty() 
{
    return ((this->m_size[0] == 0) && (this->m_size[1] == 0) && (this->m_num_elems == 0) && (this->m_data == nullptr));
}

/* checks if matrix is square */
bool Matrix::is_square() 
{
    return (this->m_size[0] == this->m_size[1]);
}

/* checks if matrix is symmetric */
bool Matrix::is_symmetric()
{
    if (this->m_size[0] != this->m_size[1])
        return false;
    
    for (unsigned int i = 0; i < this->m_size[0]; ++i)
    {
        for (unsigned int j = 0; j < this->m_size[0]; ++j)
        {
            if (this->m_data[this->ind(i, j)] != this->m_data[this->ind(j, i)])
            {
                return false;
            }
        }
    }
    return true;
}

/* ************************************************************************* 
Matrix types */

/* returns identity matrix with given size */
Matrix Matrix::identity(const unsigned int& size)
{
    Matrix ident({size, size}, 0.0);
    for (unsigned int i = 0; i < size; ++i)
    {
        ident.m_data[ident.ind(i, i)] = 1.0;
    }
    return ident;
}

/* creates diagonal matrix given col vector of values */
Matrix Matrix::diag(const Matrix& vals)
{
    assert(vals.size()[0] > 0);
    assert(vals.size()[1] == 1);
    Matrix result({vals.size()[0], vals.size()[0]}, 0.0);
    for (unsigned int i = 0; i < vals.size()[0]; ++i)
    {
        result.m_data[result.ind(i, i)] = vals.m_data[vals.ind(i, 0)];
    }
    return result;
}

/* returns empty matrix */
Matrix Matrix::empty()
{
    Matrix result;
    return result;
}

/* ************************************************************************* 
Gauss-Jordan elimination */

/* swap two rows */
void Matrix::swap_rows(const unsigned int& r0, const unsigned int& r1) 
{
    for (unsigned int i = 0; i < this->m_size[1]; ++i)
    {
        double tmp = this->m_data[this->ind(r0, i)];
        this->m_data[this->ind(r0, i)] = this->m_data[this->ind(r1, i)];
        this->m_data[this->ind(r1, i)] = tmp;
    }
}

/* swap two cols */
void Matrix::swap_cols(const unsigned int& c0, const unsigned int& c1) 
{
    for (unsigned int i = 0; i < this->m_size[0]; ++i)
    {
        double tmp = this->m_data[this->ind(i, c0)];
        this->m_data[this->ind(i, c0)] = this->m_data[this->ind(i, c1)];
        this->m_data[this->ind(i, c1)] = tmp;
    }
}

/* scale a row */
void Matrix::scale_row(const unsigned int& r, const double& a)
{
    for (unsigned int i = 0; i < this->m_size[1]; ++i)
    {
        if (this->m_data[this->ind(r, i)] != 0.0)
        {
            this->m_data[this->ind(r, i)] *= a;
        }
    }
}

/* row axpy */
void Matrix::axpy_row(const unsigned int& ry, const unsigned int& rx, const double& a)
{
    for (unsigned int i = 0; i < this->m_size[1]; ++i)
    {
        this->m_data[this->ind(ry, i)] += a * this->m_data[this->ind(rx, i)];
    }
}

/* get leading entry col for a given row index */
unsigned int Matrix::leading_entry_col(const unsigned int& r)
{
    for (unsigned int i = 0; i < this->m_size[1]; ++i)
    {
        if (this->m_data[this->ind(r, i)] != 0.0)
            return i;
    }
    return this->m_size[1];
}

/* gauss-jordan elimination to get row-reduced echelon form 
    returns tuple of reduced matrix, pivot indexes and type of solution available
    firstly moves all zero rows to the bottom of the matrix
    then locates pivots by iterating through cols and identifying the first row with a non-zero element then bringing this to the highest possible row
    below the previous pivot's row
    then making sure all elements in the target col below this element are zero using axpy operation
    then all elements above each pivot are eliminated by iterating across the columns in the reverse direction

    returns tuple containing the row-reduced echelon form matrix, the pivot location indexes and a code:
        0 = inconsistent system
        1 = unique solution available
        2 = infinite solutions
*/
std::tuple<Matrix, std::vector<std::vector<unsigned int>>, unsigned int> Matrix::row_reduced_echelon_form() 
{
    unsigned int solutionType = 1;
    Matrix mat(*this);

    // get indexes of zero rows
    std::vector<unsigned int> zeroRowArr = {};
    for (unsigned int i = 0; i < mat.m_size[0]; ++i)
    {
        bool zero = true;
        for (unsigned int j = 0; j < mat.m_size[1]; ++j)
        {
            if (mat.m_data[mat.ind(i, j)] != 0.0)
                zero = false;
        }
        if (zero == true)
            zeroRowArr.push_back(i);
    }

    // swap out each row index with bottom rows to make sure all zero rows are on the bottom
    unsigned int zeroTargetRow = mat.m_size[0] - 1;
    for (unsigned int i = 0; i < zeroRowArr.size(); ++i)
    {
        mat.swap_rows(zeroRowArr[i], zeroTargetRow);
        zeroTargetRow--;
    }

    // locate pivots
    // loop down each row for each col to locate pivots and place them at the top most 'target row', which increments each time a pivot is found
    // then loop down each row in the pivot's column and axpy each row to make sure there are zeros below each pivot in the pivot's column
    // then increment the target col to find the next pivot
    unsigned int targetCol = 0;
    unsigned int targetRow = 0;
    std::vector<std::vector<unsigned int>> pivots = {};

    while ((targetRow < mat.m_size[0]) && (targetCol < mat.m_size[1]))
    {
        bool found = false;
        for (unsigned int i = targetRow; i < mat.m_size[0]; ++i)
        {
            if (mat.m_data[mat.ind(i, targetCol)] != 0.0)
            {
                found = true;
                mat.swap_rows(targetRow, i); // once finding the row with the non zero term in the target column, swap this to be at the target row
                mat.scale_row(targetRow, 1.0 / mat.m_data[mat.ind(targetRow, targetCol)]); // then scale the target row s.t. the target element is 1
                pivots.push_back(std::vector<unsigned int>({targetRow, targetCol})); // this is a pivot so add it to the pivot vector
                if (targetCol >= mat.m_size[0]) // this condition might need to be corrected... assumes that any col beyond the row size (assuming matrix is square) shouldn't have a pivot for consistency
                {
                    solutionType = 0; // pivot located outside on LH square region of interest
                }
                break;
            }
        }

        // check if a pivot was found in this column
        if (found == true)
        {
            for (unsigned int i = targetRow + 1; i < mat.m_size[0]; ++i)
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
    for (int i = pivots.size() - 1; i >= 0; --i)
    {
        for (unsigned int j = 0; j < pivots[i][0]; ++j)
        {
            if (mat.m_data[mat.ind(j, pivots[i][1])] != 0.0)
            {
                mat.axpy_row(j, pivots[i][0], -mat.m_data[mat.ind(j, pivots[i][1])] / mat.m_data[mat.ind(pivots[i][0], pivots[i][1])]);
            }
        }
    }

    if ((solutionType != 0) && (pivots.size() < mat.m_size[0])) // check that sufficient pivots are found (= number of rows) otherwise there are infinite solutions 
    {
        solutionType = 2;
    }

    return std::tuple<Matrix, std::vector<std::vector<unsigned int>>, unsigned int>(mat, pivots, solutionType);
}

/* solve system of equations 
    requires that A is square
*/
Matrix Matrix::solve_GJ(const Matrix& A, const Matrix& b) 
{
    Matrix augMat(A);
    augMat.insert_cols(A.m_size[1], b);

    std::tuple<Matrix, std::vector<std::vector<unsigned int>>, unsigned int> tup = augMat.row_reduced_echelon_form();
    Matrix reducedAugMat = std::get<0>(tup);
    unsigned int solutionType = std::get<2>(tup);

    if (solutionType == 1)
        // if there is a single solution then return this
        return reducedAugMat.get_region({0, reducedAugMat.m_size[1] - 1}, {reducedAugMat.m_size[0], 1});
    else
        // otherwise return empty matrix
        return Matrix::empty();
}

/* inverts a matrix by augmenting on right with identity and reducing to row echelon form thereby obtaining identity on the left and inverse on the right */
Matrix Matrix::inverse_GJ()
{ 
    if (this->size()[0] != this->size()[1])
        return Matrix::empty(); // if not a square matrix then return empty matrix as obviously not invertible
    
    Matrix augMat(*this);
    augMat.insert_cols(this->size()[1], Matrix::identity(this->size()[0]));

    std::tuple<Matrix, std::vector<std::vector<unsigned int>>, unsigned int> tup = augMat.row_reduced_echelon_form();
    Matrix reducedAugMat = std::get<0>(tup);
    unsigned int solutionType = std::get<2>(tup);

    if (solutionType == 1) // correct number of pivots found, all are located in the LH square region of interest
    {
        return reducedAugMat.get_region({0, this->size()[0]}, {this->size()[0], this->size()[0]});
    }
    else
        return Matrix::empty();
}

/* finds the null space of a matrix and returns its basis vectors 
    - obtains the row reduced form
    - then rearranges columns (and thus variables) to get it in the form
    { { I_[r * r] , C_[r * (k-r)] }, { 0_[(n-r) * r], 0_[(n-r) * (k-r)] } }
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
Matrix Matrix::nullBasis()
{
    // get the row echelon form of this matrix
    std::tuple<Matrix, std::vector<std::vector<unsigned int>>, unsigned int> tup = this->row_reduced_echelon_form();
    Matrix mat = std::get<0>(tup);
    std::vector<std::vector<unsigned int>> pivots = std::get<1>(tup);
    unsigned int solutionType = std::get<2>(tup);
    unsigned int rank = pivots.size();

    if (rank == this->m_size[1])
        return Matrix::empty(); // there are same number of pivots as variables therefore nullspace is { 0.0 }

    // create a mapping for variable rearrangement to ensure pivots are contained in a identity matrix on top LH corner of matrix
    std::vector<unsigned int> rearrangeMap = {};
    for (unsigned int i = 0; i < this->m_size[1]; ++i)
        rearrangeMap.push_back(i);
    
    // go through pivot positions to determine if swaps need to be made
    unsigned int targetCol = 0;
    for (unsigned int i = 0; i < rank; ++i)
    {
        if (pivots[i][1] != targetCol)
        {
            // need to swap columns, record this in the rearragement mapping
            mat.swap_cols(pivots[i][1], targetCol);
            rearrangeMap[pivots[i][1]] = targetCol;
            rearrangeMap[targetCol] = pivots[i][1];
        }
        targetCol++;
    }
    
    // construct null basis matrix
    mathlib::Matrix nullBasis({this->m_size[1], this->m_size[1] - rank}); // create matrix for null basis vectors
    nullBasis.set_region({0, 0}, -1.0 * mat.get_region({0, rank}, {rank, this->m_size[1] - rank})); // adds the -C arbitrary matrix
    nullBasis.set_region({rank, 0}, Matrix::identity(this->m_size[1] - rank)); // adds the identity 

    // rearrange rows to get original ordering of variables
    for (unsigned int i = 0; i < this->m_size[1]; ++i)
    {
        if (rearrangeMap[i] != i)
        {
            nullBasis.swap_rows(rearrangeMap[i], i);
            rearrangeMap[rearrangeMap[i]] = rearrangeMap[i];
            rearrangeMap[i] = i;
        }
    }

    return nullBasis;
}

/* ************************************************************************* 
Determinants */

/* returns the minor of the current matrix */
Matrix Matrix::minor(const unsigned int& r, const unsigned int& c) 
{
    Matrix result(*this);
    result.remove_rows(r, 1);
    result.remove_cols(c, 1);
    return result;
}

/* calculates the determinant using recursive formula - note this is not computationally efficient, just a proof of concept */
double Matrix::determinant()
{
    if ((this->m_size[0] == 1) && (this->m_size[1] == 1))
    {
        return this->m_data[0];
    }
    else
    {
        double det = 0.0;
        double sign = 1.0;
        for (unsigned int i = 0; i < this->m_size[1]; ++i)
        {
            Matrix min = this->minor(0, i);
            det += sign * this->m_data[this->ind(0, i)] * min.determinant();
            sign *= -1.0;  
        }
        return det;
    }
}

} // namespace mathlib