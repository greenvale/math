/*
Matrix object class
William Denny, 3rd September 2022
*/

#pragma once

#include <assert.h>

template <class T>
class Matrix
{
private:
    T* m_data;
    int m_numRows;
    int m_numCols;
private:
    int index(int row, int col) const;
public:
    Matrix();
    Matrix(int numRows, int numCols);
    Matrix(int numRows, int numCols, const T& value);
    Matrix(int numRows, int numCols, const T* data);
    Matrix(const Matrix<T>& matrix); // copy ctor
    ~Matrix();
    
    // getters/setters
    T get(int row, int col) const;
    Matrix<T> getRegion(int row0, int row1, int col0, int col1);
    Matrix<T> getRow(int row);
    Matrix<T> getCol(int col);
    
    void set(int row, int col, const T& value);
    void setRegion(int row0, int row1, int col0, int col1, const Matrix<T>& matrix);
    void setRow(int row, const Matrix<T>& matrix);
    void setCol(int col, const Matrix<T>& matrix);
    
    int numRows() const;
    int numCols() const;
    
    // useful functions
    void print() const;
    
    // operator overloading (inside class)
    bool operator== (const Matrix<T>& rhs) const;
    
    Matrix<T>& operator= (const Matrix<T>& rhs);
    void operator+= (const Matrix<T>& rhs);
    void operator-= (const Matrix<T>& rhs);
    void operator*= (const T& rhs);
    void operator/= (const T& rhs);
    
    // matrix functions
    void resize(int numRows, int numCols);
    
    // row/col operations
    void addRow(int row, const Matrix<T>& matrix);
    void addCol(int col, const Matrix<T>& matrix);
    void subtractRow(int row, const Matrix<T>& matrix);
    void subtractCol(int col, const Matrix<T>& matrix);
    void scaleRow(int row, const T& coef);
    void scaleCol(int col, const T& coef);
    
    // row reduction (scale row already declared)
    void replaceRow(int row0, int row1, const T& coef);
    void swapRow(int row0, int row1);
    
    Matrix<T> echelonForm();
    
};

// operator overloading (outside class)
template <class U> Matrix<U> operator+ (const Matrix<U>& lhs, const Matrix<U>& rhs);
template <class U> Matrix<U> operator+ (const U& lhs, const Matrix<U>& rhs);
template <class U> Matrix<U> operator+ (const Matrix<U>& lhs, const U& rhs);

template <class U> Matrix<U> operator- (const Matrix<U>& lhs, const Matrix<U>& rhs);
template <class U> Matrix<U> operator- (const U& lhs, const Matrix<U>& rhs);
template <class U> Matrix<U> operator- (const Matrix<U>& lhs, const U& rhs);

template <class U> Matrix<U> operator* (const Matrix<U>& lhs, const Matrix<U>& rhs);
template <class U> Matrix<U> operator* (const U& lhs, const Matrix<U>& rhs);
template <class U> Matrix<U> operator* (const Matrix<U>& lhs, const U& rhs);

template <class U> Matrix<U> operator/ (const Matrix<U>& lhs, const U& rhs);

/* 
======================================================================
CONSTRUCTOR / DESTRUCTOR
*/ 

template <class T>
Matrix<T>::Matrix()
{
    m_numRows = 0;
    m_numCols = 0;
    m_data = nullptr;
}

template <class T>
Matrix<T>::Matrix(int numRows, int numCols)
{
    m_numRows = numRows;
    m_numCols = numCols;
    m_data = new T[numRows * numCols];
}

template <class T>
Matrix<T>::Matrix(int numRows, int numCols, const T&  value)
{
    m_numRows = numRows;
    m_numCols = numCols;
    m_data = new T[numRows * numCols];
    for (int i = 0; i < numRows * numCols; ++i)
    {
        m_data[i] = value;
    }
}

template <class T>
Matrix<T>::Matrix(int numRows, int numCols, const T* data)
{
    m_numRows = numRows;
    m_numCols = numCols;
    m_data = new T[numRows * numCols];
    for (int i = 0; i < numRows * numCols; ++i)
    {
        m_data[i] = data[i];        
    }
}

template <class T>
Matrix<T>::Matrix(const Matrix<T>& matrix)
{
    //std::cout << "Copy constructor" << std::endl;
    m_numRows = matrix.numRows();
    m_numCols = matrix.numCols();
    m_data = new T[matrix.numRows() * matrix.numCols()];
    for (int i = 0; i < matrix.numRows(); ++i)
    {
        for (int j = 0; j < matrix.numCols(); ++j)
        {
            m_data[index(i, j)] = matrix.get(i, j);
        }
    }
}

template <class T>
Matrix<T>::~Matrix()
{
    //std::cout << "Destructor" << std::endl;
    if (m_data != nullptr) 
    {
        delete[] m_data;
    }
}

/* 
======================================================================
GETTER / SETTER
*/ 

// get element
template <class T>
T Matrix<T>::get(int row, int col) const
{
    assert(row <= m_numRows);
    assert(col <= m_numCols);
    
    return m_data[index(row, col)];
}

// get region
template <class T>
Matrix<T> Matrix<T>::getRegion(int row0, int row1, int col0, int col1)
{
    assert((row0 >= 0) && (row0 < m_numRows));
    assert((row1 >= 0) && (row1 < m_numRows));
    assert((row1 >= row0));
    assert((col0 >= 0) && (col0 < m_numCols));
    assert((col1 >= 0) && (col1 < m_numCols));
    assert((col1 >= col0));
    
    int numRows = row1 - row0 + 1;
    int numCols = col1 - col0 + 1;
    Matrix<T> result(numRows, numCols);
    
    for (int i = 0; i < numRows; ++i)
    {
        for (int j = 0; j < numCols; ++j)
        {
            result.set(i, j, get(i + row0, j + col0));
        }
    }
    
    return result;
}

// get row
template <class T>
Matrix<T> Matrix<T>::getRow(int row)
{
    return getRegion(row, row, 0, m_numCols - 1);
}

// get column
template <class T>
Matrix<T> Matrix<T>::getCol(int col)
{
    return getRegion(0, m_numRows - 1, col, col);
}


// set element
template <class T>
void Matrix<T>::set(int row, int col, const T& value)
{
    assert(row <= m_numRows);
    assert(col <= m_numCols);
    
    m_data[index(row, col)] = value;
}

// set region
template <class T>
void Matrix<T>::setRegion(int row0, int row1, int col0, int col1, const Matrix<T>& matrix)
{
    assert((row0 >= 0) && (row0 < m_numRows));
    assert((row1 >= 0) && (row1 < m_numRows));
    assert((row1 >= row0));
    assert((col0 >= 0) && (col0 < m_numCols));
    assert((col1 >= 0) && (col1 < m_numCols));
    assert((col1 >= col0));
    
    int numRows = row1 - row0 + 1;
    int numCols = col1 - col0 + 1;
    
    assert(numRows == matrix.numRows());
    assert(numCols == matrix.numCols());
    
    for (int i = 0; i < numRows; ++i)
    {
        for (int j = 0; j < numCols; ++j)
        {
            set(i + row0, j + col0, matrix.get(i, j));
        }
    }
}

// set row
template <class T>
void Matrix<T>::setRow(int row, const Matrix<T>& matrix)
{
    setRegion(row, row, 0, m_numCols - 1, matrix);
}

// set column
template <class T>
void Matrix<T>::setCol(int col, const Matrix<T>& matrix)
{
    setRegion(0, m_numRows - 1, col, col, matrix);
}

// get num rows
template <class T>
int Matrix<T>::numRows() const
{
    return m_numRows;
}

// get num cols
template <class T>
int Matrix<T>::numCols() const
{
    return m_numCols;
}

/* 
======================================================================
OPERATOR OVERLOADING
*/ 

// compare operator
template <class T>
bool Matrix<T>::operator== (const Matrix<T>& rhs) const
{
    assert(m_numRows == rhs.numRows());
    assert(m_numCols == rhs.numCols());
    
    for (int i = 0; i < m_numRows; ++i)
    {
        for (int j = 0; j < m_numCols; ++j)
        {
            if (m_data[index(i, j)] != rhs.get(i,j))
            {
                return false;
            }
        }
    }
    
    return true;
}

// assignment operator
template <class T>
Matrix<T>& Matrix<T>::operator= (const Matrix<T>& rhs)
{
    //std::cout << "Assignment" << std::endl;
    
    if (m_data != nullptr)
    {
        delete[] m_data;
    }
    
    m_numRows = rhs.numRows();
    m_numCols = rhs.numCols();
    
    if ((m_numRows > 0) && (m_numCols > 0))
    {
        m_data = new T[m_numRows * m_numCols];
        
        for (int i = 0; i < m_numRows; ++i)
        {
            for (int j = 0; j < m_numCols; ++j)
            {
                m_data[index(i, j)] = rhs.get(i, j);
            }
        }
    }
    
    return *(this);
}

// += operator
template <class T>
void Matrix<T>::operator+= (const Matrix<T>& rhs)
{
    assert(m_numRows == rhs.numRows());
    assert(m_numCols == rhs.numCols());

    for (int i = 0; i < m_numRows; ++i)
    {
        for (int j = 0; j < m_numCols; ++j)
        {
            m_data[index(i, j)] += rhs.get(i, j);
        }
    }
}

// -= operator
template <class T>
void Matrix<T>::operator-= (const Matrix<T>& rhs)
{
    assert(m_numRows == rhs.numRows());
    assert(m_numCols == rhs.numCols());
    
    for (int i = 0; i < m_numRows; ++i)
    {
        for (int j = 0; j < m_numCols; ++j)
        {
            m_data[index(i, j)] -= rhs.get(i, j);
        }
    }
}

// *= operator
template <class T>
void Matrix<T>::operator*= (const T& rhs)
{
    assert(m_numRows == rhs.numRows());
    assert(m_numCols == rhs.numCols());
    
    for (int i = 0; i < m_numRows; ++i)
    {
        for (int j = 0; j < m_numCols; ++j)
        {
            m_data[index(i, j)] *= rhs;
        }
    }
}

// /= operator
template <class T>
void Matrix<T>::operator/= (const T& rhs)
{
    assert(m_numRows == rhs.numRows());
    assert(m_numCols == rhs.numCols());
    
    for (int i = 0; i < m_numRows; ++i)
    {
        for (int j = 0; j < m_numCols; ++j)
        {
            m_data[index(i, j)] /= rhs;
        }
    }
}

// + operator
// matrix + matrix
template <class T>
Matrix<T> operator+ (const Matrix<T>& lhs, const Matrix<T>& rhs)
{
    assert(lhs.numRows() == rhs.numRows());
    assert(lhs.numCols() == rhs.numCols());
    
    int numRows = lhs.numRows();
    int numCols = lhs.numCols();
    Matrix<T> result(numRows, numCols);
    
    for (int i = 0; i < numRows; ++i)
    {
        for (int j = 0; j < numCols; ++j)
        {
            result.set(i, j, lhs.get(i, j) + rhs.get(i, j));
        }
    }

    return result;    
}

// scalar + matrix
template <class T>
Matrix<T> operator+ (const T& lhs, const Matrix<T>& rhs)
{
    int numRows = rhs.numRows();
    int numCols = rhs.numCols();
    Matrix<T> result(numRows, numCols);
    
    for (int i = 0; i < numRows; ++i)
    {
        for (int j = 0; j < numCols; ++j)
        {
            result.set(i, j, lhs + rhs.get(i, j));
        }
    }
    
    return result;
}

// matrix + scalar
template <class T>
Matrix<T> operator+ (const Matrix<T>& lhs, const T& rhs)
{
    int numRows = lhs.numRows();
    int numCols = lhs.numCols();
    Matrix<T> result(numRows, numCols);
    
    for (int i = 0; i < numRows; ++i)
    {
        for (int j = 0; j < numCols; ++j)
        {
            result.set(i, j, lhs.get(i, j) + rhs);
        }
    }
    
    return result;
}

// - operator
// matrix - matrix
template <class T>
Matrix<T> operator- (const Matrix<T>& lhs, const Matrix<T>& rhs)
{
    assert(lhs.numRows() == rhs.numRows());
    assert(lhs.numCols() == rhs.numCols());
    
    int numRows = lhs.numRows();
    int numCols = lhs.numCols();
    Matrix<T> result(numRows, numCols);
    
    for (int i = 0; i < numRows; ++i)
    {
        for (int j = 0; j < numCols; ++j)
        {
            result.set(i, j, lhs.get(i, j) - rhs.get(i, j));
        }
    }

    return result;    
}

// scalar - matrix
template <class T>
Matrix<T> operator- (const T& lhs, const Matrix<T>& rhs)
{
    int numRows = rhs.numRows();
    int numCols = rhs.numCols();
    Matrix<T> result(numRows, numCols);
    
    for (int i = 0; i < numRows; ++i)
    {
        for (int j = 0; j < numCols; ++j)
        {
            result.set(i, j, lhs - rhs.get(i, j));
        }
    }
    
    return result;
}

// matrix - scalar
template <class T>
Matrix<T> operator- (const Matrix<T>& lhs, const T& rhs)
{
    int numRows = lhs.numRows();
    int numCols = lhs.numCols();
    Matrix<T> result(numRows, numCols);
    
    for (int i = 0; i < numRows; ++i)
    {
        for (int j = 0; j < numCols; ++j)
        {
            result.set(i, j, lhs.get(i, j) - rhs);
        }
    }
    
    return result;
}

// * operator
// matrix * matrix
template <class T>
Matrix<T> operator* (const Matrix<T>& lhs, const Matrix<T>& rhs)
{
    assert(lhs.numCols() == rhs.numRows());
    
    int numRows = lhs.numRows();
    int numCols = rhs.numCols();
    Matrix<T> result(numRows, numCols);
    
    for (int i = 0; i < numRows; ++i)
    {
        for (int j = 0; j < numCols; ++j)
        {   
            T entry = 0.0;
            for (int k = 0; k < lhs.numCols(); ++k)
            {
                entry += lhs.get(i, k) * rhs.get(k, j);
            }
            result.set(i, j, entry);
        }
    }
    
    return result;
}

// scalar * matrix
template <class T>
Matrix<T> operator* (const T& lhs, const Matrix<T>& rhs)
{
    int numRows = rhs.numRows();
    int numCols = rhs.numCols();
    Matrix<T> result(numRows, numCols);
    
    for (int i = 0; i < numRows; ++i)
    {
        for (int j = 0; j < numCols; ++j)
        {
            result.set(i, j, lhs * rhs.get(i, j));
        }
    }
    
    return result;
}

// matrix * scalar
template <class T>
Matrix<T> operator* (const Matrix<T>& lhs, const T& rhs)
{
    int numRows = lhs.numRows();
    int numCols = lhs.numCols();
    Matrix<T> result(numRows, numCols);
    
    for (int i = 0; i < numRows; ++i)
    {
        for (int j = 0; j < numCols; ++j)
        {
            result.set(i, j, lhs.get(i, j) * rhs);
        }
    }
    
    return result;
}

// matrix / scalar
template <class T>
Matrix<T> operator/ (const Matrix<T>& lhs, const T& rhs)
{
    int numRows = lhs.numRows();
    int numCols = lhs.numCols();
    Matrix<T> result(numRows, numCols);
    
    for (int i = 0; i < numRows; ++i)
    {
        for (int j = 0; j < numCols; ++j)
        {
            result.set(i, j, lhs.get(i, j) / rhs);
        }
    }
    
    return result;
}

/* 
======================================================================
USEFUL FUNCTIONS
*/ 

template <class T>
void Matrix<T>::print() const
{
    for (int i = 0; i < m_numRows; ++i)
    {
        for (int j = 0; j < m_numCols; ++j) 
        {
            std::cout << get(i, j) << "\t";
        }
        std::cout << std::endl;
    }
}

/* 
======================================================================
MATRIX FUNCTIONS
*/ 

// resize matrix to new shape
template <class T>
void Matrix<T>::resize(int numRows, int numCols)
{
    assert(numRows * numCols == m_numRows * m_numCols);
    
    m_numRows = numRows;
    m_numCols = numCols;
}

// Row/col operations

template <class T>
void Matrix<T>::addRow(int row, const Matrix<T>& matrix)
{
    assert((row >= 0) && (row < m_numRows));
    assert(matrix.numRows() == 1);
    assert(matrix.numCols() == m_numCols);
    
    setRow(row, getRow(row) + matrix);
}

template <class T>
void Matrix<T>::addCol(int col, const Matrix<T>& matrix)
{
    assert((col >= 0) && (col < m_numCols));
    assert(matrix.numCols() == 1);
    assert(matrix.numRows() == m_numRows);
    
    setCol(col, getCol(col) + matrix);
}

template <class T>
void Matrix<T>::subtractRow(int row, const Matrix<T>& matrix)
{
    assert((row >= 0) && (row < m_numRows));
    assert(matrix.numRows() == 1);
    assert(matrix.numCols() == m_numCols);
    
    setRow(row, getRow(row) - matrix);
}

template <class T>
void Matrix<T>::subtractCol(int col, const Matrix<T>& matrix)
{
    assert((col >= 0) && (col < m_numCols));
    assert(matrix.numCols() == 1);
    assert(matrix.numRows() == m_numRows);
    
    setCol(col, getCol(col) - matrix);
}

template <class T>
void Matrix<T>::scaleRow(int row, const T& coef)
{
    assert((row >= 0) && (row < m_numRows));
    
    setRow(row, getRow(row) * coef);
}

template <class T>
void Matrix<T>::scaleCol(int col, const T& coef)
{
    assert((col >= 0) && (col < m_numCols));
    
    setCol(col, getCol(col) * coef);
}

// Row reduction operations

template <class T>
void Matrix<T>::replaceRow(int row0, int row1, const T& coef)
{
    addRow(row0, coef * getRow(row1));
}

template <class T>
void Matrix<T>::swapRow(int row0, int row1)
{
    Matrix<T> temp = getRow(row0);
    setRow(row0, getRow(row1));
    setRow(row1, temp);
}

// return matrix in echelon form
template <class T>
Matrix<T> Matrix<T>::echelonForm()
{
    Matrix<T> result = *this;
    
    // iterate through each col and ensure leading row at top
    for (int i = 0; i < m_numCols; ++i)
    {

    }
    
    return result;
}


/* 
======================================================================
PRIVATE FUNCTIONS
*/ 

// returns index for row-ordered matrix
template <class T>
int Matrix<T>::index(int row, int col) const
{    
    return row*m_numCols + col;
}

