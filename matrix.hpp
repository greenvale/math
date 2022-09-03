/*
Matrix object class
William Denny, 3rd September 2022
*/

#pragma once

#include <assert.h>

template <class T>
class Matrix2
{
private:
    T* m_data;
    int m_numRows;
    int m_numCols;
private:
    int index(int row, int col) const;
public:
    // ctors
    Matrix2();
    Matrix2(int numRows, int numCols);
    Matrix2(int numRows, int numCols, const T& value);
    Matrix2(int numRows, int numCols, const T* data);
    Matrix2(const Matrix2<T>& matrix); // copy ctor
    
    // dtor
    ~Matrix2();
    
    T get(int row, int col) const;
    void set(int row, int col, const T& value);
    int numRows() const;
    int numCols() const;
    
    void resize(int numRows, int numCols);
    
    void print() const;
    
    // operator overloading
    bool operator== (const Matrix2<T>& rhs);
    
    Matrix2<T>& operator= (const Matrix2<T>& rhs);
    void operator+= (const Matrix2<T>& rhs);
    void operator-= (const Matrix2<T>& rhs);
    void operator*= (const T& rhs);
    void operator/= (const T& rhs);
};

template <class U> Matrix2<U> operator+ (const Matrix2<U>& lhs, const Matrix2<U>& rhs);
template <class U> Matrix2<U> operator+ (const U& lhs, const Matrix2<U>& rhs);
template <class U> Matrix2<U> operator+ (const Matrix2<U>& lhs, const U& rhs);

template <class U> Matrix2<U> operator- (const Matrix2<U>& lhs, const Matrix2<U>& rhs);
template <class U> Matrix2<U> operator- (const U& lhs, const Matrix2<U>& rhs);
template <class U> Matrix2<U> operator- (const Matrix2<U>& lhs, const U& rhs);

template <class U> Matrix2<U> operator* (const Matrix2<U>& lhs, const Matrix2<U>& rhs);
template <class U> Matrix2<U> operator* (const U& lhs, const Matrix2<U>& rhs);
template <class U> Matrix2<U> operator* (const Matrix2<U>& lhs, const U& rhs);

template <class U> Matrix2<U> operator/ (const Matrix2<U>& lhs, const U& rhs);

/* 
======================================================================
CONSTRUCTOR / DESTRUCTOR
*/ 

template <class T>
Matrix2<T>::Matrix2()
{
    m_numRows = 0;
    m_numCols = 0;
    m_data = nullptr;
}

template <class T>
Matrix2<T>::Matrix2(int numRows, int numCols)
{
    m_numRows = numRows;
    m_numCols = numCols;
    m_data = new T[numRows * numCols];
}

template <class T>
Matrix2<T>::Matrix2(int numRows, int numCols, const T&  value)
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
Matrix2<T>::Matrix2(int numRows, int numCols, const T* data)
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
Matrix2<T>::Matrix2(const Matrix2<T>& matrix)
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
Matrix2<T>::~Matrix2()
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
T Matrix2<T>::get(int row, int col) const
{
    assert(row <= m_numRows);
    assert(col <= m_numCols);
    
    return m_data[index(row, col)];
}

// set element
template <class T>
void Matrix2<T>::set(int row, int col, const T& value)
{
    assert(row <= m_numRows);
    assert(col <= m_numCols);
    
    m_data[index(row, col)] = value;
}

// get num rows
template <class T>
int Matrix2<T>::numRows() const
{
    return m_numRows;
}

// get num cols
template <class T>
int Matrix2<T>::numCols() const
{
    return m_numCols;
}

// resize matrix to new shape
template <class T>
void Matrix2<T>::resize(int numRows, int numCols)
{
    assert(numRows * numCols == m_numRows * m_numCols);
    
    m_numRows = numRows;
    m_numCols = numCols;
}

/* 
======================================================================
OPERATOR OVERLOADING
*/ 

// compare operator
template <class T>
bool Matrix2<T>::operator== (const Matrix2<T>& rhs)
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
Matrix2<T>& Matrix2<T>::operator= (const Matrix2<T>& rhs)
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
void Matrix2<T>::operator+= (const Matrix2<T>& rhs)
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
void Matrix2<T>::operator-= (const Matrix2<T>& rhs)
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
void Matrix2<T>::operator*= (const T& rhs)
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
void Matrix2<T>::operator/= (const T& rhs)
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
Matrix2<T> operator+ (const Matrix2<T>& lhs, const Matrix2<T>& rhs)
{
    assert(lhs.numRows() == rhs.numRows());
    assert(lhs.numCols() == rhs.numCols());
    
    int numRows = lhs.numRows();
    int numCols = lhs.numCols();
    Matrix2<T> result(numRows, numCols);
    
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
Matrix2<T> operator+ (const T& lhs, const Matrix2<T>& rhs)
{
    int numRows = rhs.numRows();
    int numCols = rhs.numCols();
    Matrix2<T> result(numRows, numCols);
    
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
Matrix2<T> operator+ (const Matrix2<T>& lhs, const T& rhs)
{
    int numRows = lhs.numRows();
    int numCols = lhs.numCols();
    Matrix2<T> result(numRows, numCols);
    
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
Matrix2<T> operator- (const Matrix2<T>& lhs, const Matrix2<T>& rhs)
{
    assert(lhs.numRows() == rhs.numRows());
    assert(lhs.numCols() == rhs.numCols());
    
    int numRows = lhs.numRows();
    int numCols = lhs.numCols();
    Matrix2<T> result(numRows, numCols);
    
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
Matrix2<T> operator- (const T& lhs, const Matrix2<T>& rhs)
{
    int numRows = rhs.numRows();
    int numCols = rhs.numCols();
    Matrix2<T> result(numRows, numCols);
    
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
Matrix2<T> operator- (const Matrix2<T>& lhs, const T& rhs)
{
    int numRows = lhs.numRows();
    int numCols = lhs.numCols();
    Matrix2<T> result(numRows, numCols);
    
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
Matrix2<T> operator* (const Matrix2<T>& lhs, const Matrix2<T>& rhs)
{
    assert(lhs.numCols() == rhs.numRows());
    
    int numRows = lhs.numRows();
    int numCols = rhs.numCols();
    Matrix2<T> result(numRows, numCols);
    
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
Matrix2<T> operator* (const T& lhs, const Matrix2<T>& rhs)
{
    int numRows = rhs.numRows();
    int numCols = rhs.numCols();
    Matrix2<T> result(numRows, numCols);
    
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
Matrix2<T> operator* (const Matrix2<T>& lhs, const T& rhs)
{
    int numRows = lhs.numRows();
    int numCols = lhs.numCols();
    Matrix2<T> result(numRows, numCols);
    
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
Matrix2<T> operator/ (const Matrix2<T>& lhs, const T& rhs)
{
    int numRows = lhs.numRows();
    int numCols = lhs.numCols();
    Matrix2<T> result(numRows, numCols);
    
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
PRINT FUNCTIONS
*/ 

template <class T>
void Matrix2<T>::print() const
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
PRIVATE FUNCTIONS
*/ 

// returns index for row-ordered matrix
template <class T>
int Matrix2<T>::index(int row, int col) const
{    
    return row*m_numCols + col;
}

