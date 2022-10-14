/*
Matrix object class
William Denny, 3rd September 2022
*/

#pragma once

#include <assert.h>
#include <vector>

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
    Matrix(int numRows, int numCols, const std::vector<std::vector<T>>& data);
    Matrix(int dir, const std::vector<T>& data); // row or col vector constructor
    Matrix(const Matrix<T>& matrix); // copy ctor
    ~Matrix();
    
    // getters/setters
    T get(int row, int col) const;
    Matrix<T> getRegion(int row0, int row1, int col0, int col1) const;
    Matrix<T> getRow(int row) const;
    Matrix<T> getCol(int col) const;
    T scalar() const;
    
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
    
    void insertRows(int row, int numRows, const Matrix<T>& matrix);
    void insertCols(int col, int numCols, const Matrix<T>& matrix);
    void insertRow(int row, const Matrix<T>& matrix);
    void insertCol(int col, const Matrix<T>& matrix);
    void appendRow(const Matrix<T>& matrix);
    void appendCol(const Matrix<T>& matrix);
    void removeRows(int row0, int row1);
    void removeCols(int col0, int col1);
    void removeRow(int row);
    void removeCol(int col);
    
    void addRow(int row, const Matrix<T>& matrix);
    void addCol(int col, const Matrix<T>& matrix);
    void subtractRow(int row, const Matrix<T>& matrix);
    void subtractCol(int col, const Matrix<T>& matrix);
    void scaleRow(int row, const T& coef);
    void scaleCol(int col, const T& coef);
    
    // row reduction (scale row already declared)
    void replaceRow(int row0, int row1, const T& coef);
    void swapRow(int row0, int row1);
    int getLeadingIndexFromRow(int row) const;
    int getLeadingIndexFromCol(int col) const;
    
    void transpose();
    Matrix<T> getTranspose() const;
    void echelonForm();
    void reducedEchelonForm();
    Matrix<T> getEchelonForm() const;
    Matrix<T> getReducedEchelonForm() const;
    int getNumPivotCols() const;
    bool isConsistent() const;
    void invert();
    Matrix<T> getInverse() const;
    
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
    assert(numRows > 0);
    assert(numCols > 0);
    
    m_numRows = numRows;
    m_numCols = numCols;
    m_data = new T[numRows * numCols];
}

template <class T>
Matrix<T>::Matrix(int numRows, int numCols, const T&  value)
{
    assert(numRows > 0);
    assert(numCols > 0);
    
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
    assert(numRows > 0);
    assert(numCols > 0);
    
    m_numRows = numRows;
    m_numCols = numCols;
    m_data = new T[numRows * numCols];
    for (int i = 0; i < numRows * numCols; ++i)
    {
        m_data[i] = data[i];        
    }
}

// constructor that takes data in std::vector form
template <class T>
Matrix<T>::Matrix(int numRows, int numCols, const std::vector<std::vector<T>>& data)
{
    assert(numRows > 0);
    assert(numCols > 0);
    assert(data.size() == numRows);
    
    m_numRows = numRows;
    m_numCols = numCols;
    m_data = new T[numRows * numCols];
    
    for (int i = 0; i < numRows; ++i)
    {
        assert(data[i].size() == numCols);
        
        for (int j = 0; j < numCols; ++j)
        {
            m_data[index(i, j)] = data[i][j];
        }        
    }
}

// constructor that takes data in std::vector form
template <class T>
Matrix<T>::Matrix(int dir, const std::vector<T>& data)
{
    assert(data.size() > 0);
    
    m_numRows = data.size();
    m_numCols = 1;
    m_data = new T[m_numRows * m_numCols];
    for (int i = 0; i < m_numRows; ++i)
    {
        m_data[index(i, 0)] = data[i];        
    }
    
    if (dir == 0)
    {
        transpose();
    }
}

// copy constructor
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

// destructor
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
Matrix<T> Matrix<T>::getRegion(int row0, int row1, int col0, int col1) const
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
Matrix<T> Matrix<T>::getRow(int row) const
{
    return getRegion(row, row, 0, m_numCols - 1);
}

// get column
template <class T>
Matrix<T> Matrix<T>::getCol(int col) const
{
    return getRegion(0, m_numRows - 1, col, col);
}

// if the matrix is 1x1 then allow scalar term to be returned
template <class T>
T Matrix<T>::scalar() const
{
    assert(m_numRows == 1);
    assert(m_numCols == 1);
    
    return get(0, 0);
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

template <class T>
void Matrix<T>::insertRows(int row, int numRows, const Matrix<T>& matrix)
{
    assert((m_numRows > 0) && (m_numCols > 0) && (m_data != nullptr)); // check this matrix isnt null
    assert((row >= 0) && (row <= m_numRows));
    assert(numRows > 0);
    assert(matrix.numRows() > 0);
    assert(matrix.numCols() == m_numCols);
    
    // copy matrix into tmp matrix to store data
    Matrix<T> tmp = *this;
    
    m_numRows += numRows;
    
    delete[] m_data;
    m_data = new T[m_numRows * m_numCols];
    
    for (int i = 0; i < m_numRows; ++i)
    {
        if (i < row)
        {
            setRow(i, tmp.getRow(i));
        }
        else if ((i >= row) && (i < row + numRows))
        {
            setRow(i, matrix.getRow(i - row));
        }
        else
        {
            setRow(i, tmp.getRow(i - numRows));
        }
    }
}

template <class T>
void Matrix<T>::insertCols(int col, int numCols, const Matrix<T>& matrix)
{
    assert((m_numRows > 0) && (m_numCols > 0) && (m_data != nullptr)); // check this matrix isn't null
    assert((col >= 0) && (col <= m_numCols));
    assert(numCols > 0);
    assert(matrix.numCols() > 0);
    assert(matrix.numRows() == m_numRows);
    
    // copy matrix into tmp matrix to store data
    Matrix<T> tmp = *this;
    
    m_numCols += numCols;
    
    delete[] m_data;
    m_data = new T[m_numRows * m_numCols];
    
    for (int i = 0; i < m_numCols; ++i)
    {
        if (i < col)
        {
            setCol(i, tmp.getCol(i));
        }
        else if ((i >= col) && (i < col + numCols))
        {
            setCol(i, matrix.getCol(i - col));
        }
        else
        {
            setCol(i, tmp.getCol(i - numCols));
        }
    }
}

template <class T>
void Matrix<T>::removeRows(int row, int numRows)
{
    assert((m_numRows > 0) && (m_numCols > 0) && (m_data != nullptr)); // check this matrix isn't null
    assert((row >= 0) && (row < m_numRows));
    assert(numRows > 0);
    
    // copy matrix into tmp matrix to store data
    Matrix<T> tmp = *this;
    
    m_numRows -= numRows;
    
    delete[] m_data;
    m_data = new T[m_numRows * m_numCols];
    
    for (int i = 0; i < m_numRows; ++i)
    {
        if (i < row)
        {
            setRow(i, tmp.getRow(i));
        }
        else
        {
            setRow(i, tmp.getRow(i + numRows));
        }
    }
}

template <class T>
void Matrix<T>::removeCols(int col, int numCols)
{
    assert((m_numRows > 0) && (m_numCols > 0) && (m_data != nullptr)); // check this matrix isn't null
    assert((col >= 0) && (col < m_numCols));
    assert(numCols > 0);
    
    // copy matrix into tmp matrix to store data
    Matrix<T> tmp = *this;
    
    m_numCols -= numCols;
    
    delete[] m_data;
    m_data = new T[m_numRows * m_numCols];
    
    for (int i = 0; i < m_numCols; ++i)
    {
        if (i < col)
        {
            setCol(i, tmp.getCol(i));
        }
        else
        {
            setCol(i, tmp.getCol(i + numCols));
        }
    }
}

template <class T>
void Matrix<T>::insertRow(int row, const Matrix<T>& matrix)
{
    insertRows(row, 1, matrix);
}

template <class T>
void Matrix<T>::insertCol(int col, const Matrix<T>& matrix)
{
    insertCols(col, 1, matrix);
}

template <class T>
void Matrix<T>::appendRow(const Matrix<T>& matrix)
{
    insertRow(m_numRows, matrix);
}

template <class T>
void Matrix<T>::appendCol(const Matrix<T>& matrix)
{
    insertCol(m_numCols, matrix);
}

template <class T>
void Matrix<T>::removeRow(int row)
{
    removeRows(row, 1); 
}

template <class T>
void Matrix<T>::removeCol(int col)
{
    removeCols(col, 1);
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

template <class T>
int Matrix<T>::getLeadingIndexFromRow(int row) const
{
    for (int i = 0; i < m_numCols; ++i)
    {
        if (get(row, i) != 0.0)
        {
            return i;
        }   
    }
    return -1;
}

template <class T>
int Matrix<T>::getLeadingIndexFromCol(int col) const
{
    for (int i = 0; i < m_numRows; ++i)
    {
        if (getLeadingIndexFromRow(i) == col)
        {
            return i;
        }
    }
    return -1;
}

// transpose matrix
template <class T>
void Matrix<T>::transpose()
{
    // create copy of matrix
    Matrix<T> tmp = *this;
    
    int tmp2 = m_numCols;
    m_numCols = m_numRows;
    m_numRows = tmp2;
    
    delete[] m_data;
    m_data = new T[m_numRows * m_numCols];
    
    for (int i = 0; i < m_numRows; ++i)
    {
        Matrix<T> tmp3 = tmp.getCol(i);
        tmp3.resize(1, m_numCols);
        setRow(i, tmp3);
    }
}

// return transpose matrix (note this is highly inefficient as it involves several copies of matrix, but simple implementation)
template <class T>
Matrix<T> Matrix<T>::getTranspose() const
{
    Matrix<T> tmp = *this;
    tmp.transpose();
    return tmp;
}

// change matrix into echelon form
// assumes matrix is in augmented form
template <class T>
void Matrix<T>::echelonForm()
{
    // not all pivots will on diagonal - therefore keep running var of row index as col changes
    int pivotRow = 0;
    
    // iterate through each col
    for (int i = 0; i < m_numCols; ++i)
    {        
        int skipRow = 0;
        
        // if this pivot row is zero, check if all terms below are zero too
        if (get(pivotRow, i) == 0.0)
        {
            int replaceIndex = -1;
            
            for (int j = pivotRow; j < m_numRows; ++j)
            {
                if ((get(j, i) != 0.0) && (replaceIndex == -1))
                {
                    replaceIndex = j;
                }
            }
            
            if (replaceIndex == -1)
            {
                // all terms on this col below row val are zero due to cancellation - skip this row as it's already sorted
                skipRow = 1;
            }
            else
            {
                replaceRow(pivotRow, replaceIndex, 1.0);
            }
        }
        
        // if a pivot is going to be made on this row then row reduce all rows below to get zeros below this pivot
        if (skipRow == 0)
        {
            // eliminate all zeros below
            for (int j = pivotRow + 1; j < m_numRows; ++j)
            {
                if (get(j, i) != 0.0)
                {
                    replaceRow(j, pivotRow, -1.0 * get(j, i) / get(pivotRow, i));
                }
            }
            
            pivotRow++;
        }
    }
}

// assumes matrix is in augmented form
template <class T>
void Matrix<T>::reducedEchelonForm()
{
    // put matrix in echelon form
    echelonForm();
    
    // assumes matrix is augmented matrix - loop through columns to make sure terms above are zero
    for (int i = m_numCols - 2; i > 0; --i)
    {
        // check if there is pivot on this row
        int pivotRow = getLeadingIndexFromCol(i);
        
        if (pivotRow != -1)
        {
            for (int j = pivotRow - 1; j >= 0; --j)
            {
                replaceRow(j, pivotRow, -1.0 * get(j, i) / get(pivotRow, i));
            }
            
            // normalise pivot
            scaleRow(pivotRow, 1.0 / get(pivotRow, i));
        }
    }
}

// assumes matrix is in augmented form
template <class T>
Matrix<T> Matrix<T>::getEchelonForm() const
{
    Matrix<T> tmp = *this;
    tmp.echelonForm();
    return tmp; 
}

// assumes matrix is in augmented form
template <class T>
Matrix<T> Matrix<T>::getReducedEchelonForm() const
{
    Matrix<T> tmp = *this;
    tmp.reducedEchelonForm();
    return tmp;
}

// assumes matrix is in augmented form
template <class T>
int Matrix<T>::getNumPivotCols() const
{
    Matrix<T> tmp = getReducedEchelonForm();
    int numPivots = 0;
    
    for (int i = 0; i < m_numCols; ++i)
    {
        if (tmp.getLeadingIndexFromCol(i) != -1)
        {
            numPivots++;
        }
    }
    
    return numPivots;
}

// assumes matrix is in augmented form
template <class T>
bool Matrix<T>::isConsistent() const
{
    Matrix<T> tmp = getEchelonForm();
    
    for (int i = 0; i < m_numRows; ++i)
    {
        if ((tmp.getLeadingIndexFromRow(i) == -1) && (get(i, m_numCols - 1) != 0.0))
        {
            return false;
        }
    }
    
    return true;
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

