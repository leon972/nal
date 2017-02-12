/*****************************************************************************

    Copyright (c) Leonardo Berti 2007-2013
    
    This file is part of NAL Numerical Analysis Library

    NAL is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    NAL is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with NAL.  If not, see <http://www.gnu.org/licenses/>.

 ******************************************************************************/
/*!
 * \file tmatrix.h
 * \brief Generic NxM (N rows M columns) template-base matrix class.
 */

#ifndef TMATRIX_H_INCLUDED
#define TMATRIX_H_INCLUDED

#include "math/mathdefs.h"
#include "libmath.h"
#include "tvector.h"
#include <stddef.h>
#include <iostream>

namespace mathengine {

    /**
     * unused constant:doxygen trick
     */
    const int zero = 0;

    /**
     * Template-based nxm matrix.
     * Implements all basic operators (addition,subtrction,multiplication)
     * and elements accessing methods.
     * Example:
     * ------------------------------------------------------------------
     * \code
     *  TMatrix<double> m(3, 3); //creates a 3x3 numeric matrix of type double *
     *
     *  m[0][0]=1; //sets some elements
     *  m[0][1]=2;
     *  m[0][2]=3;
     *
     *  cout<<"m="<<endl<<m<<endl; //display the matrix
     *
     * TMatrix<double> m1(3,3);
     *
     *  m1.setDiagonal(1); //sets all elements on the diagonal to the value 1
     *
     *  TMatrix<double> res(3,3);
     *
     *  res=m*m1; //matrix multiplication
     *  m+=res;   //matrix addition
     * 
     * \endcode
     * ------------------------------------------------------------------
     *
     */
    template <class T> class TMatrix {
    public:

        /**
         * Creates a marix with nrows rows and ncolumns columns
         * @param nrows number of rows
         * @param ncolumns number of columns
         */
        TMatrix(size_t nrows, size_t ncolumns) : data(0), nRows(nrows), nColumns(ncolumns), nsize(0) {
            resize(nrows, ncolumns, true);
        }

        /**
         * Creates a marix with nrows rows and ncolumns columns
         * @param nrows number of rows
         * @param ncolumns number of columns
         * @param clearMemory if set, reset the internal memory of the matrix
         */
        TMatrix(size_t nrows, size_t ncolumns, bool clearMemory) : data(0), nRows(nrows), nColumns(ncolumns), nsize(0) {
            resize(nrows, ncolumns, clearMemory);
        }

        virtual ~TMatrix() {
            if (data) {
                delete[] data;
                data = 0;
            }
        }

        /**
         * Assignment operator.
         * If m ha the same size of this matrix,elements of m are copied thos this matrix
         * otherwise,this matrix is reallocated to fit the size of m and elements are copied from m.
         * @param m
         * @return
         */
        TMatrix<T> & operator=(const TMatrix<T> &m) {

            if (m.nsize == nsize) {

                if (&m == this) return *this;

                memcpy(data, m.data, sizeof (T) * nsize);

            } else {
                resize(m.getRows(), m.getColumns(), false);
                memcpy(data, m.data, sizeof (T) * nsize);
            }

            return *this;
        }

        /**
         * Sets all elements to v
         * @param v value
         * @return this matrix
         */
        TMatrix<T> & operator=(const T &v) {

            for (size_t i = 0; i < nsize; i++) {
                data[i] = v;
            }
        }

        /**
         * Number of rows
         * @return the number of rows
         */
        size_t getRows() const {
            return nRows;
        }

        /**
         * Number of columns
         * @return the number of columns
         */
        size_t getColumns() const {
            return nColumns;
        }

        /**
         * Clear the matrix and reallocates the memory.
         * Use this method to change the matrix size
         * @param nrows new number of rows
         * @param ncolumns new number of columns
         * @param zero_element initial value of all elements.
         * @param resetAllelements if true,all elements are set to zero_element
         */
        inline void resize(int nrows, int ncolumns, bool resetAllElements = true) {

            if (data) {
                delete[] data;
                data = 0;
            }
            nsize = nrows*ncolumns;
            nRows = nrows, nColumns = ncolumns;
            data = new T[nsize];

            if (resetAllElements) {
                zeroMemory();
            }
        }

        /**
         * This operator is used to access the components of the matrix.
         * eg. m.[idxrow][idxcolumn].
         * \warning this method doesn't perform bounds check for optimization reasons.
         * @param row index of the row (0..n)
         * @return the pointer to the row.
         */
        inline T * operator [] (size_t row) const throw () {
            return &data[nColumns * row];
        }

        /**
         * Subscript operator, can be used to access the matrix elements.
         * This method can be used to modify elements in the matrix.
         * eq. T &p=matrix(i,j);
         * p->someMethod(...)
         * \warning this method doesn't perform bounds check for optimization reasons.
         * @param r row (first row has index 0)
         * @param c column (first column has index 0)
         * @return the reference of the element at (r,c)
         *
         */
        inline T & operator ()(size_t r, size_t c) throw () {
            return data[r * nColumns + c];
        }

        /**
         * Subscript operator (const version), can be used to access the matrix elements
         * \warning this method doesn't perform bounds check for optimization reasons.
         * @param r row
         * @param c column
         * @return the element at (r,c)
         *
         */
        inline T operator()(size_t r, size_t c) const throw () {
            return data[r * nColumns + c];
        }

        /**
         * Get the reference of an element with indices bound check.
         * @param row index of row (first row has index zero)
         * @param column index of column (first column has index zero)
         * @return the reference of the element at the specified row and column.
         * @throw MathException if row and/or column are out of bounds.
         */
        T& getElementAt(size_t row, size_t column) const throw (MathException) {

            if (row >= nRows) {
                MathException ex("Error accessing matrix element:row %d out of bounds.", row);
                throw ex;
            }

            if (column >= nColumns) {
                MathException ex("Error accessing matrix element:column %d out of bounds.", column);
                throw ex;
            }

            return data[row * nColumns + column];
        }

        /**
         * Set an element of the matrix with indices bound check.
         * @param row index of row (first row has index zero)
         * @param column index of column (first column has index zero)
         * @param value the new value of the element
         * @throw MathException if row and/or column are out of bounds.
         */
        void setElementAt(size_t row, size_t column, const T& value) throw (MathException) {

            if (row >= nRows) {
                MathException ex("Error setting matrix element:row %d out of bounds.", row);
                throw ex;
            }

            if (column >= nColumns) {
                MathException ex("Error setting matrix element:column %d out of bounds.", column);
                throw ex;
            }

            data[row * nColumns + column] = value;
        }

        /**
         * Checks if this matrix has the same rows and columns as m
         * @param m a matrix
         */
        inline bool hasSameSize(const TMatrix<T> &m) const throw () {
            return getRows() == m.getRows() && getColumns() == m.getColumns();
        }

        /**
         * Return true if this matrix is a square matrix (same number of rows and columns)
         * @return
         */
        inline bool isSquare() const throw () {
            return nColumns == nRows;
        }

        /**
         * Clear memory (all bytes set to zero)
         */
        inline void zeroMemory() throw () {
            memset(data, 0, sizeof (T) * nsize);
        }

        /**
         * Set all elements to value
         */
        inline void setAllElements(const T &value) throw () {
            for (size_t i = 0; i < nsize; i++) {
                data[i] = value;
            }
        }

        /**
         * Set the elements of the diagonal to the same value.
         */
        inline void setDiagonal(const T &value) throw (MathException) {

            if (nColumns != nRows) throw MathException("Cannot set the diagonal:matrix is not square.");

            for (size_t i = 0; i < nColumns; i++) {
                data[i * nColumns + i] = value;
            }
        }

        //--------------- addition --------------------------------------------

        /**
         * Add m to this matrix
         * @param m another matrix
         * @thorw MathException if m ha different size from this matrix
         */
        inline void add(const TMatrix<T> &matrix) throw (MathException) {

            if (nRows != matrix.nRows || nColumns != matrix.nColumns) throw MathException("Error adding matrix:invalid matrix size");

            const T* pdata = matrix.data;

            for (size_t i = 0; i < nsize; i++) {
                data[i] += pdata[i];
            }
        }

        /**
         * Fast addition.
         * Adds the specified matrix to this matrix and places the result in the result matrix.
         * \warning For optimization reasons this method doesn't perform arguments and bounds check.
         * @param m the matrix to add.This matrix must have the same size as this matrix.
         * @param result The matrix where the result of the addition is places.Points to an already allocated matrix of the same size as this matrix.
         */
        inline void add(TMatrix<T> *result, const TMatrix<T> &matrix) const throw () {

            T* dest_data = result->data;

            const T* src_data = matrix.data;

            for (size_t i = 0; i < nsize; i++) {

                dest_data[i] = data[i] + src_data[i];
            }
        }

        /**
         * Addition operator
         * 
         * @param m a matrix with the same number of rows and coloumns of this matrix
         * @throw MathException if matrices have different size.
         * @return the sum of this matrix and the matrix m
         */
        inline TMatrix<T> operator +(const TMatrix<T> &m) const throw (MathException) {

            if (nRows != m.nRows || nColumns != m.nColumns) throw MathException("Error adding matrix:invalid matrix size");

            TMatrix<T> t;
            t.resize(nRows, nColumns, false);

            t = *this;

            t.add(m);

            return t;

        }

        /**
         * Addition and assignment operator.
         * Add to this matrix another matrix.
         * @param matrix a matrix with the same number of rows and coloumns of this matrix
         * @throw MathException if matrices have different size.
         */
        inline TMatrix & operator +=(const TMatrix<T> &matrix) throw (MathException) {
            add(matrix);
            return *this;
        }

        //---------------- subtraction ---------------------------------------

        /**
         * Subtract  matrix from this matrix.
         *
         * @param matrix a matrix with the same number of rows and coloumns of this matrix
         * @throw MathException if matrices have different size.
         */
        inline void sub(const TMatrix<T> &matrix) throw (MathException) {

            if (nRows != matrix.nRows || nColumns != matrix.nColumns) throw MathException("Error sbtracting matrix:invalid matrix size.A matrix of size %d x %d is required as argument", nRows, nColumns);

            const T* pdata = matrix.data;

            for (size_t i = 0; i < nsize; i++) {
                data[i] -= pdata[i];
            }
        }

        /**
         * Fast subtraction.
         * Subtracts the specified matrix from this matrix and places the result in the result matrix.
         * \warning For optimization reasons this method doesn't perform arguments and bounds check.
         * @param m the matrix to add.This matrix must have the same size as this matrix.
         * @param result The matrix where the result of the addition is places.Points to an already allocated matrix of the same size as this matrix.
         */
        inline void sub(TMatrix<T> *result, const TMatrix<T> &matrix) throw () {

            T* dest_data = result->data;

            const T* src_data = matrix->data;

            for (size_t i = 0; i < nsize; i++) {

                dest_data[i] = data[i] - src_data[i];
            }
        }

        /**
         * Subtraction operator.
         * 
         * @param m a matrix with the same number of rows and coloumns of this matrix
         * @return the subtraction between this matrix and  m
         */
        inline TMatrix operator -(const TMatrix<T> &m) const throw (MathException) {

            if (nRows != m.nRows || nColumns != m.nColumns) throw MathException("Error sbtracting matrix:invalid matrix size.A matrix of size %d x %d is required as argument", nRows, nColumns);

            TMatrix<T> t;
            t.resize(nRows, nColumns, false);

            t = *this;

            t.sub(m);

            return t;
        }

        /**
         * Subtraction and assignment operator.
         * @param matrix
         * @return this matrix
         */
        inline TMatrix & operator -=(const TMatrix<T> &matrix) throw (MathException) {
            sub(matrix);
            return *this;
        }

        //multiplication

        /**
         * Matrix multiplication.
         * Multiplies this matrix by another matrix m.
         * @param m a matrix with the same number of rows as the columns of this matrix
         * and a number of colums equals to the number of rows of this matrix.
         */
        inline void mul(const TMatrix<T> &m) throw (MathException) {

            TMatrix<T> temp;
            temp.resize(nRows, nColumns, false);
            temp = (*this) * m;
            data = temp.data;
            temp.data = 0;
        }

        /**
         * Fast matrix multiplication.
         * Multiplies this matrix by m and places the result into an existing matrix.
         * 
         * This method doesn't allocate a temporary matrix and doesn't perform bounds check
         * for optimization reasons.
         * If this matrix is nxm and m is mxq the result matrix is nxq (n=numer of rows q=number of columns)
         * @param result an user supplied matrix of appropriate size to store the result of the multiplication.
         * @param m a matrix with m.nRows=this->nColumns
         */
        inline void mul(TMatrix<T> *result, const TMatrix<T>& m) const throw () {

            const size_t mcols = m.nColumns; //n. colonne

            const T* mdata = m.data;

            T* dest_data = result->data;

            T s;

            size_t k;

            size_t offs=0;
            size_t offs1=0;

            for (size_t row = 0; row < nRows; ++row) {

                offs=row*mcols;
                offs1=row * nColumns;

                for (size_t col = 0; col < mcols; ++col) {

                    s = 0;

                    //prodotto scalare fra la riga row dalla matrice corrente e la colonna col della matrice m
                    for (k = 0; k < nColumns; ++k) s += data[offs1 + k] * mdata[k * mcols + col];

                    dest_data[offs + col] = s;

                }
            }
        }

        /**
         * Matrix multiplication operator.
         * @return this matrix multiplied by another matrix m.
         * @param m a matrix with the same number of columns as this matrix rows.
         * If this matrix is nxm and m is mxq the result matrix is nxq (n=numer of rows q=number of columns)
         */
        TMatrix<T> operator *(const TMatrix<T> &m) throw (MathException) {

            //il numero di colonne di questa matrice deve essere uguale al numero di righe di m
            if (m.nRows != nColumns) throw MathException("Matrix multiplication error:invalid matrix size. A matrix with %d rows is required as argument.", nColumns);

            const size_t mcols = m.nColumns; //n. colonne

            const T* mdata = m.data;

            TMatrix<T> res;

            res.resize(nRows, m.nColumns, false);

            T* dest_data = res.data;

            T s;

            size_t k;

            size_t offs=0;
            size_t offs1=0;

            for (size_t row = 0; row < nRows; ++row) {

                offs=row*mcols;
                offs1=row * nColumns;

                for (size_t col = 0; col < mcols; ++col) {

                    s = 0;

                    //prodotto scalare fra la riga row dalla matrice corrente e la colonna col della matrice m
                    for (k = 0; k < nColumns; ++k) s += data[offs1 + k] * mdata[k * mcols + col];

                    dest_data[offs + col] = s;

                }
            }

            return res;
        }

        /**
         * Multiplication and assignment operator
         * Multiplies this matrix by another matrix m.
         * @param m
         * @return 
         */
        TMatrix<T> & operator *=(const TMatrix<T> &m) throw (MathException) {
            mul(m);
            return *this;
        }

        /**
         * Tests if this matrix equals another matrix.
         * @param m
         * @return true if the matrices have the same size and each element a[i][j] equals the corresponding element.
         */
        bool operator ==(const TMatrix<T> &m) const throw () {

            if (m.getRows() != nRows || m.getColumns() != nColumns) return false;

            const T *pdata = m.data;

            for (size_t i = 0; i < nsize; i++) {
                if (pdata[i] != data[i]) return false;
            }

            return true;
        }

        /**
         * Tests if this matrix is not equal to another matrix.
         * @param m
         * @return
         */
        bool operator !=(const TMatrix<T> &m) throw () {
            return !(m == *this);
        }

        /**
         * Unary minus operator.
         * Negate all elements of the matrix
         * @return
         */
        TMatrix<T> & operator -() throw () {
            for (size_t i = 0; i < nsize; i++) data[i] = -data[i];
            return *this;
        }

        /**
         *Returns the transpose matrix of this matrix.
         *
         */
        TMatrix<T> getTranspose() const throw () {

            TMatrix<T> t;

            t.resize(nColumns, nRows, false);

            for (size_t i = 0; i < nRows; i++) {
                for (size_t j = 0; j < nColumns; j++) {
                    t[j][i] = data[i * nColumns + j];
                }
            }

            return t;
        }

        /**
         * Matrix by vector multiplication.
         * @param v a vector whose size equals this matrix columns count.
         * @throw MathException if the vector v has size different fron nColumns.
         * @return a vector of size nRows resulting from the multiplication.
         */
        TVector<T> operator *(const TVector<T> &v) const throw (MathException) {

            if (v.size() != nColumns) throw MathException("Error in matrix by vector multiplication:invalid vector size.The vector must have size %d", nColumns);

            TVector<T> res(nRows,false);

            const T *src = v._getDataPtr();
            T* dest=res._getDataPtr();

            size_t offs=0;
            
            T temp;

            for (size_t row = 0; row < nRows; row++) {

                temp=T();

                offs=row * nColumns;

                for (size_t col = 0; col < nColumns; col++) {

                    temp += data[offs + col] * src[col];
                }

                dest[row]=temp;
            }

            return res;
        }

        /**
         * Fast matrix by vector multiplication.
         * This method returns the result of the multiplication of this matrix by a vector
         * v and places the result in a user-supplied vector.
         * @param result the vector already allocated that holds the result.Must have nRows size.
         * @param v a vector whose size equals this matrix columns count.
         * \warning for optimization reasons, this method doesn't perform bounds check and arguments check.
         */
        void mul(TVector<T> *result, const TVector<T> &v) const throw () {

            const T *src = v._getDataPtr();
            T* dest=result->_getDataPtr();

            size_t offs=0;
            T temp;

            for (size_t row = 0; row < nRows; row++) {

                T *r = result->_getDataPtr();

                temp=T();

                offs=row*nColumns;

                for (size_t col = 0; col < nColumns; col++) {

                    temp += data[offs + col] * src[col];
                }

                dest[row]=temp;
            }
        }

        /**
         *Extracts a row from this matrix into dest_row.
         *@param dest_row a vector of nRows size already allocated.
         *@param row the index of the row to extract.(first row has index zero)
         *\warning for optimization reasons, this method doesn't perform bounds check and arguments check.
         */
        void getRow(TVector<T> *dest_row, size_t row) const throw () {
            T* dest = dest_row->_getDataPtr();
            memcpy(dest, &data[row * nColumns], sizeof (T) * nColumns);
        }

        /**
         *Extracts a row from this matrix into dest_row.
         *@param dest_column a vector of nColumns size already allocated.
         *@param col the index of the column to extract.
         *\warning for optimization reasons, this method doesn't perform bounds check and arguments check.
         */
        void getColumn(TVector<T> *dest_column, size_t column) const throw () {

            T* dest = dest_column->_getDataPtr();

            for (size_t i = 0; i < nColumns; i++) {

                dest[i] = data[i * nColumns + column];

            }
        }

        /**
         * Sets a row of the matrix using an array of elements.
         * @param row row index to set
         * @param row_data and array,not null,containing nColumns elements.
         * \warning this method doesn't perform arguments check
         */
        void setRow(size_t row, const T* row_data) throw () {

            memcpy(&data[row * nColumns], row_data, sizeof (T) * nColumns);

        }

        /**
         * Sets the firs nelements of a row
         * @param row row index (first row has index zero)
         * @param row_data and array of size nelements
         * @param nelements number of elements to set, must be less or equal to nColumns
         * \warning this method doesn't perform arguments check
         */
        void setRow(size_t row, const T* row_data, size_t nelements) {

            memcpy(&data[row * nColumns], row_data, sizeof (T) * nelements);
        }

        /**
         * Sets a row of the matrix using a vector
         * @param row row index to set
         * @param row_data a vector of size nColumns
         * @throw MathException if row_data is not of size nColumns or if row is out of bounds
         */
        void setRow(size_t row, const TVector<T> &row_data) throw (MathException) {

            if (row_data.size() != nColumns) throw MathException("Error setting row %d:invalid array size, expected size %d, actual %d", row, nColumns, row_data.size());

            if (row < 0 || row >= getRows()) throw MathException("Error setting row %d:row index out of bounds (must be between 0 and %d)", row, getRows() - 1);

            memcpy(&data[row * nColumns], row_data._getDataPtr(), sizeof (T) * nColumns);

        }

        /**
         * Sets a column of the matrix using a vector.
         * @param column index of column (first column has index=zero)
         * @param column_data a vector of size nRows
         * @throw MathException if column_data has size different from nRows of column index is out of bounds
         */
        void setColumn(size_t column, const TVector<T> &column_data) throw (MathException) {

            if (column_data.size() != getRows()) throw MathException("Error setting column %d:invalid array size, expected size %d, actual %d", column, getRows(), column_data.size());

            if (column < 0 || column >= getColumns()) throw MathException("Error setting column %d:column index out of bounds (must be between 0 and %d)", column, getColumns() - 1);

            const T* cdata = column_data._getDataPtr();

            for (size_t i = 0; i < nRows; i++) {
                data[i * nColumns + column] = cdata[i];
            }
        }

        /**
         * Sets a column of the matrix using an array of elements.
         * @param column index of column (first column has index=zero)
         * @param column_data array of nRows elements
         * \warning this method doesn't perform arguments check
         */
        void setColumn(size_t column, const T* column_data) throw () {

            for (size_t i = 0; i < nRows; i++) {
                data[i * nColumns + column] = column_data[i];
            }
        }

        /**
         * @todo da testare!
         * Swaps two rows
         * @param row1 first row
         * @param row2 second row
         * \waring this method doesn't perform row index validation.
         */
        void swapRows(size_t row1, size_t row2) throw () {

            if (row1 == row2) return;

            T temp;
            T* ptr;

            size_t dest_offs = nColumns*row2; //offset destinazione
            size_t src_offs = nColumns*row1; //offset sorgente

            for (size_t i = 0; i < nColumns; i++) {

                temp = data[dest_offs];
                data[dest_offs] = data[src_offs];
                data[src_offs] = temp;

                ++src_offs;
                ++dest_offs;
            }
        }

        /**
         * Swaps two columns
         * @param col1 first column
         * @param col2 second column
         * \waring this method doesn't perform row index validation.
         */
        void swapColumns(size_t col1, size_t col2) throw () {

            if (col1 == col2) return;

            T temp;
            T* ptr;

            size_t dest_offs = col2; //offset destinazione
            size_t src_offs = col1; //offset sorgente

            for (size_t i = 0; i < nRows; i++) {

                temp = data[dest_offs];
                data[dest_offs] = data[src_offs];
                data[src_offs] = temp;

                dest_offs += nColumns;
                src_offs += nColumns;

            }
        }

        /**
         * Copy a submatrix of this matrix to a destination matrix
         * \warning this method doesn't perform bounds check if clip is not set
         * @param dest_matrix the destination matrix,must have appropriate size
         * @param dest_row destination row where data are copied
         * @param dest_col destination col where data are copied
         * @param src_row first source row
         * @param src_col first source column
         * @param end_src_row last source row of the submatrix
         * @param end_src_col last source column of the submatrtix
         * @param clip if set,clip the row indices according to the metrices' size
         */
        void copySubMatrixTo(TMatrix<T>& dest_matrix, int dest_row, int dest_col, int src_row, int src_col, int end_src_row, int end_src_col, bool clip) throw () {

            if (clip) {

                //clipping
                if (dest_row < 0) dest_row = 0;
                else if (dest_row >= dest_matrix.nRows) dest_row = dest_matrix.nRows - 1;

                if (dest_col < 0) dest_col = 0;
                else if (dest_col >= dest_matrix.nColumns) return;

                if (src_row<0) src_row=0;
                else if (src_row>= nRows) src_row=nRows-1;

                if (src_col<0) src_col=0;
                else if (src_col>=nColumns) src_col=nColumns-1;

                if (end_src_row>=nRows) end_src_row=nRows-1;
                else if (end_src_row<0) return;

                if (end_src_col>=nColumns) end_src_col=nColumns-1;
                else if (end_src_col<0) return;

                if ((dest_row+end_src_row-src_row)>=dest_matrix.nRows) {
                    end_src_row=dest_matrix.nRows-1+src_row-dest_row;
                }

                if ((dest_col+end_src_col-src_col)>=dest_matrix.nColumns) {
                    end_src_col=dest_matrix.nColumns-1+src_col-dest_col;
                }
            }

            int wc = end_src_col - src_col;
            size_t sz = sizeof (T)*(wc + 1);

            T* psrc = 0;
            T* pdest = 0;

            int dr = dest_row;
            int dncols = dest_matrix.getColumns();

            for (int r = src_row; r <= end_src_row; r++) {

                psrc = &data[r * nColumns];
                pdest = &dest_matrix.data[dr * dncols + dest_col];

                memcpy(pdest, psrc, sz);

                ++dr;

            }

        }

        /**
         * Extract a submatrix of size m_rows x n_cols from this matrix
         * @param src_row index of source row in this matrix
         * @param src_col index of source coloumn in this matrix
         * @param n_rows number of rows of the submatrix
         * @param n_cols number of columns of the submatrix
         * @return a submatrix m_rows x n_cols
         * @thorw MathException if indices and/or n_rows and n_cols are not compatible with
         * the size of this matrix.
         */
        TMatrix<T> getSubMatrix(int src_row,int src_col,size_t n_rows,size_t n_cols) throw (MathException) {

            if (src_row<0 || src_row>=nRows || src_col<0 || src_col>=nColumns) throw MathException("Error extracting a submatrix:invalid row or column indices.");
            if ((src_row+n_rows)>nRows || (src_col+n_cols)>nColumns) throw MathException("Error extracting a submatrix of size %d x %d:invalid size ",n_rows,n_cols);

            TMatrix<T> dest(n_rows,n_cols);

            copySubMatrixTo(dest,0,0,src_row,src_col,src_row+n_rows-1,src_col+n_cols-1,false);

            return dest;
        }

        /**
         * Makes this matrix a band matrix.
         * @param lowerDiags number of diagonals below the main diagonal (must be >=0)
         * @param upperDiags number of diagonals above the main diagonal (must be >=0)
         * @throws MathException if this matrix is not a square matrix
         */
        void makeBandMatrix(const int lowerDiags,const int upperDiags) throw (MathException) {

            if (!isSquare()) throw MathException("Error making the band matrix:matrix must be square.");

            if (lowerDiags>nRows-1 || lowerDiags<0 || upperDiags>nRows-1 || upperDiags<0) throw MathException("Error making the band matrix:invalid params.");

            size_t offs=0;

            for (int i=0;i<nRows;i++) {

                int u=i-lowerDiags;

                offs=i*nColumns; //matrix memory offset

                for (int j=0;j<u;j++) {
                    data[offs + j]=T();
                }

                for (int j=i+upperDiags+1;j<nColumns;j++) {
                    data[offs+j]=T();
                }
            }
        }

        /**
         * Makes this matrix a lower triangular matrix
         * @throws MathException if this matrix is not a square matrix
         */
        void makeLowerTriangular() throw (MathException) {

            if (!isSquare()) throw MathException("Error making the triangular matrix:matrix must be square.");

            size_t offs=0;

            for (int i=0;i<nRows;i++) {

                offs=i*nColumns;

                for (int j=i+1;j<nColumns;j++) {
                    data[offs+j]=T();
                }
            }

        }

        /**
         * Makes thie matrix an upper triangular matrix
         * @throws MathException if this matrix is not a square matrix
         */
        void makeUpperTriangular() throw (MathException) {

            if (!isSquare()) throw MathException("Error making the triangular matrix:matrix must be square.");

            size_t offs=0;

            for (int i=0;i<nRows;i++) {

                offs=i*nColumns;

                for (int j=0;j<i;j++) {
                    data[offs+j]=T();
                }
            }
        }

        /**
         * Internal data pointer.
         * \warning don't use this method to access the elements.
         * @return the pointer to the internal matrix data arranged in row major format.
         */
        T* _getDataPtr() const throw () {
            return data;
        }

        /**
         * Pointer to row data
         * \warning don't use this method to access the elements.
         * @param row index of row
         * @return the pointer to row data
         */
        T* _getRowPtr(size_t row) const throw () {
            return &data[row * nColumns];
        }

        /**
         * Set directly the internal data pointer.
         * \warning This method is unsafe.Don't use this method in normal computation.
         * This is a special method used to share data among matricesand avoid data duplication
         * with very large matrices.
         * @param pdata a pointer to an array with nRows*nColumns elements (row major format).
         */
        void _setDataPtr(T* pdata) throw () {
            data = pdata;
        }

        

    protected:

        TMatrix() : data(0) {

        }

        /**
         *Number of rows.
         */
        size_t nRows;
        /**
         *Number of columns.
         */
        size_t nColumns;

        /**
         *Number of elements in the matrix (nxm)
         */
        size_t nsize;

        /**
         * TMatrix components (row major format)
         * Data store the elements of the matrix arranged
         * by rows
         */
        T* data;

    };

    ////////////////////////////////////////////////////////////////////////////   

    /**
     * Output operator
     * @param ostrm
     * @param m
     * @return
     */
    template<class T> ostream & operator <<(ostream& ostrm, const TMatrix<T> &m) {

        size_t r = m.getRows();
        size_t c = m.getColumns();

        for (size_t i = 0; i < r; i++) {

            for (size_t j = 0; j < c; j++) {
                ostrm << m[i][j] << '\t';
            }
            ostrm << endl;
        }
        return ostrm;
    }


};

#endif


