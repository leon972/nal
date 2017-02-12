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
 * \file bandmatrix.h
 * \brief Generic square template-base band matrix.
 */

#ifndef BANDMATRIX_H
#define	BANDMATRIX_H

#include "math/mathdefs.h"
#include "libmath.h"
#include "tvector.h"
#include "tmatrix.h"
#include <stddef.h>
#include <iostream>

namespace mathengine {

    /**
     * unused constant:doxygen trick
     */
    const int band_matrix_zero = 0;

    /**
     * A square nxn band matrix.     
     */
    template<typename T> class BandMatrix {
    public:

        /**
         * Create a square nxn band matrix.
         * To create a tridiagonal matrix h and l must be set to 1
         * @param n size of the matrix (number of rows and columns)
         * @param l number of diagonals below the main diagonal
         * @param h number of diagonals above the main diagonal
         * @param reset_memory if set, the matrix internal memory is set to zero
         * @throws MathException if n,l,h are not congruent.
         */
        BandMatrix(size_t n, size_t l, size_t h, bool reset_memory) throw (MathException) {

            init(n, l, h, reset_memory);

        }

        /**
         * Create a square nxn band matrix.
         * To create a tridiagonal matrix h and l must be set to 1
         * @param n size of the matrix (number of rows and columns)
         * @param l number of diagonals below the main diagonal
         * @param h number of diagonals above the main diagonal
         * @throws MathException if n,l,h are not congruent.
         */
        BandMatrix(size_t n, size_t l, size_t h) throw (MathException) {
            init(n, l, h, true);
        }

        virtual ~BandMatrix() {
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
        BandMatrix<T> & operator=(const BandMatrix<T> &m) {

            if (m.n == n) {

                if (&m == this) return *this;

                memcpy(data, m.data, sizeof (T) * nsize);

            } else {

                resize(m.getRows(), m.getLdiag(), m.getHdiag(), false);
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
         * Resize this matrix
         * @param n size of the matrix (number of rows and columns)
         * @param l number of diagonals below the main diagonal
         * @param h number of diagonals above the main diagonal
         * @param reset_memory if set, the matrix internal memory is set to zero
         * @throws MathException if n,l,h are not congruent.
         */
        void resize(size_t n, size_t l, size_t h, bool reset_memory) {

            checkParams(n, l, h);

            if (data) {
                delete[] data;
                data = 0;
            }

            init(n, l, h, reset_memory);

        }

        /**
         * Number of rows
         * @return the number of rows
         */
        size_t getRows() const {
            return n;
        }

        /**
         * Number of columns
         * @return the number of columns
         */
        size_t getColumns() const {
            return n;
        }

        /**
         * Get the number of diagonals below the main diagonal
         * @param r
         * @param c
         * @return a value>=0
         */
        size_t getLdiag() const {
            return l;
        }

        /**
         * Get the number of diagonals above the main diagonal
         * @return a value>=0
         */
        size_t getHdiag() const {
            return h;
        }

        /**
         * Subscript operator, can be used to access the matrix elements.
         * This method can be used to modify elements in the matrix.
         * eq. T &p=matrix(i,j);
         * p->someMethod(...)
         * \warning this method doesn't perform bounds check for optimization reasons.
         * Calling this method with r and c outside the valid area of the matrix will cause
         * unexpected behaviour
         * @param r row (first row has index 0)
         * @param c column (first column has index 0)
         * @return the reference of the element at (r,c)
         *
         */
        inline T & operator ()(size_t r, size_t c) throw () {
            //posizione nell'area di memoria
            return data[r * pitch + c - r + l];
        }

        /**
         * Subscript operator (const version), can be used to access the matrix elements
         * \warning this method doesn't perform bounds check for optimization reasons.
         * Calling this method with r and c outside the valid area of the matrix will cause
         * unexpected behaviour
         * @param r row
         * @param c column
         * @return the element at (r,c)
         *
         */
        inline T operator()(size_t r, size_t c) const throw () {
            return data[r * pitch + c - r + l];
        }

        /**
         * Get an element at the specified row and column position.
         *
         * @param row index of row (first row has index zero)
         * @param column index of column (first column has index zero)
         * @return the reference of the element at the specified row and column.
         * If r and c are outside the band, returns the default value of type T
         * @throw MathException if row and/or column are >=n or <0
         */
        T getElementAt(size_t r, size_t c) const throw (MathException) {

            int x = c - r + l;

            if (r >= n) {
                MathException ex("Error accessing matrix element:row %d out of bounds.", r);
                throw ex;
            }

            if (c >= n) {
                MathException ex("Error accessing matrix element:column %d out of bounds.", c);
                throw ex;
            }

            if (c > r + h || r > c + l) return T();

            return data[r * pitch + c - r + l];
        }

        /**
         * Set an element of the matrix with indices bound check.
         * If r and c are outside the band, this method returns without throwing an exception.
         * @param r index of row (first row has index zero) must be < c+ l
         * @param c index of column (first column has index zero) must be < r+ h
         * @param value the new value of the element
         * @throw MathException if row and/or column are out of bounds.
         */
        void setElementAt(size_t r, size_t c, const T& value) throw (MathException) {

            int x = c - r + l;

            if (r >= n) {
                MathException ex("Error setting matrix element:row %d out of bounds.", r);
                throw ex;
            }

            if (c >= n) {
                MathException ex("Error setting matrix element:column %d out of bounds.", c);
                throw ex;
            }

            if (c > r + h || r > c + l) return; //non si possono settare elementi fuori dalla banda di valori

            data[r * pitch + c - r + l] = value;
        }

        /**
         * Determine if this matrix has the same number of rows and columns of m
         * @param m a matrix
         */
        inline bool hasSameSize(const BandMatrix<T> &m) const throw () {
            return n == m.n && m.getHdiag() == h && m.getLdiag() == l;
        }

        /**
         * Determine if this matrix has the same number of rows and columns of m
         */
        inline bool hasSameSize(const TMatrix<T> &m) const throw () {
            return m.getRows() == n && m.getColumns() == n;
        }

        /**
         * Clear internal memory (all bytes set to zero)
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

            for (size_t i = 0; i < n; i++) {
                data[i * pitch + l] = value;
            }
        }

        /**
         * Copy the band elements from a square matrix
         * @param n
         * @param l
         * @param h
         * @param reset_memory
         */
        void extractFrom(const TMatrix<T>& m) throw (MathException) {

            if (!hasSameSize(m)) throw MathException("Cannot copy band matrix from generic matrix:different matrix sizes");

            for (int i = 0; i < nsize; i++) {

                int r = i / pitch;
                int c = i - r * pitch + r - l;

                if (c >= 0 && c < n) {
                    data[i] = m[r][c];
                }
            }

        }

        /**
         * Copies this matrix to a square matrix
         * @param m a square matrix with the same size as this matrix
         * @throws Math exception if m is not square or has different size.
         */
        void copyTo(TMatrix<T>& m) const throw (MathException) {

            if (m.getRows() < n || m.getColumns()<n) throw MathException("Cannot copy band matrix:invalid destination matrix size");

            m.zeroMemory();

            for (int i = 0; i < nsize; i++) {

                int r = i / pitch;
                int c = i - r * pitch + r - l;

                if (c >= 0 && c < n) {
                    m[r][c] = data[i];
                }
            }
        }

        //--------------- addition --------------------------------------------

        /**
         * Add m to this matrix
         * @param m another matrix
         * @thorw MathException if m ha different size from this matrix
         */
        inline void add(const BandMatrix<T> &m) throw (MathException) {

            if (!hasSameSize(m)) throw MathException("Error adding band matrix:invalid matrix size");

            const T* pdata = m.data;

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
        inline void add(BandMatrix<T> *result, const BandMatrix<T> &m) const throw () {

            T* dest_data = result->data;

            const T* src_data = m.data;

            for (size_t i = 0; i < nsize; i++) {

                dest_data[i] = data[i] + src_data[i];
            }
        }

        /**
         * Addition operator
         *
         * @param m a matrix with the same number of rows and coloumns of this matrix
         * @throw MathException if matrices have different size or different band sizes.
         * @return the sum of this matrix and the matrix m
         */
        inline BandMatrix<T> operator +(const BandMatrix<T> &m) const throw (MathException) {

            if (!hasSameSize(m)) throw MathException("Error adding matrix:invalid matrix size");

            BandMatrix<T> t(n, l, h, false);

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
        inline BandMatrix<T> & operator +=(const BandMatrix<T> &matrix) throw (MathException) {
            add(matrix);
            return *this;
        }

        //---------------- subtraction -----------------------------------------

        /**
         * Subtract m from this matrix
         * @param m another matrix
         * @thorw MathException if m ha different size from this matrix
         */
        inline void sub(const BandMatrix<T> &m) throw (MathException) {

            if (!hasSameSize(m)) throw MathException("Error subtractiong band matrix:invalid matrix size");

            const T* pdata = m.data;

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
        inline void sub(BandMatrix<T> *result, const BandMatrix<T> &m) const throw () {

            T* dest_data = result->data;

            const T* src_data = m.data;

            for (size_t i = 0; i < nsize; i++) {

                dest_data[i] = data[i] - src_data[i];
            }
        }

        /**
         * Subtraction operator
         *
         * @param m a matrix with the same number of rows and coloumns of this matrix
         * @throw MathException if matrices have different size or different band sizes.
         * @return the difference of this matrix and the matrix m
         */
        inline BandMatrix<T> operator -(const BandMatrix<T> &m) const throw (MathException) {

            if (!hasSameSize(m)) throw MathException("Error subtracting matrix:invalid matrix size");

            BandMatrix<T> t(n, l, h, false);

            t = *this;

            t.sub(m);

            return t;

        }

        /**
         * Addition and assignment operator.
         * Add to this matrix another matrix.
         * @param matrix a matrix with the same number of rows and coloumns of this matrix
         * @throw MathException if matrices have different size.
         */
        inline BandMatrix<T> & operator -=(const BandMatrix<T> &matrix) throw (MathException) {
            sub(matrix);
            return *this;
        }


        //-------------------- multiplication ----------------------------------

        /**
         * Multiplies this band matrix by band matrix m
         * \warning this method is slow.Band matrices are converted to normal
         * square matrix before multiplication.
         * @param m a band matrix
         * @return a square matrix
         */
        TMatrix<T> operator *(const BandMatrix<T>& m) const throw (MathException) {

            if (m.getRows() != n) {
                throw MathException("Error in band matrix multiplication:invalid matrix size");
            }

            TMatrix<T> m1(n, n, false);
            TMatrix<T> m2(n, n, false);

            this->copyTo(m1);
            m.copyTo(m2);

            m1.mul(m2);

            return m1;
        }

        /**
         * Band matrix by vector multiplication.
         * Multiplies this matrix by vector v
         * @param v
         * @return
         */
        TVector<T> operator *(const TVector<T>& v) const throw (MathException) {

            if (v.size() != n) throw MathException("Error in band matrix by vector multiplication:invalid vector size");

            int s = 0;
            int f = 0;

            TVector<T> vr(n, false);

            size_t offs = 0;

            for (int i = 0; i < n; i++) {

                s = i - l;
                if (s < 0) s = 0;
                f = i + h + 1;
                if (f > n) f = n;

                T temp = T();

                offs = i * pitch;

                for (int j = s; j < f; j++) {
                    temp = temp + v[j] * data[offs + j - i + l];
                }

                vr[i] = temp;
            }

            return vr;
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
         * Number of elements in the internal marix memory
         * @return
         */
        size_t _getDataSize() const throw () {
            return nsize;
        }

        /**
         * Width of data memory (number of elements per row)
         * @return
         */
        size_t _getPitch() const throw() {
            return pitch;
        }

        /**
         * Write a band matrix to file
         * Use cout.precision(n) to adjust the precision of the values printed
         * @param mat a band matrix to write to an output file
         * @param file_name destination file full path
         * @param separator_char separator char used to separate the elements of the matrix in each row.
         * @param quot if set,values are quoted
         * @param comment if not null,writes comment at the beginning of the file.
         * @param writeCompact if set,zero elements are omitted,writes only the band elements.
         * @param prec precision used for output floating point values (internally uses file.precision(prec))
         */
        static void writeToFile(const BandMatrix<T>& mat, const char* file_name, const char separator_char, bool quot, const char* comment, bool writeCompact,streamsize prec) throw (MathException) {

            fstream file;

            //apre in scrittura,il file è ricreato
            file.open(file_name, fstream::out);

            const T* pdata = mat._getDataPtr();

            if (file.good()) {

                if (prec>0) file.precision(prec);

                size_t nr = mat.getRows();
                size_t nc = mat.getColumns();
                size_t pitch = mat.pitch;

                if (comment) {
                    file << '#' << comment << endl;
                }

                if (!quot) {

                    if (writeCompact) {
                        //omette gli elementi fuori dalla banda
                        for (int i = 0; i < nr; i++) {

                            size_t offs = i*pitch;

                            for (int j = 0; j < pitch - 1; j++) {
                                file << pdata[offs + j] << separator_char;
                            }

                            file << pdata[offs + pitch - 1] << endl;

                        }

                    } else {

                        for (int i = 0; i < nr; i++) {

                            for (int j = 0; j < nc - 1; j++) {
                                file << mat.getElementAt(i, j) << separator_char;
                            }

                            file << mat.getElementAt(i, nc - 1) << endl;
                        }

                    }

                } else {

                    const char quot = '"';

                    if (writeCompact) {

                        //valori separati da "
                        for (int i = 0; i < nr; i++) {

                            size_t offs = i*pitch;

                            for (int j = 0; j < pitch - 1; j++) {
                                file << quot << pdata[offs + j] << quot << separator_char;
                            }

                            file << quot << pdata[offs + pitch - 1] << quot << endl;
                        }

                    } else {

                        for (int i = 0; i < nr; i++) {

                            for (int j = 0; j < nc - 1; j++) {
                                file << quot << mat.getElementAt(i, j) << quot << separator_char;
                            }

                            file << quot << mat.getElementAt(i, nc - 1) << quot << endl;
                        }
                    }

                    file.close();
                }

            } else {
                file.close();
                throw MathException("Error writing the %d x %d matrix to file %s:file open or creation error.", mat.getRows(), mat.getColumns(), file_name);
            }

        }

        /**
         * Writes this matrix to a text file.Each row of the matrix is written in a separate row and the elements are separated by the separator_char
         * This method will overwrite the destination file if it exists.
         * @param separator_char the character used to separate the elements, eq. ',' for comma separated values.
         * @param quot if set,each element is enclosed in double quote (")
         * @param comment if not null, the comment will be written in the first row
         * @param writeCompact if set,writes this matrix in compact form (elements outside the band are omitted)
         * @param prec precision used for output floating point values (internally uses file.precision(prec))
         * @throw MathException if an I/O error occurs
         */
        void writeToFile(const char* file_name, const char separator_char, bool quot, const char* comment, bool writeCompact,streamsize prec) throw (MathException) {

            writeToFile(*this, file_name, separator_char, quot, comment, writeCompact,prec);

        }

        /**
         * Reads this matrix from a text file
         * @param separator_char the character used to separate the elements, eq. ',' for comma separated values.
         * @param quot if set,each element is enclosed in double quote (")
         * @param readCompact must be set to true if the matrix is stored in compact form.
         */
        void readFromFile(const char* file_name, const char separator_char, bool quot, bool readCompact) throw (MathException) {

            fstream file;

            //apre in lettura
            file.open(file_name, fstream::in);

            int line_cnt = 1;
            int crow = 0;
            int ccol = 0;

            if (file.good()) {

                char ch;

                while (!file.eof()) {

                    file.get(ch);

                    if (ch == '#' || ch == '\n') {

                        while (file.get() != '\n' && !file.eof()) {
                        }

                        line_cnt++;

                    } else {

                        file.unget();

                        T val;

                        if (quot) {

                            //valori fra apici ("")
                            if (readCompact) {

                                for (int i = 0; i < pitch - 1; i++) {

                                    file >> ch;

                                    if (ch != '"') {
                                        throw MathException("Error reading value from file %s at line %d at column %d:missing opening quotes", file_name, line_cnt, i);
                                    }

                                    if (!(file >> val)) {

                                        throw MathException("Error reading value from file %s at line %d at column %d:invalid format", file_name, line_cnt, i);
                                    }

                                    file >> ch;

                                    if (ch != '"') {
                                        throw MathException("Error reading value from file %s at line %d at column %d:missing closing quotes", file_name, line_cnt, i);
                                    }


                                    //setta il valore letto nella posizione corrente
                                    data[crow * pitch + i] = val;

                                    file >> ch;

                                    if (ch != separator_char) {
                                        throw MathException("Error reading value from file %s at line %d at column %d:invalid separator", file_name, line_cnt, i);
                                    }
                                }

                            } else {

                                size_t cn = getColumns() - 1;

                                for (int i = 0; i < cn; i++) {

                                    file >> ch;

                                    if (ch != '"') {
                                        throw MathException("Error reading value from file %s at line %d at column %d:missing opening quotes", file_name, line_cnt, i);
                                    }

                                    if (!(file >> val)) {

                                        throw MathException("Error reading value from file %s at line %d at column %d:invalid format", file_name, line_cnt, i);
                                    }

                                    file >> ch;

                                    if (ch != '"') {
                                        throw MathException("Error reading value from file %s at line %d at column %d:missing closing quotes", file_name, line_cnt, i);
                                    }

                                    if (i >= crow - l && i <= crow + h) {
                                        //setta il valore letto nella posizione corrente
                                        data[crow * pitch + i - crow + l] = val;
                                    }

                                    file >> ch;

                                    if (ch != separator_char) {
                                        throw MathException("Error reading value from file %s at line %d at column %d:invalid separator", file_name, line_cnt, i);
                                    }
                                }//fine for
                            }

                        } else {

                            //valori senza apici

                            if (readCompact) {

                                for (int i = 0; i < pitch - 1; i++) {

                                    if (!(file >> val)) {

                                        throw MathException("Error reading value from file %s at line %d at column %d:invalid format", file_name, line_cnt, i);
                                    }

                                    //setta il valore letto nella posizione corrente
                                    data[crow * pitch + i] = val;

                                    file >> ch;

                                    if (ch != separator_char) {
                                        throw MathException("Error reading value from file %s at line %d at column %d:invalid separator", file_name, line_cnt, i);
                                    }
                                }

                            } else {

                                //valori senza apici formato non compatto
                                int cn = getColumns() - 1;

                                for (int i = 0; i < cn; i++) {

                                    if (!(file >> val)) {

                                        throw MathException("Error reading value from file %s at line %d at column %d:invalid format", file_name, line_cnt, i);
                                    }

                                    if (i >= crow - l && i <= crow + h) {
                                        //setta il valore letto nella posizione corrente
                                        data[crow * pitch + i - crow + l] = val;
                                    }

                                    file >> ch;

                                    if (ch != separator_char) {
                                        throw MathException("Error reading value from file %s at line %d at column %d:invalid separator", file_name, line_cnt, i);
                                    }
                                }//fine for
                            }
                        }

                        int i = readCompact ? pitch - 1 : getColumns() - 1;

                        //ultimo valore riga corrente
                        if (quot) {

                            file >> ch;
                            if (!(file >> val)) {
                                throw MathException("Error reading value from file %s at line %d at column %d:invalid format", file_name, line_cnt, i);
                            }
                            file >> ch;

                        } else {
                            if (!(file >> val)) {
                                throw MathException("Error reading value from file %s at line %d at column %d:invalid format", file_name, line_cnt, i);
                            }
                        }

                        ch = file.get();

                        while (ch != '\n' && !file.eof()) {
                            ch = file.get();
                        }

                        if (readCompact) {
                            data[crow * pitch + i] = val;
                        } else {
                            if (i >= crow - l && i <= crow + h) {
                                //setta il valore letto nella posizione corrente
                                data[crow * pitch + i - crow + l] = val;
                            }
                        }

                        crow++;
                        line_cnt++;

                        if (crow == getRows()) break;
                    }
                }

                if (crow != getRows()) {
                    throw MathException("Error reading band matrix from file %s:invalid row count, expected %d actual %d", file_name, getRows(), crow);
                }

                file.close();

            } else {

                file.close();
                throw MathException("Error reading the %d x %d band matrix from file %s:file open or creation error.", getRows(), getColumns(), file_name);
            }
        }

    private:

        /**
         * Initialise the matrix memory and parameters
         *
         * Meaning of l and h parameters:
         *
         *  d,x,x,0,0,0
         *  y,d,x,x,0,0
         *  0,y,d,x,x,0
         *  0,0,y,d,x,x
         *  0,0,0,y,d,x
         *  0,0,0,0,y,d
         *
         * in the above example, l=1 and h=2
         *
         * @param n number of rows and columns
         * @param l number of diagonal rows below the main diagonal
         * @param h number of diagonal rows above the main diagonal
         * @param reset_memory if set sets all elements to zero
         */
        void init(size_t n, size_t l, size_t h, bool reset_memory) throw (MathException) {

            checkParams(n, l, h);

            this->n = n;
            this->l = l;
            this->h = h;

            pitch = 1 + h + l;

            nsize = pitch*n;

            data = new T[nsize];

        }

        /**
         *Check size parameters
         */
        void checkParams(size_t n, size_t l, size_t h) throw (MathException) {

            if (n < 1) throw MathException("Band matrix bad size");

            if (h > n - 1 || l > n - 1) throw MathException("Band matrix bad params.");


        }

        size_t n; //size of the matrix
        int l; //number of diagonals below the main diagonal
        int h; //number of diagonals above the main diagonal
        size_t pitch; //pitch of data memory
        size_t nsize; //mem size

        T* data; //internal data memory (matrix is stored in compact form to save memory)
    };

    ////////////////////////////////////////////////////////////////////////////

    /**
     * Output operator
     * @param ostrm
     * @param m
     * @return
     */
    template<class T> ostream & operator <<(ostream& ostrm, const BandMatrix<T> &m) {

        size_t r = m.getRows();

        for (size_t i = 0; i < r; i++) {

            for (size_t j = 0; j < r; j++) {
                ostrm << m.getElementAt(i, j) << '\t';
            }
            ostrm << endl;
        }
        return ostrm;
    }

};



#endif	/* BANDMATRIX_H */

