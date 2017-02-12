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
 * \file tmatrixex.h
 * \brief Extended TMatrix<T> template-base matrix class.
 * This class implements additional methods to compute the inverse matrix,
 * the determinant and others method to read and save the matrix on files.
 */

#ifndef TMATRIXEX_H
#define	TMATRIXEX_H

#include "math/mathdefs.h"
#include "libmath.h"
#include "mathutil.h"
#include "tmatrix.h"
#include "crout.h"
#include <stddef.h>
#include <iostream>
#include <fstream>
#include <sstream>

namespace mathengine {

    using namespace std;

    /**
     * unused constant:doxygen trick
     */
    const int TMatrixEx_zero = 0;

    /**
     * Extended version of TMatrix<T> class.
     * This class implements additional methods to compute the inverse matrix,
     * the determinant and others method to read and save the matrix on files.
     * @see TMatrix
     */
    template <class T> class TMatrixEx : public TMatrix<T> {
    public:

        /**
         * Creates a matrix nrows x ncolumns
         * @param nrows
         * @param ncolumns
         */
        TMatrixEx(size_t nrows, size_t ncolumns) : TMatrix<T>(nrows, ncolumns) {


        }

        /**
         * Creates a matrix nrows x ncolumns
         * @param nrows
         * @param ncolumns
         * @param clearMemory if set,all elements are set to zero after creation
         */
        TMatrixEx(size_t nrows, size_t ncolumns, bool clearMemory) : TMatrix<T>(nrows, ncolumns, clearMemory) {

        }

        /**
         * Destructor
         */
        virtual ~TMatrixEx() {
            if (TMatrix<T>::data) {
                delete[] TMatrix<T>::data;
                TMatrix<T>::data = 0;
            }
        }

        /**
         *Computes the determinant of the matrix.
         *@throw MathException if this matrix is not square
         */
        T getDeterminant() throw (MathException) {
            return LinearEqCrout<T>::getDeterminant(*this);
        }

        /**
         * Returns the inverse matrix.
         * @param dest destination matrix that will stores the inverse matrix.
         */
        void getInverse(TMatrix<T> &dest) throw (MathException) {
            LinearEqCrout<T>::invertMatrix(dest, *this);
        }

        /**
         * Write a matrix to file
         * @param mat
         * @param file_name
         * @param separator_char
         * @param quot
         * @param comment
         * @param prec precision of the floating point values
         */
        static void writeToFile(const TMatrix<T>& mat,const char* file_name, const char separator_char, bool quot, const char* comment,streamsize prec) throw (MathException) {

            fstream file;

            //apre in scrittura,il file è ricreato
            file.open(file_name, fstream::out);

            const T* pdata=mat._getDataPtr();

            if (file.good()) {

                if (prec>0) file.precision(prec);

                size_t nr = mat.getRows();
                size_t nc = mat.getColumns();

                if (comment) {
                    file << '#' << comment << endl;
                }

                if (!quot) {

                    for (int i = 0; i < nr; i++) {

                        size_t offs = i*nc;

                        for (int j = 0; j < nc - 1; j++) {
                            file << pdata[offs + j] << separator_char;
                        }

                        file <<pdata[offs + nc - 1] << endl;

                    }

                } else {

                    const char quot = '"';

                    //valori separati da "
                    for (int i = 0; i < nr; i++) {

                        size_t offs = i*nc;

                        for (int j = 0; j < nc - 1; j++) {
                            file << quot << pdata[offs + j] << quot << separator_char;
                        }

                        file << quot << pdata[offs + nc - 1] << quot << endl;
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
         * This method will overwrite the destination file if it exists and uses default precision.
         *
         * @param file_name
         * @param separator_char
         * @param quot
         * @param comment
         * @param prec floating point precision (the maximum number of digits used to express floating point values)
         */
        void writeToFile(const char* file_name, const char separator_char, bool quot, const char* comment,streamsize prec) throw (MathException) {

            writeToFile(*this,file_name,separator_char,quot,comment,prec);

        }

        /**
         * Reads this matrix from a text file
         * @param separator_char the character used to separate the elements, eq. ',' for comma separated values.
         * @param quot if set,each element is enclosed in double quote (")
         */
        void readFromFile(const char* file_name, const char separator_char, bool quot) throw (MathException) {

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

                        size_t cn = TMatrix<T>::nColumns - 1;

                        if (quot) {
                            //valori fra apici ("")

                            for (size_t i = 0; i < cn; i++) {

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
                                TMatrix<T>::data[crow * TMatrix<T>::nColumns + i] = val;

                                file >> ch;

                                if (ch != separator_char) {
                                    throw MathException("Error reading value from file %s at line %d at column %d:invalid separator", file_name, line_cnt, i);
                                }
                            }

                        } else {

                            //valori senza apici
                            for (size_t i = 0; i < cn; i++) {

                                if (!(file >> val)) {

                                    throw MathException("Error reading value from file %s at line %d at column %d:invalid format", file_name, line_cnt, i);
                                }

                                //setta il valore letto nella posizione corrente
                                TMatrix<T>::data[crow * TMatrix<T>::nColumns + i] = val;

                                file >> ch;

                                if (ch != separator_char) {
                                    throw MathException("Error reading value from file %s at line %d at column %d:invalid separator", file_name, line_cnt, i);
                                }
                            }

                        }

                        int i = TMatrix<T>::nColumns - 1;
                        //ultimo valore della riga

                        if (quot) {

                            file>>ch;
                            if (!(file >> val)) {
                                throw MathException("Error reading value from file %s at line %d at column %d:invalid format", file_name, line_cnt, i);
                            }
                            file>>ch;

                        } else {
                            if (!(file >> val)) {
                                throw MathException("Error reading value from file %s at line %d at column %d:invalid format", file_name, line_cnt, i);
                            }
                        }
                        //legge l'ultimo valore della riga
                        TMatrix<T>::data[crow * TMatrix<T>::nColumns + i] = val;

                        ch = file.get();

                        while (ch != '\n' && !file.eof()) {
                            ch = file.get();
                        }

                        crow++;
                        line_cnt++;

                        if (crow == TMatrix<T>::nRows) break;

                    }
                }

                if (crow != TMatrix<T>::nRows) {
                    throw MathException("Error reading matrix from file %s:invalid row count, expected %d actual %d", file_name, TMatrix<T>::nRows, crow);
                }

                file.close();

            } else {

                file.close();
                throw MathException("Error reading the %d x %d matrix from file %s:file open or creation error.", TMatrix<T>::getRows(), TMatrix<T>::getColumns(), file_name);

            }
        }

    };

};

#endif	/* TMATRIXEX_H */

