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
 * \file tvector.h
 * \brief Generic template-base vector class.

 */

#ifndef TVECTOR_H_INCLUDED
#define TVECTOR_H_INCLUDED

#include "math/mathdefs.h"
#include "libmath.h"
#include <stddef.h>
#include <iostream>
#include <fstream>
#include <sstream>

namespace mathengine {

    /**
     * Template based n size vector.
     */
    template <class T> class TVector {
    public:

        /**
         * Creates a marix with nrows rows and ncolumns columns
         * @param nrows
         * @param ncolumns
         */
        TVector(size_t size) : data(0), nsize(0) {
            resize(size, true);
        }

        TVector(size_t size,bool clear_memory) : data(0), nsize(0) {
            resize(size, clear_memory);
        }

        virtual ~TVector() {
            if (data) {
                delete[] data;
                data = 0;
            }
        }

        /**
         * Assignment operator.
         * If v ha the same size of this vector,elements of v are copied to this vector
         * otherwise,this vector is reallocated to fit the size of v and elements are copied from v.
         * @param m
         * @return
         */
        TVector<T> & operator=(const TVector<T> &v) {

            if (v.nsize == nsize) {

                if (&v == this) return *this;

                memcpy(data, v.data, sizeof (T) * nsize);

            } else {

                resize(v.size(), false);
                memcpy(data, v.data, sizeof (T) * nsize);

            }

            return *this;
        }

        /**
         * Copy the elements of this vector into dest vector dest.
         * If destination vector has a differente size from this vector,the copy takes place only for
         * some elements.
         * @param dest
         */
        void copyTo(TVector<T>& dest) const {

            int sz=dest.size();

            if (sz>size()) sz=size();

            if (sz==0) return;

            memcpy(dest.data,data,sizeof(T)*sz);
        }

        /**
         * Sets all elements to v
         * @param v value
         * @return this vector
         */
        TVector<T> & operator=(const T &v) {

            for (size_t i = 0; i < nsize; i++) {
                data[i] = v;
            }

            return *this;
        }

        /**
         * Clear the vector and reallocates the memory.
         * Use this method to change the vector size
         * @param nrows new number of rows
         * @param ncolumns new number of columns
         * @param zero_element initial value of all elements.
         * @param resetAllelements if true,all elements are set to zero_element
         */
        inline void resize(size_t new_size, bool resetAllElements = true) {

            if (data) {
                delete[] data;
                data = 0;
            }
            nsize = new_size;
            data = new T[nsize];

            if (resetAllElements) {
                zeroMemory();
            }
        }

        /**
         * Size of this vector.
         * @return the actual number of elements in this vector.
         */
        inline size_t size() const throw () {
            return nsize;
        }

        /**
         * This operator is used to access the components of the vector.
         * eg. v.[index].The first element of the vector has index 0.
         * \warning this method doesn't perform bounds check form optimization reasons.
         * @param index index of the element
         * @return the element with the specified index.
         */
        inline T operator [] (size_t index) const {
            return data[index];
        }

        /**
         * Returns the reference to an element.
         * @param index
         * @return 
         */
        inline T & operator[] (size_t index) {
            return data[index];
        }

        /**
         * Subscript operator, can be used to access the vector elements.
         * This method can be used to modify elements in the vector.
         * eq. T &p=vector(i,j);
         * p->someMethod(...)
         * \warning this method doesn't perform bounds check form optimization reasons.
         * @param index index of the element
         * 
         * @return the reference of the element with index index
         *
         */
        inline T & operator ()(size_t index) throw () {
            return data[index];
        }

        /**
         * Subscript operator (const version), can be used to access the vector elements
         * \warning this method doesn't perform bounds check form optimization reasons.
         * @param index index of the element
         * @return the copy of the element at index.
         *
         */
        inline T operator()(size_t index) const throw () {
            return data[index];
        }

        /**
         * Get the reference of an element with indices bound check.
         * @param index index of the element (firs element has index zero)
         * 
         * @return the reference of the element at the specified index.
         * @throw MathException index is out of bounds.
         */
        T& getElementAt(size_t index) throw (MathException) {

            if (index >= nsize) {
                MathException ex("Error accessing vector element:index %d out of bounds.", index);
                throw ex;
            }

            return data[index];
        }

        void setElementAt(T value, size_t index) throw (MathException) {

            if (index >= nsize) {
                MathException ex("Error accessing vector element:index %d out of bounds.", index);
                throw ex;
            }

            data[index] = value;


        }

        /**
         * Checks if this vector has the same rows and columns as v
         * @param v a vector
         */
        inline bool hasSameSize(const TVector<T> &v) throw () {
            return nsize == v.size();
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

        //--------------- addition --------------------------------------------

        /**
         * Add vector to this vector
         * @param vector another vector
         * @thorw MathException if vectors have different sizes.
         */
        inline void add(const TVector<T> &vector) throw (MathException) {

            if (nsize != vector.size()) throw MathException("Error adding vectors:vectors have different sizes.");

            const T* pdata = vector.data;

            for (size_t i = 0; i < nsize; i++) {
                data[i] += pdata[i];
            }
        }

        /**
         * Addition operator
         *
         * @param v a vector with the same number of rows and coloumns of this vector
         * @throw MathException if matrices have different size.
         * @return the sum of this vector and the vector v
         */
        inline TVector<T> operator +(const TVector<T> &v) const throw (MathException) {

            if (nsize != v.size()) throw MathException("Error adding vectors:vectors have different sizes.");

            TVector<T> t;
            t.resize(nsize, false);

            t = *this;

            t.add(v);

            return t;

        }

        /**
         * Addition and assignment operator.
         * Add to this vector another vector.
         * @param vector a vector with the same number of rows and coloumns of this vector
         * @throw MathException if matrices have different size.
         */
        inline TVector & operator +=(const TVector<T> &vector) throw (MathException) {
            add(vector);
            return *this;
        }

        //---------------- subtraction ---------------------------------------

        /**
         * Subtracts vector from this vector.
         *
         * @param vector a vector with the same size of this vector.
         * @throw MathException if matrices have different size.
         */
        inline void sub(const TVector<T> &vector) throw (MathException) {

            if (nsize != vector.size()) throw MathException("Error subtracting vectors:vectors have different sizes.");

            const T* pdata = vector.data;

            for (size_t i = 0; i < nsize; i++) {
                data[i] -= pdata[i];
            }
        }

        /**
         * Subtraction operator.
         *
         * @param v a vector with the same size of this vector
         * @return the subtraction between this vector and  v
         */
        inline TVector operator -(const TVector<T> &v) const throw (MathException) {

            if (nsize != v.size()) throw MathException("Error subtracting vectors:invalid vectors size.");

            TVector<T> t;

            t.resize(nsize, false);

            t = *this;

            t.sub(v);

            return t;
        }

        /**
         * Subtraction and assignment operator.
         * @param vector
         * @return this vector
         */
        inline TVector & operator -=(const TVector<T> &vector) throw (MathException) {
            sub(vector);
            return *this;
        }

        //multiplication

        /**
         * Scalar product.
         * @param v a vector with the same size of this vector.
         * \waring for speed optimizations reasons this method doesn't perform bounds check,
         * make sure vectors have the same size.
         */
        inline T scalar_prod(const TVector<T> &v) const throw () {

            T r = 0;

            T* ptr = v.data;

            for (size_t i = 0; i < nsize; i++) {
                r = r + ptr[i] * data[i];
            }

            return r;
        }

        /**
         * Scalar product.
         * @param v a vector with the same size of this vector.         
         * @throw MathException if the size of v is different from this vector's size.
         */
        T operator *(const TVector<T> &v) throw (MathException) {

            if (v.size() != nsize) throw MathException("Scalar product error:vector must have size %d", nsize);

            T r = 0;

            T* ptr = v.data;

            for (size_t i = 0; i < nsize; i++) {
                r = r + ptr[i] * data[i];
            }

            return r;
        }

        /**
         * Tests if this vector equals another vector.
         * @param m
         * @return true if the vectors have the same size and each element equals the corresponding element.
         */
        bool operator ==(const TVector<T> &v) throw () {

            if (v.size() != nsize) return false;

            const T &pdata = v.data;

            for (size_t i = 0; i < nsize; i++) {
                if (pdata[i] != data[i]) return false;
            }

            return true;
        }

        /**
         * Tests if this vector is not equal to another vector.
         * @param v
         * @return
         */
        bool operator !=(const TVector<T> &v) throw () {
            return !(v == *this);
        }

        /**
         * Unary minus operator.
         * Negate all elements of the vector
         * @return this vector.
         */
        TVector<T> & operator -() throw () {
            for (size_t i = 0; i < nsize; i++) data[i] = -data[i];
            return *this;
        }

        /**
         * Internal data pointer.
         * \warning don't use this method to access the elements.
         * @return the pointer to the internal vector data.
         */
        T* _getDataPtr() const {
            return data;
        }

        void writeToFile(const char* file_name, const char separator_char, bool quot, const char* comment,streamsize prec) throw (MathException) {
            writeToFile(*this, file_name, separator_char, quot, comment,prec);
        }

        /**
         * Write a vector to file
         * @param vec the vector to write
         * @param file_name the output file name (fully qualified path)
         * @param separator_char separator character for components
         * @param quot if set,all values are quoted
         * @param comment comment added to the file header
         */
        static void writeToFile(const TVector<T>& vec, const char* file_name, const char separator_char, bool quot, const char* comment,streamsize prec) throw (MathException) {

            fstream file;

            //apre in scrittura,il file è ricreato
            file.open(file_name, fstream::out);

            if (file.good()) {

                if (prec>0) file.precision(prec);
                
                size_t nv = vec.size();

                if (comment) {
                    file << '#' << comment << endl;
                }

                if (!quot) {

                    for (int i = 0; i < nv - 1; i++) {
                        file << vec[i] << separator_char;
                    }

                    file << vec[nv - 1] << endl;

                } else {

                    const char quot = '"';

                    for (int i = 0; i < nv - 1; i++) {
                        file << quot << vec[i] << quot << separator_char;
                    }

                    file << quot << vec[nv - 1] << quot << endl;

                    file.close();
                }

            } else {
                file.close();
                throw MathException("Error writing vector to file %s:file open or creation error.", file_name);
            }
        }


    protected:

        TVector() : data(0) {

        }

        /**
         *Number of elements in the vector
         */
        size_t nsize;

        /**
         * TVector components (row major format)
         * Data store the elements of the vector arranged
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
    template<class T> ostream & operator <<(ostream& ostrm, const TVector<T> &v) {

        size_t s = v.size() - 1;

        ostrm << '(';

        for (size_t i = 0; i < s; i++) {

            ostrm << v(i) << ',';
        }

        ostrm << v(v.size() - 1) << ')';

        return ostrm;
    }

}; //end namespace

#endif

