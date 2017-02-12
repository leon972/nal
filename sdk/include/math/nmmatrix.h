/***************************************
Matrice nxm
Code by L.Berti (c) 2009
 ****************************************/
/*!
 * \file nmmatrix.h
 * \brief NxM (N rows M columns) matrix class.
 
 */

#ifndef NMMATRIX_H_INCLUDED
#define NMMATRIX_H_INCLUDED

#include "math/mathdefs.h"
#include "libmath.h"
#include "math/nvector.h"
#include <stddef.h>
#include <iostream>

/**
 *Identity value for LoadIdentity
 */
#define IDENTITY_VALUE 1.0

namespace mathengine {

    typedef enum {
        ABOVE_DIAGONAL = 0, //elementi sopra la diagonale
        UNDER_DIAGONAL = 1 //elementi sotto la diagonale

    } MATRIX_ELEMS;

    using namespace std;

    /**
     *Generic rectangular matrix of REAL.
     *
     */
    class MATH_EXPORT NMMatrix {
    protected:

        /**
         *Number of rows.
         */
        int n;
        /**
         *Number of columns.
         */
        int m;

        /**
         *Size of the data used to store the matrix.
         */
        size_t matrix_size;

        /**
         * Matrix components (row major format)
         * Data store the elements of the matrix arranged
         * by rows 
         */
        REAL* data;

        NMMatrix();

    public:

        /**
         * Create a matrix with nrows and ncolumns.
         * The indices of rowsand columns starts at zero.
         * To access the first component of the diagonal
         *,for example, use get_value(0,0)
         * @param nrows number of rows.
         * @param ncolumns number of columns.
         * @param set_to_zero if if set, all components are set to zero on creation.
         */
        NMMatrix(int nrows, int ncolumns, bool set_to_zero = true);
        virtual ~NMMatrix();

        /**
         * Assignment operator.
         * @param m
         * @return
         */
        NMMatrix & operator =(const NMMatrix& m);

        /**
         * Recreate the matrix changing the number of rows an columns.
         * All data contained are erase.
         * @param nrows number of rows.
         * @param ncolumns number of columns.
         * @param set_to_zero if is set all components are set to zero.
         */
        void Reset(int nrows, int ncolumns, bool set_to_zero);

        /**
         * Get a component.
         * \warning this method doesn't perform bounds check form optimization reasons.
         * @param row index of row (0..N-1)
         * @param col index of column (0..M-1)
         * @return  the component.
         */
        REAL get_value(int row, int col) const;

        /**
         * Set a component.
         * \warning this method doesn't perform bounds check form optimization reasons.
         * @param row index of the row (0..N-1)
         * @param col index of the column (0..M-1)
         * @param value the value.
         */
        void set_value(int row, int col, REAL value);

        /**
         * Set the values of a row.
         * The values are copied to the location of the row.
         * \warning this method doesn't perform bounds check form optimization reasons.
         * @param row index of row (0..N-1)
         * @param values the address of an array of at least M values.Must not be null.
         */
        void set_row(int row, const REAL* values);

        /**
         * Set the the values at the beginning of a row.
         * @param row index of the row (0..N-1)
         * @param values an array with row_size elements at least.
         * @param row_size the number of elements to be set.
         */
        void set_row(int row, const REAL* values, int row_size) throw (MathException);

        /**
         * Get a row
         * @param row index of the row (0..N-1)
         * @return the read-only pointer to the row.
         */
        const REAL* get_row(int row) const;

        /**
         * Set the values of a column
         * \warning this method doesn't perform bounds check form optimization reasons.
         * @param col index of the column (0..M-1)
         * @param values and array of at least N values.(Must not be null!)
         */
        void set_column(int col, const REAL* values);

        /**
         * Set the values of a column using a vector.
         * @param col index of the column (0..M-1)
         * @param col_data a vector with the column values.
         */
        void set_column(int col, const NVector& col_data);

        /**
         * Swaps two rows
         * @param row1 index of the 1st row (0..N-1)
         * @param row2 index of the 2nd row (0..N-1)
         */
        void SwapRows(int row1, int row2);

        /**
         * Swaps two columns.
         * @param col1 index of the 1st column (0..M-1)
         * @param col2 index of the 2st column (0..M-1)
         */
        void SwapColumns(int col1, int col2);

        /**Restituisce una riga
        //nota � possibile usando questo operatore accedere agli elementi in questo modo m[i][j]
         */
        /**
         * This operator is used to access the components of the matrix.
         * eg. m.[idxrow][idxcolumn].
         * \warning this method doesn't perform bounds check form optimization reasons.
         * @param row index of the row (0..n)
         * @return the pointer to the row.
         */
        REAL * operator [] (size_t row);

        /**
         * Number of rows.
         * @return the number of rows.
         */
        int GetRows() const;

        /**
         * Number of columns.
         * @return the number of columns.
         */
        int GetColumns() const;
        /**
         * Get the size of the matrix components if bytes.
         * @return the memory in bytes used to store the elements.
         */
        size_t GetDataSize() const;

        /**
         * Matrix multiplication.
         * Multiply this matrix by another matrix.
         * The multipliation can only be done only if this matrix is NxM and m is MxQ
         * the result matrix will be NxQ
         *
         * @param m a rectangular matrix whose number of rows equals the number of columns
         * of this matrix.
         * @param result the result of this matrix multiplied by m.Must have the same number of rows of this matrix.
         */
        void mul(NMMatrix* result, const NMMatrix& m) const;

        /**
         * Matrix addition.
         * Adds this matrix to a matrix of the same size.
         * @param result the sum of this matrix to m, result must have the same size of this matrix.
         * @param m
         */
        void add(NMMatrix* result, const NMMatrix& m) const;

        /**
         * Matrix subctraction.
         * Subctracts a matrix m from this matrix.
         * @param result the result matrix, must have the same size of this matrix.
         * @param m a matrix with the same size of this matrix.
         */
        void sub(NMMatrix* result, const NMMatrix& m) const;

        /**
         * Multiply this matrix by a vector.
         * @param res a result vector whose size equals the number of rows of this matrix.
         * @param v a vector whose size equals the number of columns of this matrix.
         * @throws MathException if the vector ha a wrong size
         */
        void mul(NVector* res, const NVector& v) const throw (MathException);


        /**
         * Matrix multiplication.
         * Multiply this matrix by another matrix.
         * The multipliation can only be done only if this matrix is NxM and m is MxQ
         * the result matrix will be NxQ
         *
         * @param m a rectangular matrix whose number of rows equals the number of columns
         * of this matrix.
         * @return the result of this matrix multiplied by m.Must have the same number of rows of this matrix.
         */
        NMMatrix operator *(const NMMatrix& m) const;

        /**
         * Multiply this matrix by a vector.
         * @return a vector whose size equals the number of rows of this matrix.
         * @param v a vector whose size equals the number of columns of this matrix.
         * @throws MathException if the vector ha a wrong size
         */
        NVector operator *(const NVector& v) const throw (MathException);

        /**
         * Matrix subctraction.
         * Subctracts a matrix m from this matrix.
         * @return the result matrix.
         * @param m a matrix with the same size of this matrix.
         */
        NMMatrix operator -(const NMMatrix& m) const;

        /**
         * Logical equal operator
         * @param m
         * @return true if m equals this matrix
         */
        bool operator ==(const NMMatrix& m);

        /**
         * Inequality operator
         * @param m
         * @return true if m is different from this matrix
         */
        bool operator !=(const NMMatrix& m);

        /**
         * Matrix by scalar
         * @param a
         * @return
         */
        NMMatrix operator *(const REAL a) const;

        /**
         * Matrix by scalar and assignment operator
         * @param a
         * @return
         */
        NMMatrix & operator*=(const REAL a);

        /**
         * Set all elements to zero.
         */
        void ZeroMatrix();

    };

    /* nota:
       � possibile evitare di copiare il valore di ritorno dall'operatore
       definendo l'operatore come
       NMMatrix& operator + (const NMMatrix& m1,const NMMatrix& m2);
       restituendo cio� un reference.In questo modo pero' non si puo' restituire
       un reference ad un oggetto locale.Occorre definire uno stack di matrici
       temporanee e restituire un elemento di quello stack.Questo metodo pero' non funziona
       se i calcoli coinvolgono piu' matrici temporanee di quanti sono gli elementi dello stack
     */

    //NMMatrix operator + (const NMMatrix& m1,const NMMatrix& m2);

    /**
     *A square matrix.
     *This matrix has the same number of rows and columns.
     */
    class MATH_EXPORT Matrix : public NMMatrix {
    protected:

    public:

        /**
         * Constructs a square matrix with n rows and n columns.
         * @param n matrix size
         * @param set_to_zero if is set, all elements are set to zero.
         */
        Matrix(int n, bool set_to_zero = true);

        /**
         * Assignment operator
         * copy elements from matrix m to this marix.
         * @param m
         * @return 
         */
        Matrix& operator = (const Matrix& m);

        /**
         * Memory deallocation.
         */
        virtual ~Matrix();

        // Matrix& Matrix::operator = (const Matrix& mat);

        /**
         * Transform this matrix to the  identity matrix.
         */
        void LoadIdentity();

        /**
         * Get the diagonal of this matrix.
         * \warning this method doesn't perform bounds check for optimization reasons.
         * @param diag the array already sized with n elements where the diagonal elements are copied.(Must not be null!)
         */
        void GetDiagonal(REAL* diag);

        /**
         * Matrix multiplication
         * @param m
         * @return
         */
        Matrix operator *(const NMMatrix& m) const;

        /**
         * Multiplication and assignment operator
         * @param m matrix to add
         * @return this matrix after the multiplication with m
         */
        Matrix & operator *=(const NMMatrix& m);

        /**
         * Matrix addition.
         * @param m
         * @return
         */
        Matrix operator +(const NMMatrix& m) const;

        /**
         * Addition and assignment operator
         * @param m matrix to add
         * @return this matrix after adding m
         */
        Matrix & operator +=(const NMMatrix& m);

        /**
         * Matrix subtraction.
         * @param m
         * @return
         */
        Matrix operator -(const NMMatrix& m) const;


        /**
         * Matrix subtraction and assignment operator.
         * @param m matrix to subtract
         * @return  this matrix after subtracting m.
         */
        Matrix & operator -=(const NMMatrix& m);

        /**
         * Matrix by vector multiplication.
         * @param v a vector with the same size of this matrix.
         * @return a vector with the result of the multiplication.
         */
        NVector operator *(const NVector& v) const throw (MathException);


        /**
         * Matrix by scalar multiplication
         * @param a
         * @return
         */
        Matrix operator *(const REAL a) const;

        /**
         * Matrix by scalar multiplication and assignment operator
         * @param elems
         */
        Matrix & operator *=(const REAL a);

        /**
         * Logical equal operator
         * @param m
         * @return true if m equals this matrix
         */
        bool operator ==(const NMMatrix& m);

        /**
         * Inequality operator
         * @param m
         * @return true if m is different from this matrix
         */
        bool operator !=(const NMMatrix& m);

        /**
         * Makes this matrix symmetric.
         * @param elems ABOVE_DIAGONAL to keep the elements above the diagonal,
         * UNDER_DIAGONAL to keep the elements under the diagonal.
         */
        void MakeSymmetric(MATRIX_ELEMS elems);

        /**
         * Makes this matrix triangular.
         * @param elems ABOVE_DIAGONAL to keep the elements above the diagonal,
         * UNDER_DIAGONAL to keep the elements under the diagonal.
         */
        void MakeTriangular(MATRIX_ELEMS elems);


        /**
         * Makes a band matrix.
         * @param width size of the band.
         */
        void MakeDiagonal(int width);

    };

    /**
     * Print a matrix to output stream
     * @param s
     * @param m
     * @return
     */

    ostream & operator<<(ostream& s, const NMMatrix& m);
};

#endif // NMMATRIX_H_INCLUDED
