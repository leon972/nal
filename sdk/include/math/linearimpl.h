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
 \file linearimpl.h
 \brief Base class for linear equations set.  
 */

#ifndef LINEARIMPL_H_INCLUDED
#define LINEARIMPL_H_INCLUDED

#include "ilinear.h"
#include "tvector.h"
#include "tmatrix.h"
#include "libmath.h"


namespace mathengine {

    /**
     * \brief Base class for a set of linear algebraic equations.
     */
    template<class T>class MATH_EXPORT LinearImpl : public virtual LinearEquationsSet<T> {
    protected:

        bool solved;

        TMatrix<T> mat; //matrice dei coefficienti + colonna dei termini noti (matrice aumentata)
        //nota bene se n � il numero di eq. la matrice mat
        //ha n righe e n+1 colonne (comprende la colonna dei termini noti)
        TVector<T> res; //vettore dei risultati
        int rank;
        int n;
        int *eqmap; //indici delle equazioni originali

        /**
         * Swap equations with index e1 and eq2
         * @param eq1
         * @param eq2
         */
        inline void swapEquations(int eq1, int eq2) {
           mat.swapRows(eq1, eq2);
           int temp = eqmap[eq2];
           eqmap[eq2] = eqmap[eq1];
           eqmap[eq1] = temp;
        }

        inline void resetEqMap() {
            for (int i = 0; i < n; i++) eqmap[i] = i;
        }        

    public:

        /**
         * Create a set of linear equations with eq_num equations.
         * @param eq_num numer of equations in the set.
         */
        LinearImpl(int eq_num) : mat(eq_num, eq_num + 1), res(eq_num) {
            solved = false;
            rank = 0;
            n = eq_num;
            eqmap = new int[eq_num];
            resetEqMap();
        }

        /**
         * Destructor
         */
        virtual ~LinearImpl() {

            if (eqmap) {
                delete[] eqmap;
                eqmap = 0;
            }
        }

        /**
         * Set the coefficients and right-hand member for an equation of the set.
         * @param eq_index index of equation (index of row in the matrix of coefficients)
         * @param coeff vector of coefficients (known members).Must have size EquationsCount
         * @param r known term of the equation.
         * @return the equation's index eq_index
         */
        virtual int setEquation(int eq_index, const T* coeff, T r) throw (MathException) {

            if (eq_index < 0 || eq_index >= n) throw MathException("equation index out of bounds:the index must be between %d and % , actual value %d", 0, n - 1, eq_index);

            if (!coeff) throw MathException("invalid coefficients:null not allowed");

            for (int i = 0; i < n; ++i) mat[eq_index][i] = coeff[i];

            //imposta il termine noto
            mat[eq_index][n] = r;

            solved = false;

            return eq_index;

        }

        /**
         * Set the row of coefficients in the augmented matrix (known coefficients + right hand member)
         * This method copy the content of adjcoeff into the internal data array.
         * @param eq_index index of row
         * @param adjcoeff coefficients + right hand member
         * @return  the index of row
         */
        int setEquation(int eq_index, const T* adjcoeff) throw (MathException) {

            if (!adjcoeff) throw MathException("invalid coefficients:null not allowed");

            mat.setRow(eqmap[eq_index], adjcoeff);

            solved = false;

            return eq_index;

        }

        /**
         * Get the augmented marix of coefficients+known values vector
         * (const version)
         * \warning this method returns the internal matrix that can be modified by the solving method
         * Is not safe to call this method after solving the equations
         * @return a matrix of n rows and n+1 columns, where n is the number of equations.
         */
        const TMatrix<T>& getAugmentedMatrix() const {
            return mat;
        }

        /**
         * Sets a coefficient in the coefficients matrix
         * @param eq_index equation index (index of row of the matrix)
         * @param column column index
         * @param coeff value of the coefficient
         * \warning this method doesn't perform bounds check
         */
        inline void setCoefficient(T coeff, int eq_index, int column) throw () {

            mat[eqmap[eq_index]][column] = coeff;
            solved = false;
        }

        /**
         * Gets a coefficient from the coefficients matrix
         * @param eq_index equation index (index of row of the matrix)
         * @param column column index
         * @return the value of the coefficient
         * \warning this method doesn't perform bounds check
         */
        T getCoefficient(int eq_index, int column) const throw () {
            return mat[eqmap[eq_index]][column];
        }

        /**
         * Sets a value in the right hand knwown terms.
         * @param rhs value
         * @param eq_index
         * \warning this method doesn't perform bounds check
         */
        void setRhsValue(T rhs, int eq_index) throw () {
            mat[eqmap[eq_index]][n] = rhs;
            solved = false;
        }

        /**
         * Gets a value from the column of the known terms
         * @param eq_index index or equations
         * @return the value of the known term
         */
        T getRhsValue(int eq_index) const throw () {
            return mat[eqmap[eq_index]][n];
        }

        /**
         * Get the column of known terms.
         * The column of known terms is copied to rhs_column
         * @return
         */
        void getRhsColumn(TVector<T>* rhs_column) const throw (MathException) {

            if (rhs_column == 0) return;

            if (rhs_column->size() != n) {
                throw MathException("Error acquiring the column of known terms:invalid destination vector's size");
            }

            for (int i = 0; i < n; i++) {
                rhs_column->operator[](eqmap[i]) = mat[i][n];
            }
        }

        /**
         * Set the coefficient matrix and right-hand column
         * \warning This method copies coeff and rhs into the internal matrix.
         * @param coeff coefficient matrix
         * @param rhs right hand column (must have the same size as coeff)
         */
        virtual void setEquations(const TMatrix<T>& coeff, const TVector<T>& rhs) throw (MathException) {

            if (coeff.getRows() != n || coeff.getColumns() != n) throw MathException("Invalid coefficients matrix,the matrix must be %d x %d, (actual %d x %d) ", n, n, coeff.getRows(), coeff.getColumns());
            if (rhs.size() != static_cast<size_t> (n)) throw MathException("Invalid rhs column");
            //copia la memoria da coeff a mat

            for (size_t r = 0; r < n; r++) {
                mat.setRow(r, coeff._getRowPtr(r), n);
            }

            //coefficients
            mat.setColumn(n, rhs._getDataPtr());

            resetEqMap();

            solved = false;

        }
        

        /**
         * Set the right-hand coefficient column.
         * With the Crout's reduction method can be found the solution for a set of known terms
         * without recomputing the auxiliary matrix.
         * @param rhs
         */
        virtual void setRightHandColumn(const TVector<T>& rhs) throw (MathException) {

            if (rhs.size() != static_cast<size_t> (LinearImpl<T>::n)) throw MathException("Invalid rhs column expected size %d , actual %d", LinearImpl<T>::n, rhs.size());

            //quando si imposta questa colonna si deve controllare se nel calcolo della
            //matrice ausiliaria sono state scambiate delle equazioni
            for (int i = 0; i < LinearImpl<T>::n; ++i) {
                LinearImpl<T>::mat[i][LinearImpl<T>::n] = rhs[eqmap[i]];
            }

            LinearImpl<T>::solved = false;
        }

        /**
         * Get the number of equations in the set
         * @return
         */
        virtual int equationCount() const {
            return n;
        }

        /**
         * Get the rank of the coefficients matrix.
         * @return
         */
        virtual int getRank() const throw (MathException) {
            if (!solved) throw MathException("cannot evaluate rank:system not solved.");
            return rank;
        }

        /**
         * Return true if the set has been solved (after calling Solve method).
         * @return
         */
        virtual bool isSolved() const {
            return solved;
        }

        /**
         * Get the computed x value
         * @param eq_index index of result values
         * @return
         */
        virtual T getX(int eq_index) const throw (MathException) {

            if (eq_index < 0 || eq_index >= n) throw MathException("cannot get the x%d result:index %d out of bounds (expected between 0 and %d)", eq_index, eq_index, n - 1);

            if (!solved) throw MathException("cannot get x%d :equations not solved", eq_index);

            return res[eq_index];
        }

        /**
         * Get the results vector.If the equations set is not solved throw and exception.
         *
         * @param result the vector where result is copied
         */
        virtual void getResult(T* result) const throw (MathException) {

            if (!solved) throw MathException("equations not solved");
            if (!result) throw MathException("result vector is null");

            for (int i=0;i<n;i++) {
                result[i]=res[i];
            }
        }

        /**
         * Copy the results vector in res
         * @param res
         */
        virtual void getResult(TVector<T>* dest_res) const throw (MathException) {

            if (!solved) throw MathException("Error acquiring equations results:equations not solved yet");
            if (dest_res == 0) throw MathException("Error acquiring equations results:destination vector is null");
            if (dest_res->size() != n) throw MathException("Error acquiring equations results:destinazione vector must have size %d (actual %d)", n, dest_res->size());

            //nota bene:scambiare le equazioni fra di loro (e i termini noti corrispondenti)
            //non altera il vettore della soluzione quindi non è necessario usare il mapping per ottenere gli x[i]

            for (int i=0;i<n;i++) {
                dest_res->setElementAt(res[i],i);
            }
        }
        

        /**
         * Get the augmented matrix.
         * \warning returns a reference to an internal matrix.
         * @return the matrix of coefficients + right hand members. (nxn+1 matrix)
         */
        TMatrix<T>& getAugmentedMatrix() throw (MathException) {
            return LinearImpl<T>::mat;
        }

        /**
         * Copy the actual coefficients matrix to a destination matrix.
         * \warning Calling this method after solving may copy the trasformed/reduced matrix
         * and not the original coefficients matrix.
         * @param dest_matrix
         */
        void getCoefficientsMatrix(TMatrix<T>* dest_matrix) throw (MathException) {

            if (!dest_matrix) return;

            if (!dest_matrix->isSquare() || dest_matrix->getRows()!=n) {
                throw MathException("Error copying coefficients matrix:destination matrix has wrong size");
            }

            mat.copySubMatrixTo(*dest_matrix,0,0,0,0,n-1,n-1,true);
        }

        /**
         * Esegue il pivoting su tutte le equazioni
         * @param opt
         */
        void pivoting(int opt) {
            pivoting(opt,0);
        }

        /**
         * Pivoting,arrange equations in order to eliminate zero elements on diagonal.
         * @param opt if OPT_ROUNDOFFS is set performs a deep pivoting (elements with max absolute values are put on diagonal)
         */
        void pivoting(int opt,int from_index) throw (MathException) {

            TMatrix<T>& mat=LinearImpl<T>::mat;

            int n=LinearImpl<T>::n;

            if ((opt & OPT_ROUNDOFFS) > 0) {

                T temp = T();

                int idx;

                //minimize roundoffs O(n^2)
                //puts max elements on diagonal
                for (int i = from_index; i < n; i++) {

                    T max = T();

                    for (int j = from_index; j < n; j++) {

                        temp = abs(mat[j][i]);

                        if (temp > max) {
                            max = temp;
                            idx = j;
                        }
                    }

                    if (idx != i) {
                        LinearImpl<T>::swapEquations(i, idx);
                    }
                }

            } else {

                //lightweight pivoting
                for (int i = 0; i < n; i++) {

                    if (abs(mat[i][i]) < EPSILON) {

                        for (int j = 0; j < n; j++) {

                            if (abs(mat[j][i]) > EPSILON && abs(mat[i][j]) > EPSILON) {

                                LinearImpl<T>::swapEquations(i, j);
                                break;

                            }
                        }
                    }
                }
            }
        }

        /**
         * Verify the result.
         *
         * @param calc_knowns calculated known terms.cal_knowns is obtained by multiplying
         * the coefficients matrix by the result vector.
         * @param abs_errors (optional,may be null) difference between calculated known terms and actual known terms
         * @throw MathException if the system is not solved.
         * @return
         */
        void verifyResult(TVector<T>& calc_knowns, TVector<T>* abs_error, const TMatrix<T>& aij, const TVector<T>* rhs) const throw (MathException) {

            if (!isSolved()) {
                throw MathException("Cannot verify result:system is not solved yet.");
            }

            if (calc_knowns.size() != n) {
                throw MathException("Error verifying result:invalid vector size.");
            }

            TVector<T> r(n);

            getResult(&r);
           
            calc_knowns = aij*r;

            if (abs_error && rhs) {

                for (int i = 0; i < n; i++) {
                    abs_error->operator[](i) = calc_knowns[i] - rhs->operator[](i);
                }
            }
        }

        /**
         *Copy indices of equations into eqmap.
         *eqmap contains indices after scrambling
         */
        void GetEquationsIndices(TVector<int>& dest_eqmap) {

            for (int i=0;i<n;i++) {
                dest_eqmap[i]=eqmap[i];
            }
        }

        /**
         *Detect if some equations are swapped after pivoting
         */
        bool hasScrambled() {
            for (int i=0;i<n;i++) {
                if (eqmap[i]!=i) return true;
            }

            return false;
        }
    };
};


#endif // LINEARIMPL_H_INCLUDED
