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
 \file tridiag.h
 \brief Numerical solution of linear tridiagonal sets of equations with Crout's reduction method.
 */

#ifndef TRIDIAG_H
#define	TRIDIAG_H

#include <iostream>
#include "math/mathdefs.h"
#include "math/ilinear.h"
#include "math/bandmatrix.h"
#include "math/tmatrix.h"

using namespace std;

namespace mathengine {

    /**
     *Tridiagonal set of linear equations.
     *This class uses the Crout's reduction method to solve the equations.
     *The reduction may be used to solve a set for many different kwnows columns.
     *
     */
    template<typename T> class TridiagEq : public LinearEquationsSet<T> {
    public:

        /**
         * Creates a set with the specified number of equations
         * @param eq_num number of equations
         */
        TridiagEq(const int eq_num) : n(eq_num), mat(eq_num, 1, 1),
        rhs(eq_num), rhs1(eq_num), res(eq_num), invalid_aux_matrix(true), solved(false), augmentedMatrix(0) {

        }

        virtual ~TridiagEq() {

            if (augmentedMatrix) {
                delete augmentedMatrix;
                augmentedMatrix = 0;
            }
        }

        /**
         * Set the coefficients and right-hand member for an equation of the set.
         * @param eq_index index of equation (index of row in the matrix of coefficients)
         * @param coeff the non-zero elements of the equations with index eq_index.Must have size 3
         * @param rhs right-hand member
         * @return the equation's index eq_index
         */
        int setEquation(int eq_index, const T* coeff, T rhs_value) throw (MathException) {

            invalid_aux_matrix = true;
            solved = false;
            mat.setElementAt(eq_index, eq_index - 1, coeff[0]);
            mat.setElementAt(eq_index, eq_index, coeff[1]);
            mat.setElementAt(eq_index, eq_index + 1, coeff[2]);
            rhs[eq_index] = rhs_value;

        }

        /**
         * Set the row of coefficients in the augmentmat(eq_num,eq_num+1,true),res(eq_num)ed matrix (known coefficients + right hand member)
         * @param eq_index index of row
         * @param adjcoeff non zero coefficients + right hand member.Must have size 4
         * @return  the index of row
         */
        int setEquation(int eq_index, const T* adjcoeff) throw (MathException) {
            setEquation(eq_index, adjcoeff, adjcoeff[3]);
        }

        /**
         * Set the coefficient matrix and right-hand column
         * @param coeff coefficient matrix
         * @param rhs right hand column (must have the same size as coeff)
         */
        void setEquations(const TMatrix<T>& coeff, const TVector<T>& rhs) throw (MathException) {

            const char* err = "Error setting equations in the tridiagonal set:";

            if (coeff.getRows() < n) throw MathException("%s:coefficients matrix must have %d rows at least.", err, n);
            if (coeff.getColumns() < n) throw MathException("%s:coefficients matrix must have %d columns at least.", err, n);

            for (int i = 0; i < n; i++) {

                mat.setElementAt(i, i - 1, coeff[i][i - 1]);
                mat.setElementAt(i, i, coeff[i][i]);
                mat.setElementAt(i, i + 1, coeff[i][i + 1]);
            }

            solved = false;
            invalid_aux_matrix = true;

        }

        /**
         * Get the augmented matrix of coefficients+known values vector
         * @return a matrix of n rows and n+1 columns, where n is the number of equations.
         */
        const TMatrix<T>& getAugmentedMatrix() const {

            if (!augmentedMatrix) {
                augmentedMatrix = new TMatrix<T > (n, n + 1);
            }

            mat.copyTo(*augmentedMatrix);
            augmentedMatrix->setColumn(n, rhs);

            return *augmentedMatrix;
        }

        /**
         * Copy the actual coefficients matrix to a destination matrix.
         * \warning Calling this method after solving may copy the trasformed/reduced matrix
         * and not the original coefficients matrix.
         * @param dest_matrix
         */
        void getCoefficientsMatrix(TMatrix<T>* dest_matrix) throw (MathException) {

            if (dest_matrix == 0 || dest_matrix->getRows() < n || dest_matrix->getColumns() < n) throw MathException("Error acquiring coefficients matrix from tridiagonal set:invalid destination matrix.");

            
            mat.copyTo(*dest_matrix);
        }

        /**
         * Return the coefficients band matrix
         * @return
         */
        BandMatrix<T>& getCoefficientsMatrix() {
            return mat;
        }

        /**
         * Set the right-hand members column (known members)
         * @param rhs
         */
        void setRightHandColumn(const TVector<T>& rhs) throw (MathException) {

            if (rhs.size() < n) throw MathException();

            rhs.copyTo(this->rhs);

            solved = false;
        }

        /**
         * Sets a coefficient in the coefficients matrix
         * @param eq_index equation index (index of row of the matrix)
         * @param column column index
         * @param coeff value of the coefficient
         * \warning this method doesn't perform bounds check
         */
        void setCoefficient(T coeff, int eq_index, int column) throw () {

            if (eq_index < 0 || eq_index >= n || column < 0 || column >= n) return;

            mat.setElementAt(eq_index, column, coeff);

            solved = false;
            invalid_aux_matrix = true;

        }

        /**
         * Gets a coefficient from the coefficients matrix
         * @param eq_index equation index (index of row of the matrix)
         * @param column column index
         * @return the value of the coefficient
         * \warning this method doesn't perform bounds check
         */
        T getCoefficient(int eq_index, int column) const throw () {

            if (eq_index < 0 || eq_index >= n || column < 0 || column >= n) return T();

            return mat.getElementAt(eq_index, column);
        }

        /**
         * Sets a value in the right hand knwown terms.
         * @param rhs value
         * @param eq_index
         * \warning this method doesn't perform bounds check
         */
        void setRhsValue(T rhs, int eq_index) throw () {

            if (eq_index >= 0 && eq_index < n) {
                this->rhs[eq_index] = rhs;
                solved = false;
            }

        }

        /**
         * Get a value from the column of the known terms
         * @param eq_index index or equations
         * @return the value of the known term
         */
        T getRhsValue(int eq_index) const throw () {

            if (eq_index >= 0 && eq_index < n) {
                return rhs[eq_index];
            }

            return T();

        }

        /**
         * Get the column of known terms.
         * The column of known terms is copied to rhs_column
         * @return
         */
        void getRhsColumn(TVector<T>* rhs_column) const throw (MathException) {

            if (!rhs_column || rhs_column->size() != n) throw MathException("Error acquiring rhs column:invalid destination vector.");

            rhs.copyTo(*rhs_column);

        }

        /**
         * Gets the reduced rhs column
         * @return
         */
        void getReducedRhsColumn(TVector<T>* rhs_reduced) const throw (MathException) {

            if (!solved) throw MathException("Cannot acquire rhs reduced column:not solved yet");

            rhs1.copyTo(*rhs_reduced);

        }

        /**
         * Get the number of equations in the set
         * @return
         */
        int equationCount() const {
            return n;
        }

        /**
         * Get the rank of the coefficients matrix.
         * @return
         */
        int getRank() const throw (MathException) {
            return n;
        }

        /**
         * Solve the equations.
         * @throw MathException if the equations can not be solved.
         */
        void solve() throw (MathException) {
            solve(0);
        }

        /**
         * Solve the equations.
         * For tridiagonal equations this method is equivalent to Solve()
         * @param opt optimizations flags, ignored for tridiagonal sets.
         */

        void solve(int opt) throw (MathException) {

            if (invalid_aux_matrix) {

                //internal matrix memory
                T* pdata = mat._getDataPtr();

                const int pitch = mat._getPitch();
                int l = mat.getLdiag();
                int h = mat.getHdiag();

                //apply the reduction to matrix mat
                //first row
                T d1 = mat(0, 0);
                T d2;

                if (d1 == 0) throw MathException("Cannot solve:diagonal element equals zero.");

                mat(0, 1) = mat(0, 1) / d1;
                rhs1[0] = rhs[0] / d1;

                size_t offs = 0;
                size_t offs_prev = l;


                for (int i = 1; i < n; i++) {

                    offs = i * pitch + l;

                    //d2=mat(i,i)-mat(i,i-1)*mat(i-1,i);
                    d2 = pdata[offs] - pdata[offs - 1] * pdata[offs_prev + 1];
                    //mat(i,i)=d2;
                    pdata[offs] = d2;
                    //mat(i,i+1)=mat(i,i+1)/d2;
                    pdata[offs + 1] /= d2;
                    //rhs1[i]=(rhs[i]-mat(i,i-1)*rhs1[i-1])/d2;
                    rhs1[i] = (rhs[i] - pdata[offs - 1] * rhs1[i - 1]) / d2;

                    offs_prev = offs;
                }

                invalid_aux_matrix = false;

            } else {
                //matrice già ridotta, deve solo aggiornare la colonna di appoggio
                //rhs1

                const int pitch = mat._getPitch();

                int l = mat.getLdiag();

                const T* pdata = mat._getDataPtr();

                const T d1 = mat(0, 0);

                T d2;

                rhs1[0] = rhs[0] / d1;

                size_t offs = 0;
                size_t offs_prev = l;

                for (int i = 0; i < n; i++) {

                    offs = i * pitch + l;

                    d2 = pdata[offs] - pdata[offs - 1] * pdata[offs_prev + 1];

                    rhs1[i] = (rhs[i] - pdata[offs - 1] * rhs1[i - 1]) / d2;

                    offs_prev = offs;
                }

            }

            res[n - 1] = rhs1[n - 1]; //xn

            for (int i = n - 2; i >= 0; i--) {
                res[i] = rhs1[i] - mat(i, i + 1) * res[i + 1];
            }

            solved = true;
        }

        /**
         * Return true if the set has been solved (after calling Solve method).
         * @return
         */
        bool isSolved() const {
            return solved;
        }

        /**
         * Get the computed unknown x value
         * @param eq_index index of result values
         * @return
         */
        T getX(int eq_index) const throw (MathException) {

            if (eq_index < 0 || eq_index >= n) throw MathException("cannot get the x%d result:index %d out of bounds (expected between 0 and %d)", eq_index, eq_index, n - 1);

            if (!solved) throw MathException("cannot get x%d :equations not solved", eq_index);

            return res[eq_index];

        }

        /**
         * Get the results vector.If the equations set is not solved throw and exception.         
         * @param result the vector where result is copied
         */
        void getResult(T* result) const throw (MathException) {

            if (!solved) throw MathException("Cannot acquire results:equations not solved");

            for (int i = 0; i < n; i++) {
                result[i] = res[i];
            }

        }

        /**
         * Copy the results vector in res
         * @param res
         */
        void getResult(TVector<T>* res) const throw (MathException) {

            if (!solved) throw MathException("Cannot acquire results:equations not solved");

            if (res == 0 || res->size() < n) throw MathException("Error acquiring results:destination bad size");

            this->res.copyTo(*res);
        }

        /**
         * Verify the result.
         *
         * @param calc_knowns calculated known terms.cal_knowns is obtained by multiplying
         * the coefficients matrix by the result vector.
         * @param abs_errors difference between calculated known terms and actual known terms
         * @throw MathException if the system is not solved.
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
         *Method used for solving the system
         */
        const char* getMethodName() const {
            return "Crout reduction method (tridiagonal set)";
        }

        /**
         * Number of non-zero elements right of the diagonal
         * @return
         */
        int getHdiag() {
            return mat.getHdiag();
        }

        /**
         *Number of non-zero elements left of the diagonal
         */
        int getLdiag() {
            return mat.getLdiag();
        }

    protected:

        //number of equations
        int n;
        //coefficients matrix (is a tridiagonal matrix)
        BandMatrix<T> mat;
        //rhs column of the reduced matrix
        TVector<T> rhs1;
        //actual column of knowns terms
        TVector<T> rhs;
        //results column
        TVector<T> res;
        //this flag is set if the reduction is needed
        bool invalid_aux_matrix;
        //this flag is set after solving the equations
        bool solved;
        //used as return value
        mutable TMatrix<T>* augmentedMatrix;

    };
};

#endif	/* TRIDIAG_H */

