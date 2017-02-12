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
 \file ilinear.h
 * \brief Interface for a linear equations set
 */


#ifndef ILINEAR_H_INCLUDED
#define ILINEAR_H_INCLUDED

#include "mathdefs.h"
#include "libmath.h"
#include "mexception.h"
#include "tvector.h"
#include "tmatrix.h"

/**
 \brief math engine objects
 */
namespace mathengine {

    class Matrix;
    class MathException;

    /**
     *Solution optimizations
     */
    typedef enum {
        OPT_NONE = 0, /*!<- no optimization*/
        OPT_ROUNDOFFS = 1, /*!<- minimize roundoff errors*/
        OPT_DONTARRANGE_MATRIX = 2 /*!<- do not arange equations before solving*/

    } OPTIMIZE_FLAGS;

    /**
     * Interface for linear algebraic set of equations
     */
    template<class T> class LinearEquationsSet {
    public:
        
        /**
         * Destructor
         * @param eq_index
         * @param adjcoeff
         * @return
         */
        virtual ~LinearEquationsSet() {
        };

        /**
         * Set the coefficients and right-hand member for an equation of the set.
         * @param eq_index index of equation (index of row in the matrix of coefficients)
         * @param coeff vector of coefficients (known members).Must have size EquationsCount
         * @param rhs right-hand member
         * @return the equation's index eq_index
         */
        virtual int setEquation(int eq_index, const T* coeff, T rhs) throw (MathException) = 0;


        /**
         * Set the row of coefficients in the augmentmat(eq_num,eq_num+1,true),res(eq_num)ed matrix (known coefficients + right hand member)
         * @param eq_index index of row
         * @param adjcoeff coefficients + right hand member
         * @return  the index of row
         */
        virtual int setEquation(int eq_index, const T* adjcoeff) throw (MathException) = 0;

        /**
         * Set the coefficient matrix and right-hand column
         * @param coeff coefficient matrix
         * @param rhs right hand column (must have the same size as coeff)
         */
        virtual void setEquations(const TMatrix<T>& coeff, const TVector<T>& rhs) throw (MathException) = 0;

        /**
         * Get the augmented matrix of coefficients+known values vector
         * @return a matrix of n rows and n+1 columns, where n is the number of equations.
         */
        virtual const TMatrix<T>& getAugmentedMatrix() const = 0;


        /**
         * Copy the actual coefficients matrix to a destination matrix.
         * \warning Calling this method after solving may copy the trasformed/reduced matrix
         * and not the original coefficients matrix.
         * @param dest_matrix
         */
        virtual void getCoefficientsMatrix(TMatrix<T>* dest_matrix) throw (MathException)=0;
        

        /**
         * Set the right-hand members column (known members)
         * @param rhs
         */
        virtual void setRightHandColumn(const TVector<T>& rhs) throw (MathException) = 0;

        /**
         * Sets a coefficient in the coefficients matrix
         * @param eq_index equation index (index of row of the matrix)
         * @param column column index
         * @param coeff value of the coefficient
         * \warning this method doesn't perform bounds check
         */
        virtual void setCoefficient(T coeff, int eq_index, int column) throw () = 0;

        /**
         * Gets a coefficient from the coefficients matrix
         * @param eq_index equation index (index of row of the matrix)
         * @param column column index
         * @return the value of the coefficient
         * \warning this method doesn't perform bounds check
         */
        virtual T getCoefficient(int eq_index, int column) const throw () = 0;

        /**
         * Sets a value in the right hand knwown terms.
         * @param rhs value
         * @param eq_index
         * \warning this method doesn't perform bounds check
         */
        virtual void setRhsValue(T rhs, int eq_index) throw () = 0;

        /**
         * Get a value from the column of the known terms
         * @param eq_index index or equations
         * @return the value of the known term
         */
        virtual T getRhsValue(int eq_index) const throw () = 0;

        /**
         * Get the column of known terms.
         * The column of known terms is copied to rhs_column
         * @return
         */
        virtual void getRhsColumn(TVector<T>* rhs_column) const throw (MathException) = 0;

        /**
         * Get the number of equations in the set
         * @return
         */
        virtual int equationCount() const = 0;

        /**
         * Get the rank of the coefficients matrix.
         * @return
         */
        virtual int getRank() const throw (MathException) = 0;

        /**
         * Solve the equations.
         * @throw MathException if the equations can not be solved.
         */
        virtual void solve() throw (MathException) = 0;

        /**
         * Solve the equations.
         * @param opt optimizations flags (OPT_NONE,OPT_ROUNDOFFS=1,OPT_DONTARRANGE_MATRIX)
         */
        virtual void solve(int opt) throw (MathException) = 0;

        /**
         * Return true if the set has been solved (after calling Solve method).
         * @return
         */
        virtual bool isSolved() const = 0;

        /**
         * Get the computed unknown x value
         * @param eq_index index of result values
         * @return
         */
        virtual T getX(int eq_index) const throw (MathException) = 0;

        /**
         * Get the results vector.If the equations set is not solved throw and exception.
         *
         * @param result the vector where result is copied
         */
        virtual void getResult(T* result) const throw (MathException) = 0;

        /**
         * Copy the results vector in res
         * @param res
         */
        virtual void getResult(TVector<T>* res) const throw (MathException) = 0;
        

        /**
         * Verify the result.
         *
         * @param calc_knowns calculated known terms.cal_knowns is obtained by multiplying
         * the coefficients matrix with by the result vector.
         * @param abs_errors difference between calculated known terms and actual known terms
         * @throw MathException if the system is not solved.        
         */
        virtual void verifyResult(TVector<T>& calc_knowns,TVector<T>* abs_error,const TMatrix<T>& aij,const TVector<T>* rhs) const throw (MathException) =0;

        /**
         *Method used for solving the system
         */
        virtual const char* getMethodName() const=0;

    };


};


#endif // ILINEAR_H_INCLUDED
