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
 \file crout.h
 \brief Numerical solution of linear equations and matrix inversion with the Crout reduction method.
 */


#ifndef CROUT_H_INCLUDED
#define CROUT_H_INCLUDED

#include "libmath.h"
#include "linearimpl.h"
#include "tmatrix.h"
#include <vector>

namespace mathengine {

    using namespace std;

    /**
     * Solve a set of linear algebraic equations using Crout reduction.
     * The Crout reduction is useful to investigate the variation of the result vector
     * when the right-hand members changes.
     * The solution can be found for different values of the known column without recomputing
     * the auxiliry matrix.
     * This method may also be used to compute the inverse matrix of aij and the determinant of aij.
     */

    template<class T> class LinearEqCrout : public LinearImpl<T> {
    public:

        /**
         * Compute the Crout's reduction "in place" on the matrix a
         * The matrix a must have ncols>=nrows, the reduction is performed only on the coefficients matrix.
         * \warning this method doesn't perform pivoting, matrix a should be pivoted before the reduction
         * (all diagonal elements must be != 0)
         * @param a
         */
        static void croutReduction(TMatrix<T>& a) {

            size_t n = a.getRows();
            size_t nc=a.getColumns();
            size_t offset;

            register int i, j, k;

            T* pdata=a._getDataPtr();

            for (k = 0; k < n; ++k) {

                T akk = a[k][k];

                offset=nc*k;      //data[nColumns * row]

                for (i = k + 1; i < n; i++) {
                    //a[k][i] /= akk; //normalizza elementi sopra diagonale
                    pdata[offset+i] /= akk;
                }

                for (i = k + 1; i < n; i++) {

                    offset=nc*k;

                    for (j = k + 1; j < n; j++) {
                        //a[j][i] -= a[k][i] * a[j][k]; //prodotto scalare fra colonna e riga della sottomatrice
                        a[j][i] -= pdata[offset+i]*a[j][k];
                    }
                }
            }
        }
        

        /**
         * Computes the inverse matrix
         * @param dest destination matrix where the invers of src is copied.
         * @param src matrix to be inverted.
         */
        static void invertMatrix(TMatrix<T>& dest, const TMatrix<T>& src) throw (MathException) {

            /*
        La matrice inversa di una matrice A � una matrice A~ tale che A*A~=I
        si dimostra che la matrice A~ � costituita dagli elementi A~ij
        per cui A~ij = cofattore dell'elemento ij in A / determinante di A
        Il cofattore di un elemento Aij � il determinante della matrice che si
        ottiene eliminando la riga i e la colonna j da A moltiplicato per (-1) ^ (i+j)
        Ogni colonna Xk della matrice inversa A~ si ottiene come colonna soluzione del sistema
        A*Xk=Ck dove Ck � un vettore che il k-esmimo elemento settato a 1 e tutti gli altri settati a zero*/

            if (src.getRows() != dest.getRows()) throw MathException("Cannot determinate the inverse matrix.Matrices have different size.");

            int n = src.getRows();

            LinearEqCrout eqset(n);

            TVector<T> d(n);

            eqset.setEquations(src, d);

            for (int i = 0; i < n; ++i) {
                //ricava la colonna i della matrice inversa
                d.zeroMemory();

                d[i] = 1;

                eqset.setRightHandColumn(d);

                eqset.solve();

                eqset.getResult(&d);

                dest.setColumn(i, d);
            }

        }

        /**
         * Calculates the determinant of matrix m
         * @param m the matrix
         * \warning this method modifies the matrix m.
         * @return the determinant of square matrix m computed with Crout reduction
         */
        static T getDeterminant(TMatrix<T>& m) throw (MathException) {

            /*il determinante della matrice � il prodotto degli elementi
              della diagonale della matrice ridotta con il metodo di Crout*/

            int n = m.getRows();           

            T det = 0;

            try {

                //ottiene la matrice ridotta
                croutReduction(m);

                det = m[0][0];

                for (int i = 1; i < n; ++i) det *= m[i][i];

                return det;

            } catch (MathException& e) {
                //un elemento sulla diagonale � zero quindi � zero anche il determinante
                return 0;
            }
        }

        /**
         * Construct a set of linear equations to be solved with Crout method.
         * @param eq_num number of equations
         */
        LinearEqCrout(int eq_num) : LinearImpl<T>(eq_num),invalid_aux_matrix(true) {

        }

        virtual ~LinearEqCrout() {

        }

        /**
         * Set the coefficients and right-hand member for an equation of the set.
         * @param eq_index index of equation (index of row in the matrix of coefficients)
         * @param coeff vector of coefficients (known members).Must have size EquationsCount
         * @param rhs right-hand member
         * @return the equation's index eq_index
         */
        virtual int setEquation(int eq_index, const T* coeff, T rhs) throw (MathException) {
            invalid_aux_matrix = true;
            LinearImpl<T>::solved = false;
            return LinearImpl<T>::setEquation(eq_index, coeff, rhs);
        }

        /**
         * Set a coefficient of the matrix
         * @param coeff
         * @param eq_index
         * @param column
         */
        inline void setCoefficient(T coeff, int eq_index, int column) throw () {
            invalid_aux_matrix = true;
            LinearImpl<T>::setCoefficient(coeff, eq_index, column);
        }

        /**
         * Force the execution of Crout's reduction to the next call to Solve.
         */
        void invalidateAuxMatrix() {
            invalid_aux_matrix = true;
        }

        /**
         * Checks if the current coefficients matrix has been reduced.
         * If the matrix is reduced,and the rhs changed, the system can be solved without
         * recomputing the reduced matrix.
         *
         * @return true if the matrix has benne reduced
         */
        bool isMatrixReduced() const throw () {
            return !invalid_aux_matrix;
        }

        /**
         * Set the row of coefficients in the augmented matrix (known coefficients + right hand member)
         * @param eq_index index of row
         * @param adjcoeff coefficients + right hand member
         * @return  the index of row
         */
        virtual int setEquation(int eq_index, const T* adjcoeff) throw (MathException) {
            invalid_aux_matrix = true;
            LinearImpl<T>::solved = false;
            return LinearImpl<T>::setEquation(LinearImpl<T>::eqmap[eq_index], adjcoeff);
        }

        /**
         * Set the coefficient matrix and right-hand column
         * @param coeff coefficient matrix
         * @param rhs right hand column (must have the same size as coeff)
         */
        virtual void setEquations(const TMatrix<T>& coeff, const TVector<T>& rhs) throw (MathException) {

            invalid_aux_matrix = true;
            LinearImpl<T>::solved = false;
            LinearImpl<T>::setEquations(coeff, rhs);

        }
        

        /**
         * Solve the equations.
         * This method leaves the coefficients matrix unchanged.
         * @throw MathException if the equations can not be solved.
         */
        virtual void solve() throw (MathException) {
            solve(OPT_NONE);
        }

        /**
         * Solve the equations.
         * This method leaves the coefficients matrix unchanged.
         * @param opt optimizations flags (OPT_NONE,OPT_ROUNDOFFS=1,OPT_DONTARRANGE_MATRIX)
         */
        virtual void solve(int opt) throw (MathException) {

            // bool opt_roundoffs;
            register int i, j, k;
            
            int _n = LinearImpl<T>::n; //numero di equazioni nel sistema

            TMatrix<T> &_mat = LinearImpl<T>::mat; //matrice coefficienti

            TVector<T> &_res = LinearImpl<T>::res; //vettore dei risultati

            //nota: con il metodo di crout in generale non occorre riposizionare
            //le righe se alcuni elementi della diagonale sono zero
            //vanno invece riposizionate se mat[0][0] � zero
            if (!(opt & OPT_DONTARRANGE_MATRIX)) {

                if (invalid_aux_matrix) {
                    
                    LinearImpl<T>::resetEqMap();
                    LinearImpl<T>::pivoting(opt);
                    
                }               
            }            

            T t;
            T inv;

            //la matrice ausiliaria viene ricalcolata solo se sono variati i coefficienti
            //questo consente di risolvere il sistema molto velocemente se varia solo la colonna dei termini noti
            //es. quando si studia coma varia la soluzione al variare dei termini noti
            //se si riarrangiano le equazioni la matrice va invalidata
            if (invalid_aux_matrix) {

                LinearEqCrout::croutReduction(_mat);
                invalid_aux_matrix = false;

            } 
            
            //esegue la riduzione sulla colonna dei termini noti
            //_mat=è la matrice ridotta secondo la regola di Crout
            for (j = 0; j < _n; ++j) {
                
                inv = 1 / _mat[j][j];

                t = 0;
                for (k = 0; k < j; ++k) t += _mat[j][k] * _mat[k][_n];
                _mat[j][_n] = _mat[j][_n] - t;
                _mat[j][_n] *= inv;

            }

            //ricava la soluzione
            for (i = _n - 1; i >= 0; --i) {

                t = 0;

                for (k = i + 1; k < _n; ++k) {
                    t += _mat[i][k] * _res[k];
                }

                _res[i] = _mat[i][_n] - t;

            }

            LinearImpl<T>::solved = true;

        }
        
        /**
         * Get the augmented matrix.
         * \warning returns a reference to an internal matrix.
         * calling this method after solving will return a transformed matrix and not the
         * orginal coefficients matrix.
         * @return the matrix of coefficients + right hand members. (nxn+1 matrix)
         */
        const TMatrix<T>& getAugmentedMatrix() throw (MathException) {
            return LinearImpl<T>::mat;
        }

        const char* getMethodName() const {
            return "Crout's reduction";
        }

    private:
        
        //flag chi indica che la matrice dei coefficienti � cambiata e quindi si deve ricalcolare la matrice
        //ausiliaria
        bool invalid_aux_matrix;
        
    };
};

#endif // CROUT_H_INCLUDED
