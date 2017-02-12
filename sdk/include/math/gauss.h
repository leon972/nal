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
 * \file gauss.h
 * \brief Numerical solution of linear algebraic equations with the Gauss method and back-substitution
 */

#ifndef GAUSS_H_INCLUDED
#define GAUSS_H_INCLUDED

#include "libmath.h"
#include "linearimpl.h"
#include "tmatrix.h"
#include <vector>


namespace mathengine {


    using namespace std;

    /**
     * Solve a set of linear algebraic equations with the Gauss method.
     */
    template<class T> class MATH_EXPORT LinearEqGauss : public LinearImpl<T> {
        
    public:

        /**
         * Construct the set
         * @param n number of eqautions
         */
        LinearEqGauss(int eq_num) throw (MathException) : LinearImpl<T>(eq_num) {

            if (eq_num <= 0) throw MathException("invalid equation number");

        }

        /**
         * Destructor
         */
        virtual ~LinearEqGauss() {
        }

        /**
         * Solve the equations with the Gauss method
         * using the default options (OPT_NONE).
         * Equations are arranged before solving in order to minimize roundoff errors.
         * This method calculates also the rank of the system.
         * \warning this method changes the coefficient matrix
         */
        void solve() throw (MathException) {
            solve(OPT_NONE);
        }

        /**
         * Solve the equations with the Gauss method with back substitution
         * \warning this method changes the coefficient matrix
         * @param opt optimizations options:
         * OPT_DONTARRANGE_MATRIX : do not arrange equations in order to minimize roundoff errors
         * OPT_ROUNDOFFS : additional optimiziation to minimize roundoff errors
         * This method calculates also the rank of the system.
         */
        void solve(int opt) throw (MathException) {

            TMatrix<T> &_mat=LinearImpl<T>::mat; //matrice coefficienti
            TVector<T> &_res=LinearImpl<T>::res; //vettore dei risultati
            size_t _n=LinearImpl<T>::n; //numero di equazioni nel sistema
            size_t nc=_mat.getColumns();
            size_t offset,offset1;
            register int i, j, k;
            T t;
            T v;
            T* pdata=_mat._getDataPtr();

            LinearImpl<T>::rank = 0;

            if (!(opt & OPT_DONTARRANGE_MATRIX)) {
                //fa in modo che gli elementi sulla diagonale non siano zero
                LinearImpl<T>::pivoting(opt);
            }            

            for (i = 0; i < _n; ++i) {

                t = 1 / _mat[i][i];
                _mat[i][i] = 1;
                //divide la riga per l'elemento della diagonale
                //usa il calcolo dell'offset per ridurre il numero di moltiplicazioni necessarie
                offset=nc*i;

                for (k = i + 1; k <= _n; ++k) {

                    //_mat[i][k] *= t;
                    pdata[offset+k] *= t;
                }

                //elimina xi dalle equazioni restanti
                for (j = i + 1; j < _n; ++j) {
                    t = _mat[j][i];

                    //il coeff. di xi nell'eq. k viene posto a zero per via dell'eliminazione
                    _mat[j][i] = 0;

                    offset1=nc*j;

                    for (k = i + 1; k <= _n; ++k) {
                        //_mat[j][k] = _mat[j][k] - _mat[i][k] * t;
                        pdata[offset1+k] -= (pdata[offset+k]*t);
                    }
                }

                //fa in modo che non ci siano zeri nella diagonale delle eq. restanti
                LinearImpl<T>::pivoting(opt,i + 1);

                ++LinearImpl<T>::rank; //calcola il rango della matrice
            }

            //a questo punto procede a ritroso per ricavare i risultati
            //(back substitution)
            for (i = _n - 1; i >= 0; --i) {
                t = 0;

                for (k = i + 1; k < _n; ++k) {
                    t += _res[k] * _mat[i][k];
                }

                _res[i] = _mat[i][_n] - t;
            }

            LinearImpl<T>::solved = true;
        }

        const char* getMethodName() const {
            return "Gauss reduction with back-substitution";
        }

    };
};

#endif // GAUSS_H_INCLUDED
