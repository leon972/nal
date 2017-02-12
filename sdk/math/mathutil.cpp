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

#include "math/mathutil.h"
#include <cmath>

#ifdef _USENAMESPACE_
namespace mathengine {
#endif

    /**
     * Determina le soluzioni reali di una equazione di secondo grado
     * del tipo ax^2+bx+c=0
     * @param x1 usato per ottenere la prima soluzione
     * @param x2 usato per ottenere la seconda soluzione
     * @param a coefficiente di x quadro
     * @param b coefficiente di x
     * @param c termine noto
     * @throw NoSolutionException se non ammette soluzioni reali
     */
    void solveQuadraticEquation(double *x1, double *x2, const double a, const double b, const double c) throw (NoSolutionException) {

        double d = b * b - 4 * a*c;

        if (d < 0) throw NoSolutionException(); //non ammette soluzioni reali

        if (a == 0) {
            //si riduce ad una eq. lineare
            if (b==0) throw NoSolutionException();

            if (x1) {
                *x1=-c/b;
            }
            if (x2) {
                *x2=-c/b;
            }
        }
        else {

            //nota: se a e c sono molto piccoli si incorre in errori di approssimazione
            //in quanto si sottrarrebbe da b un numero quasi uguale a b (essende a e c picoli)

            double sgn = (b >= 0 ? 1 : -1);
            double q = -(b + sgn * sqrt(d))*0.5;

            if (x1) {
                *x1 = q / a;
            }

            if (x2) {
                *x2 = c / q;
            }

        }
    }

    /**
     *Risolve un sistema lineare di 2 equazioni in 2 incognite usando la regola di Cramer
     *@param aij = coefficienti matrice
     *@param di = riga termini noti
     *@param x1,x2 parametri output in cui viene restituito il risultato
     */
    void solveLinearEquations(double *x1, double *x2, const double a11, const double a12, const double a21, const double a22, const double d1, const double d2) throw (NoSolutionException) {

        //deteminante matrice
        double D = a11 * a22 - a21*a12;

        //il determinante è nullo, il sistema ammette infinite soluzioni
        if (D == 0) throw NoSolutionException();

        if (x1) {
            *x1 = (a22 * d1 - a12 * d2) / D;
        }

        if (x2) {
            *x2 = (a11 * d2 - a21 * d1) / D;
        }
    }

#ifdef _USENAMESPACE_
};
#endif

