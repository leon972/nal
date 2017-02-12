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
 * \file mathutil.h
 * \brief Math utility functions.

 */


#ifndef MATHUTIL_INCLUDE
#define MATHUTIL_INCLUDE

#include "mathdefs.h"
#include "libmath.h"

#ifdef _USENAMESPACE_
/*!
 * \brief main namespace of the math engine
 */
namespace mathengine {
#endif

    /**
     * Finds the real roots of a quadratic equations
     * ax^2+bx+c=0
     * @param x1 out parameter:returns the first root
     * @param x2 out parameter:returns the second root
     * @param a coefficient of the quadratic term x*x
     * @param b coefficient of the linear term x
     * @param c constant term
     * @throw NoSolutionException if the equation doesn't have real roots
     */
    MATH_EXPORT void solveQuadraticEquation(double *x1,double *x2,const double a,const double b,const double c) throw (NoSolutionException);

    /**
     *Fast solution of a linear algebraic equations set of two equations.
     *aij = coefficients (i index of row , j =index of column)
     *@param d1 1st known term of the right hand column.
     *@param d2 2nd known term of the right hand column.
     *@param x1,x2 out parameters:returns the roots.
     */
    MATH_EXPORT void solveLinearEquations(double *x1,double *x2,const double a11,const double a12,const double a21,const double a22,const double d1,const double d2) throw (NoSolutionException);
   

    /**
     * General purpose formatter interface.
     * Used for output and input of numeric entities
     */
    template <class T> class NumericFormatter {
    public:

        /**
         * Convert a numeric entity into a string
         * @param value entity to convert (eg. vector,double,matrix ...)
         * @return the string representing the entity
         */
        virtual const char* format(const T &value) = 0;

        /**
         * Convert a string into a numeric entity of type T
         * @param dest destination numeric entity
         * @param source string representing the entity (eg. a vector,a double ...)
         * @return
         */
        virtual bool parse(T &dest, const char* source) = 0;

        
        virtual ~NumericFormatter() {
        }

    };


#ifdef _USENAMESPACE_
};
#endif

#endif


