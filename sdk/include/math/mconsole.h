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
 \file mconsole.h
 \brief Utility module for matrices and vectors I/O from console.
 */

#ifndef MCONSOLE_H_INCLUDED
#define MCONSOLE_H_INCLUDED

#include "mathdefs.h"
#include "libmath.h"
#include "vectors.h"
#include "vectors2.h"
#include "vectors3.h"
#include "matrix2.h"
#include "quat2.h"
#include <cmath>
#include <iostream>

using namespace std;

#ifdef _USENAMESPACE_
namespace mathengine {
#endif

/**
 * Read a real number form and input stream
 * @param s input stream
 * @param r real number read from s.
 * @return true if the number is valid, false otherwise.
 */
bool readReal(istream& s,REAL& r);

/**
 * Read a comma from input stream.
 * @param s input stream.
 * @return true if the char read is a comma, false otherwise.
 */
bool readComma(istream& s);

/**
 * Read a 3 components vector form an input stream.
 * The vector must have the components separated by comma, eg 3.4,1.2,1.0
 * @param s input stream
 * @param vi vector read.
 * @return
 */
istream& operator>>(istream& s,CHVector3<REAL>& vi);

/**
 * Prints a 2 components vector to an output stream.
 * @param s
 * @param v
 * @return
 */
ostream& operator<<(ostream& s,const VECTOR2D &v);

/**
 * Print a 3 components vector to an output stream.
 * @param s
 * @param v
 * @return
 */
ostream& operator<<(ostream& s,const VECTOR3D &v);

/**
 * Prints a 3 components homogenous vector to an output stream.
 * @param s
 * @param v
 * @return
 */
ostream& operator<<(ostream& s,const CHVector3<REAL>& v);

/**
 * Prints a vector.
 * @param s
 * @param v
 * @return
 */
ostream& operator<<(ostream& s,const CVector2<REAL>& v);

/**
 * Prints a transformation matrix to an output stream
 * @param s
 * @param m
 * @return
 */
ostream& operator<<(ostream& s,const CHMatrix33<REAL>& m);

/**Read a matrix row by row from input stream s*/
istream& operator>>(istream& s,CHMatrix33<REAL>& m);

/**
 * Read a quaternion from an input stream.
 * @param s
 * @param q
 * @return
 */
istream& operator>>(istream& s,CQuaternion<REAL>& q);

/**
 * Prints a quaternion.
 * @param s
 * @param q
 * @return
 */
ostream& operator<<(ostream& s,const CQuaternion<REAL>& q);

#ifdef _USENAMESPACE_
};
#endif

#endif // MCONSOLE_H_INCLUDED
