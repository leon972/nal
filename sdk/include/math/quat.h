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
 \file quat.h
 \brief C-style quaternions and quaternions algorithms. 
 */

#ifndef QUAT_H_INCLUDED
#define QUAT_H_INCLUDED

//questo modulo dipende dal modulo base vectors
//real e' un tipo di dati definito in vectors.h
#include "libmath.h"
#include "vectors.h"
#include "matrix.h"

#ifdef _USENAMESPACE_
namespace mathengine
{
#endif

/**QUAT3D - quaternion 3D*/
typedef struct QUAT3D_TYP
{
    real r; /*!<-real part*/
    real x,y,z; /*!<-imaginary part*/

} QUAT3D,*QUAT3D_PTR;

#ifdef __cplusplus
//quando si compila in c++ questa serve a fare in modo che non vengano 'decorati'
//i nomi dal linker
extern "C"
{
#endif

/**
 * Print a quaternion to file
 * eg.  (r+xi+yj+zk)
 * @param out output file already open
 * @param q quaternion to print
 */
MATH_EXPORT void printQuat3D(FILE* out,const QUAT3D& q);

/**
 * create a quaternion
 * @param r real part
 * @param x
 * @param y
 * @param z
 * @return the quaternion r+ix+yj+zk
 */
MATH_EXPORT QUAT3D getQuat3D(const real r,const real x,const real y,const real z);

/**
 * Create a quaternion using a vector for the imaginary part
 * @param v vector,
 * @return the quaternion 0+v.x*i+v.y*j+v.z*k
 */
MATH_EXPORT QUAT3D vectorToQuat3D(HVECTOR3D_PTR v);

/**
 * Create a quaternion in place
 * @param q destination quaternion
 * @param v vector used as imaginary part
 */
MATH_EXPORT void vectorToQuat3DEx(QUAT3D_PTR q,const HVECTOR3D& v);

/**
 * Create a vector from a quaternion using the imaginary part
 * @param q quaternion
 * @return the vector q.x,q.y,q.z,1
 */
MATH_EXPORT HVECTOR3D quatToVector(const QUAT3D& q);

/**
 * Create a vector from a quaternion using the imaginary part
 * @param v destination vector
 * @param q quaternion
 */
MATH_EXPORT void quatToVectorEx(HVECTOR3D_PTR v,const QUAT3D& q);

/**
 * Quaternion addition
 * @param q1 first addend
 * @param q2 second addend
 * @return  q1+q2    (q1.r+q2.r,q1.x+q2.x,q1.y+q2.y,q1.z+q2.z)
 */
MATH_EXPORT QUAT3D sumQuat3D(const QUAT3D& q1,const QUAT3D& q2);

/**
 * Quaternion subtraction
 * @param result q1-q2  (q1.r-q2.r,q1.x-q2.x,q1.y-q2.y,q1.z-q2.z)
 * @param q1
 * @param q2
 */
MATH_EXPORT void sumQuat3DEx(QUAT3D_PTR result,const QUAT3D& q1,const QUAT3D& q2);

/**
 * Quaternion subtraction
 * @param q1
 * @param q2
 * @return q1-q2 (q1.r-q2.r,q1.x-q2.x,q1.y-q2.y,q1.z-q2.z)
 */
MATH_EXPORT QUAT3D subQuat3D(const QUAT3D& q1,const QUAT3D& q2);

/**
 * Quaternion subtraction
 * @param result the destination quaternion (q1-q2) (q1.r-q2.r,q1.x-q2.x,q1.y-q2.y,q1.z-q2.z)
 * @param q1
 * @param q2
 */
MATH_EXPORT void subQuat3DEx(QUAT3D_PTR result,const QUAT3D& q1,const QUAT3D& q2);

/**
 * Quaternion product
 * @param q1
 * @param q2
 * @return q1*q2
 */
MATH_EXPORT QUAT3D mulQuat3D(const QUAT3D& q1,const QUAT3D& q2);

/**
 * Quaternion product
 * @param result result quaternion, q1*q2
 * @param q1
 * @param q2
 */
MATH_EXPORT void mulQuat3DEx(QUAT3D_PTR result,const QUAT3D& q1,const QUAT3D& q2);

/**
 * Quaternion conjugate
 * @param q quaternion
 * @return the conjugate of q  (q.r,-q.xi,-q.yi,-q.zk)
 */
MATH_EXPORT QUAT3D quatConjugate3D(const QUAT3D& q);

/**
 * Quaternion conjugate
 * @param result destination quaternion,will hold the conjugate of q
 * @param q quaternion
 */
MATH_EXPORT void quatConjugateEx(QUAT3D_PTR result,const QUAT3D& q);

/**
 * Quaternion normalisation
 * Return a quaternion of unit length
 * @param q
 * @return q normalised
 */
MATH_EXPORT QUAT3D normalizeQuat3D(const QUAT3D& q);

/**
 * Quaternion normalisation
 * @param q
 */
MATH_EXPORT void normalizeQuat3DEx(QUAT3D_PTR q);

/**
 * Compute the magnitude of a quaternion
 * @param q
 * @return the magnitude of q
 */
MATH_EXPORT real magnitudeQuat3D(const QUAT3D& q);

/**
 * Compute the magnitude of a quaternion
 * @param q
 * @return the magnitude of q
 */
MATH_EXPORT real magnitudeQuat3DEx(QUAT3D_PTR q);

/**
 * Create a quaternion representing the rotation of theta radians
 * about the axis whose versor is versAxis
 * @param versAxis a vector of unit length representin the axis of rotation
 * @param theta rotation in radians
 * @return the quaternion representing a rotation. (cos(theta/2)+versAxis*(sin(theta/2))
 */
MATH_EXPORT QUAT3D getRotQuat3D(const HVECTOR3D& versAxis,real theta);

/**
 * Create a quaternion representing the rotation of theta radians
 * about the axis whose versor is versAxis
 * @param q result quaternion (not null) (cos(theta/2)+versAxis*(sin(theta/2))
 * @param versAxis a vector of unit length representin the axis of rotation
 * @param theta rotation in radians
 */
MATH_EXPORT void getRotQuat3DEx(QUAT3D_PTR q,const HVECTOR3D& versAxis,real theta);

/**
 * Compute the rotation matrix corresponding to quaternion q
 *
 * @param m result matrix
 * @param q orientation expressed using a quaternion
 */
MATH_EXPORT void quatToRotMatrixEx(HMATRIX33_PTR m,const QUAT3D& q);

/**
 * Compute the rotation matrix corresponding to quaternion q
 * @param q orientation expressed using a quaternion
 * @return the rotation matrix
 */
MATH_EXPORT HMATRIX33 quatToRotMatrix(const QUAT3D& q);

/**
 * Convert a rotation matrix to a quaternion representing the same
 * rotation
 * @param q result quaternion
 * @param m
 */
MATH_EXPORT void rotMatToQuatEx(QUAT3D_PTR q,const HMATRIX33& m);

/**
 * Convert a rotation matrix to a quaternion 
 * @param m
 * @return the quaternion representing the same rotation of m
 */
MATH_EXPORT QUAT3D rotMatToQuat(const HMATRIX33& m);

/**
 * Quaternion SLERP (spherical linear interpolation)
 * Interpolate between two quaternions representing rotations.
 * \note This method is useful to achieve smooth interpolation between
 * any two orientations.The first quaternion represents the original orientation
 * of the axes, the second quaternion the final orientation of the axes.
 *
 * @param r interpolated quaternion between qa and qb
 * @param qa first quaternion
 * @param qb second quaternion
 * @param t a scalar between 0.0 and 1.0
 */
MATH_EXPORT void quatSlerp(QUAT3D_PTR r,const QUAT3D& qa,const QUAT3D& qb,real t);

/**
 * Quaternion triple poduct
 * \note The triple product can be used to rotate a vector about any axis.
 * If the quaternion q defines the axis of rotation and rotation angle (@see getRotQuat3D)
 * anf v the vector to be rotated, then the imaginary part of q * v * q' represents the rotated vector
 * v' about the axis. q' is the conjugate of q.
 * @param result q1*q2*q3
 * @param q1
 * @param q2 
 * @param q3
 */
MATH_EXPORT void quatTripleProduct(QUAT3D_PTR result,const QUAT3D& q1,const QUAT3D& q2,const QUAT3D& q3);

#ifdef __cplusplus
}
#endif

#ifdef _USENAMESPACE_
};
#endif

#endif // QUAT_H_INCLUDED
