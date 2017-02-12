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
 \file matrix.h
 \brief C-style 2x2,3x3 and 4x4 geometric transformations matrices and algorithms.
 
 */

#ifndef MATRIX33_H_INCLUDED
#define MATRIX33_H_INCLUDED

//questo modulo dipende dal modulo base vectors
//real � un tipo di dati definito in vectors.h
#include "libmath.h"
#include "vectors.h"

#ifdef _USENAMESPACE_
namespace mathengine
{
#endif

#ifdef __cplusplus
//quando si compila in c++ questa serve a fare in modo che non vengano 'decorati'
//i nomi dal linker
extern "C"
{
#endif

/**
 * 2x2 Matrix
 */
typedef struct MATRIX22_TYP
{
    real a00,a01;
    real a10,a11;

} MATRIX22,*MATRIX22_PTR;

/**
 * 3x3 Matrix
 */
typedef struct MATRIX33_TYP
{
    real a00,a01,a02;
    real a10,a11,a12;
    real a20,a21,a22;

} MATRIX33,*MATRIX33_PTR;

/**
 * 3D transformation matrix (rotation+translation)
 * This is a convenience method to represent a rotation+translation
 * 
 */
typedef struct HMATRIX33_TYP
{
    real a00,a01,a02,a03;
    real a10,a11,a12,a13;
    real a20,a21,a22,a23;
    real a30,a31,a32,a33;

} HMATRIX33,*HMATRIX33_PTR;

/**Transform a 2x2 to identity matrix
 *@param m destination matrix
 */
MATH_EXPORT void getIdentityMatrix22Ex(MATRIX22_PTR m);
/**
 * Transform a 3x3 to identity matrix
 * @param m destination matrix
 */
MATH_EXPORT void getIdentityMatrix33Ex(MATRIX33_PTR m);
/**
 * Transform a 3x3 to identity matrix
 * @param m destination matrix
 */
MATH_EXPORT void getIdentityHMatrix33Ex(HMATRIX33_PTR m);

/**
 * Create a 2x2 identity matrix
 * @return the 2x2 identity matrix
 */
MATH_EXPORT MATRIX22 getIdentityMatrix22(void);
/**
 * Create a 3x3 identity matrix
 * @return a 3x3 identity matrix
 */
MATH_EXPORT MATRIX33 getIdentityMatrix33(void);

/**
 * Create a 3x3 identity matrix
 * @return
 */
MATH_EXPORT HMATRIX33 getIdentityHMatrix33(void);

/**
 * Matrices addition
 * @param result the sum m1+m2
 * @param m1 first matrix
 * @param m2 second matrix
 */
MATH_EXPORT void addMatrix33(MATRIX33_PTR result,const MATRIX33& m1,const MATRIX33&m2);

/**
 * Add a 3x3 matrix to another matrix
 * @param dest destination matrix (dest=dest+m)
 * @param m matrix to add
 */
MATH_EXPORT void addMatrix33Ex(MATRIX33_PTR dest,const MATRIX33& m);

/**
 * Matrices addition
 * @param result the sum of m1 and m2
 * @param m1 first matrix
 * @param m2 second matrix
 */
MATH_EXPORT void addHMatrix33Ex(HMATRIX33_PTR result,const HMATRIX33& m1,const HMATRIX33&m2);

/**
 * Add a 3x3 (omogeneous coordinates) transform matrix to another matrix
 * @param dest destination matrix (dest=dest+m
 * @param m matrix to add
 */
MATH_EXPORT void addHMatrix33(HMATRIX33_PTR dest,const HMATRIX33& m);

/**
 * Matrix subtraction
 * @param result the difference between m1 and m2 (result=m1-m2)
 * @param m1 first matrix
 * @param m2 second matrix
 */
MATH_EXPORT void subMatrix33(MATRIX33_PTR result,const MATRIX33& m1,const MATRIX33&m2);

/**
 * Subtract a matrix from another matrix
 * @param dest destination matrix (dest=dest-m)
 * @param m matrix to subtract
 */
MATH_EXPORT void subMatrix33Ex(MATRIX33_PTR dest,const MATRIX33& m);

/**
 * Subtraction between transformation matrices (homogeneous coordinates)
 * @param result the difference between m1 and m2 (result=m1-m2)
 * @param m1 first matrix
 * @param m2 second matrix
 */
MATH_EXPORT void subHMatrix33(HMATRIX33_PTR result,const HMATRIX33& m1,const HMATRIX33&m2);

/**
 * Subtract a transformation matrix from another transformation matrix
 * @param dest destination matrix (dest=dest-m)
 * @param m transformation matrix to subtract
 */
MATH_EXPORT void subHMatrix33Ex(HMATRIX33_PTR dest,const HMATRIX33& m);

/**
 * Compute the determinant of a 3x3 matrix
 * @param m 3x3 matrix
 * @return the determinant of m
 */
MATH_EXPORT real det(const MATRIX33& m);

/**
 * Compute the inverse of a 3x3 matrix
 * @param r destination matrix (in r will be copied the inverse of m)
 * @param m matrix to invert
 * @return true if the inverse matrix exeists,false otherwise
 */
MATH_EXPORT bool inverse(MATRIX33_PTR r,const MATRIX33& m);


/**
 * Solves a system of three equations using the Cramer's rule
 * @param x result vector
 * @param aij matrix of coefficients
 * @param d knowns column
 * @return true if the system can be solved false otherwise (if the determinant of aij equals zero)
 */
MATH_EXPORT bool solveEquationsSet3(VECTOR3D_PTR x,const MATRIX33& aij,const VECTOR3D& d);


/**
 * Solves a system of two equations using the Cramer's rule
 * @param x result vector
 * @param aij matrix of coefficients
 * @param d knowns column
 * @return true if the system can be solved false otherwise (if the determinant of aij equals zero)
 */
MATH_EXPORT bool solveEquationsSet2(VECTOR2D_PTR x,const MATRIX22& aij,const VECTOR2D& d);

/**
 * Multiply a scalar by a 3x3 matrix
 * @param result m*a
 * @param a a scalar
 * @param m the matrix to be multiplied by a
 */
MATH_EXPORT void mulScalMatrix33(MATRIX33_PTR result,real a,const MATRIX33& m);


/**
 * Multiply a scalar by a 3x3 matrix in place
 * @param dest the destination matrix (m=m*a)
 * @param a a scalar
 *
 */
MATH_EXPORT void mulScalMatrix33Ex(MATRIX33_PTR dest,real a);

/**
 * Multiply a scalar by a transformation matrix
 * @param result m*a
 * @param a a scalar
 * @param m the matrix to be multiplied
 */
MATH_EXPORT void mulScalHMatrix33(HMATRIX33_PTR result,real a,const HMATRIX33& m);

/**
 * Multiply a transformation matrix by a scalar in place
 * @param dest destination matrix (dest=dest*a)
 * @param a a scalar
 */
MATH_EXPORT void mulScalHMatrix33Ex(HMATRIX33_PTR dest,real a);

/**
 * 3x3 matrix multiplication
 * @param result m1*m2
 * @param m1 first matrx
 * @param m2 second matrix
 */
MATH_EXPORT void mulMatrix33(MATRIX33_PTR result,const MATRIX33& m1,const MATRIX33& m2);

/**
 * Transformation matrices multiplication.
 * This function is used for subsequent transformations.
 * @param result m1*m2,the transformation applied from m2 to m1 referred to local axes
 * @param m1 first matrix
 * @param m2 second matrix
 */
MATH_EXPORT void mulHMatrix33(HMATRIX33_PTR result,const HMATRIX33& m1,const HMATRIX33& m2);

/**
 * Compute the transpose transformationmatrix
 * @param result the transpose of m  , result->aij=m.aji , i=0..3 j=0..4
 * @param m the original matrix
 */
MATH_EXPORT void transposeMatrixEx(HMATRIX33_PTR result,const HMATRIX33& m);

/**
 * Transpose matrix
 * @param m the matrix to transpose
 * @return the transpose of m
 */
MATH_EXPORT HMATRIX33 transposeMatrix(const HMATRIX33& m);

/**
 * Multiply a 2x2 matrix by a 2 components vector
 * @param result m*v
 * @param m 2x2 matrix
 * @param v 2 components vector
 */
MATH_EXPORT void mulMatrixVect2(VECTOR2D_PTR result,const MATRIX22& m,const VECTOR2D& v);

/**
 * Multiply a 3x3 matrix by a 3 components vector
 * @param result m*v
 * @param m 3x3 matrix
 * @param v 3 components vector
 */
MATH_EXPORT void mulMatrixVect3(VECTOR3D_PTR result,const MATRIX33& m,const VECTOR3D& v);

/**
 * Multiply a transform matrix by a vector
 * @param result m*v, the result vector is the vector expressed in the frame of reference defined by m
 * @param m transformation matrix
 * @param v 3 compoenents vector in homogeneous coordinates
 */
MATH_EXPORT void mulMatrixHVect3Ex(HVECTOR3D_PTR result,const HMATRIX33& m,const HVECTOR3D& v);

/**
 * Multiply a transform matrix by a 3 components vector in homogeneous coordinates
 * @param m transform matrix
 * @param v 3 components vector in homogeneous coordinates 
 * @return m*v
 */
MATH_EXPORT HVECTOR3D mulMatrixHVect3(const HMATRIX33& m,const HVECTOR3D& v);

/**
 * create a 3x3 matrix from a transformation matrix (copied the rotation part of to result)
 * @param result
 * @param m
 */
MATH_EXPORT void toMatrix33(MATRIX33_PTR result,const HMATRIX33& m);

/**
 * Create a transformation matrix from a 3x3 matrix
 * @param result the transformation matrix.Translation part is set to zero, resul.a33=1
 * @param m 3x3 matrix
 */
MATH_EXPORT void toHMatrix33(HMATRIX33_PTR result,const MATRIX33& m);

#ifdef __cplusplus
}
#endif

#ifdef _USENAMESPACE_
};
#endif

#endif // MATRIX33_H_INCLUDED
