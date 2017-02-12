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
 \file vectors.h
 \brief 2D,3D and 3D homogeneous C-style vectors and vectors algorithms.
 */

#include <stdio.h>
#include <math.h>
#include "libmath.h"
#include "mathdefs.h"

#ifndef _VECTORS_INCLUDE_
#define _VECTORS_INCLUDE_


// ---- engine switches ------------------------------------------------
// modificare queste definizioni per adattare l'engine al compilatore

//indica di utilizzare i namespace (per default � definito)
//#define _USENAMESPACE_

//indica di definire utilizzare il tipo double per rappresentare i numeri
//reali (per default non � definito e usa i float)
//#define _USEDOUBLE_



#ifdef _USENAMESPACE_

namespace mathengine {

#endif

#ifdef __cplusplus
    //quando si compila in c++ questa serve a fare in modo che non vengano 'decorati'
    //i nomi dal linker
    extern "C" {
#endif


#ifdef _USEDOUBLE_

        //utilizza il tipo double per rappresentare i numeri reali
        typedef double real;
        //definisce le funzioni matemetiche da usare con il tipo double
#define sinr(x) sin(x)
#define cosr(x) cos(x)
#define sqrtr(x) sqrt(x)
#define acosr(x) acos(x)

#else

        //utilizza il tipo float per rappresentare i numeri reali (default perch� +veloce)
        typedef float real;
        //usa la versione float per le seguenti funzioni matematiche
#define sinr(x) sinf(x)
#define cosr(x) cosf(x)
#define sqrtr(x) sqrt(x)
#define acosr(x) acosf(x)

#endif

        /**         
         *2D Vector type.
         */
        typedef struct VECTOR2D_TYP {

            VECTOR2D_TYP() {

            }

            VECTOR2D_TYP(real x, real y) {
                this->x = x;
                this->y = y;
            }

            real x, y;
        } VECTOR2D, *VECTOR2D_PTR, POINT2D, *POINT2D_PTR;

        /**
         *3D vector type.
         */
        typedef struct VECTOR3D_TYP {

            VECTOR3D_TYP() {
            }

            VECTOR3D_TYP(real x, real y, real z) {
                this->x = x;
                this->y = y;
                this->z = z;
            }

            real x, y, z;

        } VECTOR3D, *VECTOR3D_PTR, POINT3D, *POINT3D_PTR;

        /**
         *3D vector (homogeneous 4 components).
         *This vector type is used with transformation matrix HMATRIX33
         */
        typedef struct HVECTOR3D_TYP {
            real x; /**<- x component.*/
            real y; /**<- y component */
            real z; /**<- z component */
            real w; /**<- homogeneous component (always 1)*/

        } HVECTOR3D, *HVECTOR3D_PTR, HPOINT3D, *HPOINT3D_PTR;

        /**
         * Create a 3D vector.
         * @param x x component
         * @param y y component
         * @return il vettore 2D
         */
        MATH_EXPORT  VECTOR2D getVector2D(real x, real y);

        /**
         * Create a 3D vector.
         * @param x x component
         * @param y y component
         * @param z z component
         * @return the 3D vector (x,y,z)
         */
        MATH_EXPORT VECTOR3D  getVector3D(real x, real y, real z);

        /**
         * Create a 3D vector with homogeneous components.
         * @param x x component
         * @param y y component
         * @param z z component
         * @return the 3D vector (x,y,z,1)
         */
        MATH_EXPORT HVECTOR3D  getHVector3D(real x, real y, real z);

        /**
         * Set all components to zero.
         * @param vec2 2D vector.
         */
        MATH_EXPORT  void zeroVector2(VECTOR2D &vec2);

        /**
         * Set all components to zero.
         * @param vec3 3D vector.
         */
        MATH_EXPORT  void zeroVector3(VECTOR3D &vec3);

        /**
         * Set all components to zero. (w is set to 1)
         * @param vec3 3D vector.
         */
        MATH_EXPORT  void zeroVector3H(HVECTOR3D &vec3);

        /**
         * Copy x,y component to a 2D vector
         * @param result destination vector
         * @param x x component
         * @param y t component
         */
        MATH_EXPORT  void getVector2DEx(VECTOR2D_PTR result, real x, real y);

        /**
         * Copy x,y,z component to a 3D vector.
         */
        MATH_EXPORT  void getVector3DEx(VECTOR3D_PTR result, real x, real y, real z);


        /**
         * Copy x,y,z component to a 3D vector.
         * @param result destination vector.
         * @param x
         * @param y
         * @param z
         */
        MATH_EXPORT  void getHVector3DEx(HVECTOR3D_PTR result, real x, real y, real z);


        /**
         * Print 2D vector to file.
         * Vector is printed using the format (x,y)
         * @param out destination file.
         * @param v vector to be printed.
         */
        MATH_EXPORT void printVector2D(FILE* out, const VECTOR2D& v);


        /**
         * Print 3D vector to file.
         * Vector is printed using the format (x,y,z)
         * @param out destination file.
         * @param v vector to be printed.
         */
        MATH_EXPORT void printVector3D(FILE* out, const VECTOR3D& v);

        
        /**
         * Print 3D vector to file.
         * Vector is printed using the format (x,y,z)
         * @param out destination file.
         * @param v vector to be printed.
         */
        MATH_EXPORT void printVectorH3D(FILE* out, const HVECTOR3D& v);


        /**
         * Set all components to 1.
         * @param result destination vector.
         * result.x=1 and result.y=1 after this method call.
         */
        MATH_EXPORT  void getUnitVector2D(VECTOR2D_PTR result);

        /**
         * Set all components to 1.
         * @param result destination vector.
         * result.x=1 and result.y=1 and result.z=1 after this method call.
         */
        MATH_EXPORT  void getUnitVector3D(VECTOR3D_PTR result);

        /**
         * Set all components to 1.
         * @param result destination vector.
         * result.x=1 and result.y=1 and result.z=1 after this method call.
         */
        MATH_EXPORT  void getUnitHVector3D(HVECTOR3D_PTR result);

        /**
         * Get the opposite vector.
         * @param v original vector.
         * @return the opposite vector of v.
         */
        MATH_EXPORT  VECTOR2D getNegVector2D(const VECTOR2D& v);

        /**
         * Get the opposite vector.
         * @param v original vector.
         * @return the opposite vector of v.
         */
        MATH_EXPORT  VECTOR3D getNegVector3D(const VECTOR3D& v);

        /**
         * Get the opposite vector.
         * @param v original vector.
         * @return the opposite vector of v.
         */
        MATH_EXPORT  HVECTOR3D getNegHVector3D(const HVECTOR3D& v);

        /**
         * Create the opposite vector in place.
         * v.x=-v.x
         * v.y=-v.y ...
         * @param v destination vector.This is trasformed to its opposite.
         */
        MATH_EXPORT  void negVector2D(VECTOR2D_PTR v);

        /**
         * Create the opposite vector in place.
         * @param v destination vector.This is trasformed to its opposite.
         */
        MATH_EXPORT  void negVector3D(VECTOR3D_PTR v);

        /**
         * Create the opposite vector in place.
         * @param v destination vector.This is trasformed to its opposite.
         */
        MATH_EXPORT  void negHVector3D(HVECTOR3D_PTR v);

        /**
         * 2D vector addition.
         * Add v to dest.
         * @param dest destination vector. dest=dest+v
         * @param v vector to be added
         */
        MATH_EXPORT  void addVector2D(VECTOR2D_PTR dest, const VECTOR2D& v);

        /**
         * 3D vector addition.
         * Add v to dest.
         * @param dest destination vector. dest=dest+v
         * @param v vector to be added
         */
        MATH_EXPORT  void addVector3D(VECTOR3D_PTR dest, const VECTOR3D& v);

        /**
         * 3D vector addition.
         * Add v to dest.
         * @param dest destination vector. dest=dest+v
         * @param v vector to be added
         */
        MATH_EXPORT  void addHVector3D(HVECTOR3D_PTR dest, const HVECTOR3D& v);

        /**
         * 2D Vector addition
         * @param result sum of v1 and v2 (result=v1+v2)
         * @param v1 first vector
         * @param v2 second vector
         */
        MATH_EXPORT void sumVector2DEx(VECTOR2D_PTR result, const VECTOR2D& v1, const VECTOR2D& v2);

        /**
         * 3D Vector addition
         * @param result sum of v1 and v2 (result=v1+v2)
         * @param v1 first vector
         * @param v2 second vector
         */
        MATH_EXPORT void sumVector3DEx(VECTOR3D_PTR result, const VECTOR3D& v1, const VECTOR3D& v2);

        /**
         * 3D Vector addition
         * @param result sum of v1 and v2 (result=v1+v2)
         * @param v1 first vector
         * @param v2 second vector
         */
        MATH_EXPORT void sumHVector3DEx(HVECTOR3D_PTR result, const HVECTOR3D& v1, const HVECTOR3D& v2);

        /**
         * 2D Vector addition.
         * @param v1
         * @param v2
         * @return v1+v2
         */
        MATH_EXPORT VECTOR2D sumVector2D(const VECTOR2D& v1, const VECTOR2D& v2);

        /**
         * 3D Vector addition.
         * @param v1
         * @param v2
         * @return v1+v2
         */
        MATH_EXPORT VECTOR3D sumVector3D(const VECTOR3D& v1, const VECTOR3D& v2);

        /**
         * 3D Vector addition.
         * @param v1
         * @param v2
         * @return v1+v2
         */
        MATH_EXPORT HVECTOR3D sumHVector3D(const HVECTOR3D& v1, const HVECTOR3D& v2);

        /**
         * 2D Vector subtraction
         * @param result v1-v2
         * @param v1
         * @param v2
         */
        MATH_EXPORT void subVector2DEx(VECTOR2D_PTR result, const VECTOR2D& v1, const VECTOR2D& v2);

        /**
         * 3D Vector subtraction
         * @param result v1-v2
         * @param v1
         * @param v2
         */
        MATH_EXPORT void subVector3DEx(VECTOR3D_PTR result, const VECTOR3D& v1, const VECTOR3D& v2);

        /**
         * 3D Vector subtraction
         * @param result v1-v2
         * @param v1
         * @param v2
         */
        MATH_EXPORT void subHVector3DEx(HVECTOR3D_PTR result, const HVECTOR3D& v1, const HVECTOR3D& v2);

        /**
         * 2D Vector subtraction
         * @param v1
         * @param v2
         * @return v1-v2
         */
        MATH_EXPORT VECTOR2D subVector2D(const VECTOR2D& v1, const VECTOR2D& v2);

        /**
         * 2D Vector subtraction
         * @param v1
         * @param v2
         * @return v1-v2
         */
        MATH_EXPORT VECTOR3D subVector3D(const VECTOR3D& v1, const VECTOR3D& v2);

        /**
         * 3D Vector subtraction
         * @param v1
         * @param v2
         * @return v1-v2
         */
        MATH_EXPORT HVECTOR3D subHVector3D(const HVECTOR3D& v1, const HVECTOR3D& v2);

        /**
         * 3D Scalar by vector multiplication.
         * @param v
         * @param scalar
         * @param result destination vector where the result is copied.
         */
        MATH_EXPORT void mulHVector3DEx(HVECTOR3D_PTR result, const HVECTOR3D&v, const real scalar);

        /**
         * Scalar by vector product
         * @param v
         * @param scalar
         * @return scalar*v
         */
        MATH_EXPORT HVECTOR3D mulHVector3D(const HVECTOR3D&v, const real scalar);

        /**
         * 2D vectors dot product
         * @param v1
         * @param v2
         * @return v1 . v2
         */
        MATH_EXPORT real dotProduct2(const VECTOR2D& v1, const VECTOR2D& v2);

        /**
         * 3D vectors dot product
         * @param v1
         * @param v2
         * @return v1 . v2
         */
        MATH_EXPORT real dotProduct3(const VECTOR3D& v1, const VECTOR3D& v2);

        /**
         * 3D vectors dot product
         * @param v1
         * @param v2
         * @return v1 . v2
         */
        MATH_EXPORT real dotProduct3H(const HVECTOR3D& v1, const HVECTOR3D& v2);

        /**
         * 3D vectors cross product
         * @param v1
         * @param v2
         * @return v1 x v2
         */
        MATH_EXPORT VECTOR3D crossProduct3(const VECTOR3D& v1, const VECTOR3D& v2);

        /**
         * 3D vectors cross product
         * @param v1
         * @param v2
         * @return v1 x v2
         */
        MATH_EXPORT HVECTOR3D crossProduct3H(const HVECTOR3D& v1, const HVECTOR3D& v2);


        /**
         * 3D vectors cross product
         * @param v1
         * @param v2
         * @param result destination vector, where the cross product is copied.
         */
        MATH_EXPORT void crossProduct3Ex(VECTOR3D_PTR result, const VECTOR3D& v1, const VECTOR3D& v2);

        /**
         * 3D vectors cross product
         * @param v1
         * @param v2
         * @param result destination vector where the result of v1 x v2 is copied.
         */
        MATH_EXPORT void crossProduct3HEx(HVECTOR3D_PTR result, const HVECTOR3D& v1, const HVECTOR3D& v2);

        /**
         * Normalized cross product.
         * Compute the versor of cross product.
         * The versor returned in perpendicular to v1 and v2.
         * @param v1
         * @param v2
         * @return the versor (vector of magnitude 1) of the cross product v1 x 2.
         */
        MATH_EXPORT VECTOR3D normCrossProduct3(const VECTOR3D& v1, const VECTOR3D& v2);

        /**
         * Normalized cross product.
         * Compute the versor of cross product.
         * The versor returned in perpendicular to v1 and v2.
         * @param v1
         * @param v2
         * @return the versor (vector of magnitude 1) of the cross product v1 x 2.
         */
        MATH_EXPORT HVECTOR3D normCrossProduct3H(const HVECTOR3D& v1, const HVECTOR3D& v2);

        /**
         * Normalized cross product.
         * Compute the versor of cross product.
         * The versor returned in perpendicular to v1 and v2.
         * @param result destination vector.After this method call,result is a versor normal to v1 and v2.
         * @param v1
         * @param v2
         */
        MATH_EXPORT void normCrossProduct3HEx(HVECTOR3D_PTR result, const HVECTOR3D& v1, const HVECTOR3D& v2);

        /**
         * Normalized cross product.
         * Compute the versor of cross product.
         * The versor returned in perpendicular to v1 and v2.
         * @param result destination vector.After this method call,result is a versor normal to v1 and v2.
         * @param v1
         * @param v2
         */
        MATH_EXPORT void normCrossProduct3Ex(VECTOR3D_PTR result, const VECTOR3D& v1, const VECTOR3D& v2);

        /**
         * Get the versor (vector of magnitude 1)
         * @param v vector to be normalized.
         * @return the versor of v
         */
        MATH_EXPORT VECTOR2D normalizeVector2D(const VECTOR2D& v);


        /**
         * Get the versor (vector of magnitude 1)
         * @param v vector to be normalized.
         * @return the versor of v
         */
        MATH_EXPORT VECTOR3D normalizeVector3D(const VECTOR3D& v);

        /**
         * Get the versor (vector of magnitude 1)
         * @param v vector to be normalized.
         * @return the versor of v
         */
        MATH_EXPORT HVECTOR3D normalizeHVector3D(const HVECTOR3D& v);

        /**
         * Normalize vector in place.
         * v is a versor of unit length after this method call.
         * @param v vector to be normalized
         */
        MATH_EXPORT void normalizeVector2DEx(VECTOR2D_PTR v);

        /**
         * Normalize vector in place.
         * v is a versor of unit length after this method call.
         * @param v vector to be normalized
         */
        MATH_EXPORT void normalizeVector3DEx(VECTOR3D_PTR v);


        /**
         * Normalize vector in place.
         * v is a versor of unit length after this method call.
         * @param v vector to be normalized
         */
        MATH_EXPORT void normalizeHVector3DEx(HVECTOR3D_PTR v);

        /**
         * Compute the magnitude of a 2D vector
         * @param v
         * @return the magnitude of v
         */
        MATH_EXPORT real magnitude2(const VECTOR2D& v);

        /**
         * Compute the magnitude of a 2D vector
         * @param v
         * @return the magnitude of v
         */
        MATH_EXPORT real magniture2Ex(VECTOR2D_PTR v);

        /**
         * Compute the magnitude of a 3D vector
         * @param v
         * @return the magnitude of v
         */
        MATH_EXPORT real magnitude(const VECTOR3D& v);

        /**
         * Compute the magnitude of a 3D vector
         * @param v
         * @return the magnitude of v
         */
        MATH_EXPORT real magnitudeEx(VECTOR3D_PTR v);

        /**
         * Compute the magnitude of a 3D vector
         * @param v
         * @return the magnitude of v
         */
        MATH_EXPORT real magnitudeh(const HVECTOR3D& v);

        /**
         * Compute the magnitude of a 3D vector
         * @param v
         * @return the magnitude of v
         */
        MATH_EXPORT real magnitudehEx(HVECTOR3D_PTR v);

#ifdef __cplusplus
        //quando si compila in c++ questa serve a fare in modo che non vengano 'decorati'
        //i nomi dal linker
    }
#endif

#ifdef _USENAMESPACE_

};

#endif

#endif
