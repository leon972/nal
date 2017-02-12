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
 
 \file 3dtransf.h
 \brief Geometric transformations with matrices (rotation,translation,scaling,euler angles).
 Global functions to perform basic geometric transformation with homogeneous matrices and vectors.
 
*/

#ifndef TRANSF_H_INCLUDED
#define TRANSF_H_INCLUDED

#include "libmath.h"
#include "vectors.h"  //modulo base per i vettori
#include "matrix.h"   //matrici
#include "quat.h"     //quaternions

#ifdef _USENAMESPACE_
namespace mathengine {
#endif

#ifdef __cplusplus
    //quando si compila in c++ questa serve a fare in modo che non vengano 'decorati'
    //i nomi dal linker
    extern "C" {
#endif

        /**
         * Rotate the vector v about a generic axis defined by versor vers by theta radians.
         * This function uses the Rodrigues formula to calculate the rotated vector.
         * @param v vectors to be rotate
         * @param vers axis of rotation versor
         * @param theta rotation angle in radians
         * @return the rotated vector
         */
        MATH_EXPORT HVECTOR3D rotateVector(const HVECTOR3D& v, const HVECTOR3D& vers, real theta);


        /**
         * Compute the longitude angle in radians.
         * @param vers vector
         * @return the angle between the plane containing vers and normal to plane X-Y and the plane X-Z
         */
        MATH_EXPORT real getLongitudeZ(const HVECTOR3D& vers);

        /**
         * Compute the azimut angle in radians.
         * @param vers vector
         * @return the angle between the vector vers and the plane X-Y
         */
        MATH_EXPORT real getAzimutZ(const HVECTOR3D& vers);

        /**
         * Create the counter-clockwise rotation matrix about Z axis.
         * The rows of the obtained matrix are the versors of the
         * rotated axes.
         * @param theta radians
         * @return the Z-axis rotation matrix
         */
        MATH_EXPORT HMATRIX33 getRotZMatrix(real theta);
        /**
         * Create the rotation matrix about Y axis.
         * @param theta radians
         * @return the Z-axis rotation matrix
         */
        MATH_EXPORT HMATRIX33 getRotYMatrix(real theta);
        /**
         * Create the rotation matrix about X axis.
         * @param theta radians
         * @return the Z-axis rotation matrix
         */
        MATH_EXPORT HMATRIX33 getRotXMatrix(real theta);

        //Ottiene la matrice che esegue la rotazione di theta radianti intorno al versore vers
        //(fomula di Rodrigues)
        /**
         * Create then rotation matrix about a generic axis.
         * (Rodrigues formula)
         * @param vers versor of the axis of rotation.
         * @param theta counter-clockwise angle of rotation in radians.
         * @return the rotation matrix.
         */
        MATH_EXPORT HMATRIX33 getRotMatrix(const HVECTOR3D& vers, real theta);

        /**
         * Copy the rotation matrix about a generic axis on an existing matrix
         * @param r destination matrix (this matrix becomes the rotation matrix,translation components are deleted)
         * @param vers versor of the axis of rotation.
         * @param theta counter-clockwise angle of rotation in radians.
         */
        MATH_EXPORT void getRotMatrixEx(HMATRIX33_PTR r, const HVECTOR3D& vers, real theta);
         
        /**
         * Create the Euler rotation matrix.* 
         * @param psi first angle of rotation about Z axis.
         * @param theta the angle of rotation about the node axis  (x')
         * @param phi rotation bout rotated Z' axis
         * @return the Euler marix.The rows are the versors of the rotated axes.
         */
        MATH_EXPORT HMATRIX33 getEulerMatrix(const real psi, const real theta, const real phi);

        /**
         * Copy the Euler matrix to m
         * @param m the destination matrix
         * @param psi first angle of rotation about Z axis.
         * @param theta the angle of rotation about the node axis  (x')
         * @param phi rotation bout rotated Z' axis
         */
        MATH_EXPORT void getEulerMatrixEx(HMATRIX33* m, const real psi, const real theta, const real phi);

        /**
         * Create the translation matrix.
         * The last column of the metrix contains tx,ty,tz
         * @param tx axis translation component
         * @param ty yaxis translation component
         * @param tz zaxis translation component
         * @return  the translation matrix.
         */
        HMATRIX33 getTranslationMatrix(const real tx, const real ty, const real tz);

        /**
         * Copy the translation matrix to an existing matrix.
         * Rotation components are erased and diagonal is set to 1
         * @param m the destination matrix
         * @param tx axis translation component
         * @param ty yaxis translation component
         * @param tz zaxis translation component
         */
        MATH_EXPORT void getTranslationMatrixEx(HMATRIX33* m, const real tx, const real ty, const real tz);

        /**
         * Create a scaling matrix.
         * The diagonal components a00,a11,a22 ar set to sx,sy,sz
         * @param sx x scale
         * @param sy y scale
         * @param sz z scale
         * @return The scaling matrix.
         */
        MATH_EXPORT HMATRIX33 getScaleMatrix(const real sx, const real sy, const real sz);

        /**
         * Copy the scaling matrix to m.
         * @param m destination matrix.
         * @param sx x scale
         * @param sy y scale
         * @param sz z scale
         */
        MATH_EXPORT void getScaleMatrixEx(HMATRIX33* m, const real sx, const real sy, const real sz);

#ifdef __cplusplus
        //quando si compila in c++ questa serve a fare in modo che non vengano 'decorati'
        //i nomi dal linker
    }
#endif

#ifdef _USENAMESPACE_
};
#endif

#endif // 3DTRANSF_H_INCLUDED
