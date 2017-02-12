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

/****************************************************************
 Code by Leonardo Berti (C) 2006
 Math engine - vettori e matrici 2D 3D e 3D in coordinate omogenee

 ****************************************************************/

#include "math/vectors.h"
#include <math.h>
#include <mem.h>

#ifdef _USENAMESPACE_

namespace mathengine {

#endif

    //costruzione dei vettori

    VECTOR2D getVector2D(real x, real y) {
        auto VECTOR2D r;
        r.x = x;
        r.y = y;

        return (r); //costruisce il vettore sullo stack
    }

    VECTOR3D getVector3D(real x, real y, real z) {
        auto VECTOR3D r;

        r.x = x;
        r.y = y;
        r.z = z;

        return (r);

    }

    HVECTOR3D getHVector3D(real x, real y, real z) {
        auto HVECTOR3D r;

        r.x = x;
        r.y = y;
        r.z = z;
        r.w = 1;

        return (r);
    }

    //output

    void printVector2D(FILE* out, const VECTOR2D& v) {
        fprintf(out, "(%f,%f)", v.x, v.y);
    }

    void printVector3D(FILE* out, const VECTOR3D& v) {
        fprintf(out, "(%f,%f,%f)", v.x, v.y, v.z);
    }

    void printVectorH3D(FILE* out, const HVECTOR3D& v) {
        fprintf(out, "(%f,%f,%f)", v.x, v.y, v.z, v.w);
    }

    //creazione

    void getVector2D(VECTOR2D_PTR result, real x, real y) {
        result->x = x;
        result->y = y;
    }

    void getVector3D(VECTOR3D_PTR result, real x, real y, real z) {
        result->x = x;
        result->y = y;
        result->z = z;

    }

    void getHVector3D(HVECTOR3D_PTR result, real x, real y, real z) {
        result->x = x;
        result->y = y;
        result->z = z;
        result->w = 1;

    }

    void zeroVector2(VECTOR2D &vec2) {
        vec2.x = 0;
        vec2.y = 0;
    }

    void zeroVector3(VECTOR3D &vec3) {

        vec3.x = 0;
        vec3.y = 0;
        vec3.z = 0;

    }

    void zeroVector3H(HVECTOR3D &vec3) {

        vec3.x = 0;
        vec3.y = 0;
        vec3.z = 0;
        vec3.w = 1;

    }

    //somma di vettori v1+v2

    void sumVector2DEx(VECTOR2D_PTR result, const VECTOR2D& v1, const VECTOR2D& v2) {
        result->x = v1.x + v2.x;
        result->y = v1.y + v2.y;
    }

    VECTOR2D sumVector2D(const VECTOR2D& v1, const VECTOR2D& v2) {
        VECTOR2D r;

        r.x = v1.x + v2.x;
        r.y = v1.y + v2.y;

        return (r);
    }

    //identit�

    void getUnitVector2D(VECTOR2D_PTR result) {
        result->x = result->y = 1;
    }

    void getUnitVector3D(VECTOR3D_PTR result) {
        result->x = result->y = result->z = 1;
    }

    void getUnitHVector3D(HVECTOR3D_PTR result) {
        result->x = result->y = result->z = result->w = 1;
    }

    void sumVector3DEx(VECTOR3D_PTR result, const VECTOR3D& v1, const VECTOR3D& v2) {
        result->x = v1.x + v2.x;
        result->y = v1.y + v2.y;
        result->z = v1.z + v2.z;
    }

    VECTOR3D sumVector3D(const VECTOR3D& v1, const VECTOR3D& v2) {
        auto VECTOR3D r;

        r.x = v1.x + v2.x;
        r.y = v1.y + v2.y;
        r.z = v1.z + v2.z;

        return (r);

    }

    //somma fra vettori in coordinate omogenee

    void sumHVector3DEx(HVECTOR3D_PTR result, const HVECTOR3D& v1, const HVECTOR3D& v2) {
        result->x = v1.x + v2.x;
        result->y = v1.y + v2.y;
        result->z = v1.z + v2.z;

    }

    HVECTOR3D sumHVector3D(const HVECTOR3D& v1, const HVECTOR3D& v2) {
        auto HVECTOR3D r;

        r.x = v1.x + v2.x;
        r.y = v1.y + v2.y;
        r.z = v1.z + v2.z;

        return (r);
    }

    /**ottiene il vettore opposto*/
    VECTOR2D getNegVector2D(const VECTOR2D& v) {
        auto VECTOR2D r;
        r.x = -v.x;
        r.y = -v.y;
        return r;
    }

    VECTOR3D getNegVector3D(const VECTOR3D& v) {
        auto VECTOR3D r;

        r.x = -v.x;
        r.y = -v.y;
        r.z = -v.z;

        return r;
    }

    HVECTOR3D getNegHVector3D(const HVECTOR3D& v) {
        HVECTOR3D r;

        r.x = -v.x;
        r.y = -v.y;
        r.z = -v.z;
        r.w = 1;

        return r;
    }

    void negVector2D(VECTOR2D_PTR v) {
        v->x = -v->x;
        v->y = -v->y;
    }

    void negVector3D(VECTOR3D_PTR v) {
        v->x = -v->x;
        v->y = -v->y;
        v->z = -v->z;
    }

    void negHVector3D(HVECTOR3D_PTR v) {
        v->x = -v->x;
        v->y = -v->y;
        v->z = -v->z;
    }

    /**aggiunge a dest il vettore v dest+=v */
    void addVector2D(VECTOR2D_PTR dest, const VECTOR2D& v) {
        dest->x += v.x;
        dest->y += v.y;
    }

    void addVector3D(VECTOR3D_PTR dest, const VECTOR3D& v) {
        dest->x += v.x;
        dest->y += v.y;
        dest->z += v.z;
    }

    void addHVector3D(HVECTOR3D_PTR dest, const HVECTOR3D& v) {
        dest->x += v.x;
        dest->y += v.y;
        dest->z += v.z;
    }   
   

    //differenza fra vettori v1-v2

    void subVector2DEx(VECTOR2D_PTR result, const VECTOR2D& v1, const VECTOR2D& v2) {
        result->x = v1.x - v2.x;
        result->y = v1.y - v2.y;
    }

    VECTOR2D subVector2D(const VECTOR2D& v1, const VECTOR2D& v2) {
        VECTOR2D r;

        r.x = v1.x - v2.x;
        r.y = v1.y - v2.y;

        return (r);
    }

    void subVector3DEx(VECTOR3D_PTR result, const VECTOR3D& v1, const VECTOR3D& v2) {
        result->x = v1.x - v2.x;
        result->y = v1.y - v2.y;
        result->z = v1.z - v2.z;

    }

    VECTOR3D subVector3D(const VECTOR3D& v1, const VECTOR3D& v2) {
        VECTOR3D r;

        r.x = v1.x - v2.x;
        r.y = v1.y - v2.y;
        r.z = v1.z - v2.z;

        return (r);

    }

    void subHVector3DEx(HVECTOR3D_PTR result, const HVECTOR3D& v1, const HVECTOR3D& v2) {
        result->x = v1.x - v2.x;
        result->y = v1.y - v2.y;
        result->z = v1.z - v2.z;
    }

    HVECTOR3D subHVector3D(const HVECTOR3D& v1, const HVECTOR3D& v2) {
        HVECTOR3D r;

        r.x = v1.x - v2.x;
        r.y = v1.y - v2.y;
        r.z = v1.z - v2.z;

        return (r);
    }

    //prodotto fra uno scalare e un vettore

    void mulHVector3DEx(HVECTOR3D_PTR result, const HVECTOR3D&v, const real scalar) {
        result->x = scalar * v.x;
        result->y = scalar * v.y;
        result->z = scalar * v.z;
    }

    HVECTOR3D mulHVector3D(const HVECTOR3D&v, const real scalar) {
        auto HVECTOR3D r;

        r.x = scalar * v.x;
        r.y = scalar * v.y;
        r.z = scalar * v.z;

        return r;
    }

    //prodotto scalare

    real dotProduct2(const VECTOR2D& v1, const VECTOR2D& v2) {
        return v1.x * v2.x + v1.y * v2.y;
    }

    real dotProduct3(const VECTOR3D& v1, const VECTOR3D& v2) {
        return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
    }

    real dotProduct3H(const HVECTOR3D& v1, const HVECTOR3D& v2) {
        return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
    }

    //prodotto vettoriale (v1 X v2)

    VECTOR3D crossProduct3(const VECTOR3D& v1, const VECTOR3D& v2) {
        auto VECTOR3D r; //variabile creata sollo stack (memoria rilasciata automaticamente)

        r.x = (v1.y * v2.z - v1.z * v2.y);
        r.y = (v1.z * v2.x - v1.x * v2.z);
        r.z = (v1.x * v2.y - v1.y * v2.x);

        return (r);
    }

    void crossProduct3Ex(VECTOR3D_PTR result, const VECTOR3D& v1, const VECTOR3D& v2) {
        result->x = (v1.y * v2.z - v1.z * v2.y);
        result->y = (v1.z * v2.x - v1.x * v2.z);
        result->z = (v1.x * v2.y - v1.y * v2.x);
    }

    HVECTOR3D crossProduct3H(const HVECTOR3D& v1, const HVECTOR3D& v2) {
        auto HVECTOR3D r;

        r.x = (v1.y * v2.z - v1.z * v2.y);
        r.y = (v1.z * v2.x - v1.x * v2.z);
        r.z = (v1.x * v2.y - v1.y * v2.x);
        r.w = 1;

        return (r);

    }

    //ottiene il versore del prodotto vettoriale tra due vettori

    HVECTOR3D normCrossProduct3H(const HVECTOR3D& v1, const HVECTOR3D& v2) {
        auto HVECTOR3D r;

        r.x = (v1.y * v2.z - v1.z * v2.y);
        r.y = (v1.z * v2.x - v1.x * v2.z);
        r.z = (v1.x * v2.y - v1.y * v2.x);
        r.w = 1;

        real m = 1 / sqrtr(r.x * r.x + r.y * r.y + r.z * r.z);

        r.x *= m;
        r.y *= m;
        r.z *= m;

        return (r);

    }

    void normCrossProduct3HEx(HVECTOR3D_PTR result, const HVECTOR3D& v1, const HVECTOR3D& v2) {
        result->x = (v1.y * v2.z - v1.z * v2.y);
        result->y = (v1.z * v2.x - v1.x * v2.z);
        result->z = (v1.x * v2.y - v1.y * v2.x);
        result->w = 1;

        real m = 1 / sqrtr(result->x * result->x + result->y * result->y + result->z * result->z);

        result->x *= m;
        result->y *= m;
        result->z *= m;
    }

    VECTOR3D normCrossProduct3(const VECTOR3D& v1, const VECTOR3D& v2) {
        auto VECTOR3D r;

        r.x = (v1.y * v2.z - v1.z * v2.y);
        r.y = (v1.z * v2.x - v1.x * v2.z);
        r.z = (v1.x * v2.y - v1.y * v2.x);

        real m = 1 / sqrtr(r.x * r.x + r.y * r.y + r.z * r.z);

        r.x *= m;
        r.y *= m;
        r.z *= m;

        return (r);

    }

    void normCrossProduct3Ex(VECTOR3D_PTR result, const VECTOR3D& v1, const VECTOR3D& v2) {
        result->x = (v1.y * v2.z - v1.z * v2.y);
        result->y = (v1.z * v2.x - v1.x * v2.z);
        result->z = (v1.x * v2.y - v1.y * v2.x);

        real m = 1 / sqrtr(result->x * result->x + result->y * result->y + result->z * result->z);

        result->x *= m;
        result->y *= m;
        result->z *= m;
    }

    //ottiene il versore di un vettore 2D

    VECTOR2D normalizeVector2D(const VECTOR2D& v) {
        auto VECTOR2D r;
        real m = 1 / sqrtr(v.x * v.x + v.y * v.y);

        r.x = v.x*m;
        r.y = v.y*m;

        return r;

    }

    void normalizeVector2DEx(VECTOR2D_PTR v) {

        real m = 1 / sqrtr(v->x * v->x + v->y * v->y);
        v->x *= m;
        v->y *= m;
    }

    VECTOR3D normalizeVector3D(const VECTOR3D& v) {
        auto VECTOR3D r;

        real m = 1 / sqrtr(v.x * v.x + v.y * v.y + v.z * v.z);

        r.x = v.x*m;
        r.y = v.y*m;
        r.z = v.z*m;

        return r;
    }

    void normalizeVector3DEx(VECTOR3D_PTR v) {
        real m = 1 / sqrtr(v->x * v->x + v->y * v->y + v->z * v->z);

        v->x *= m;
        v->y *= m;
        v->z *= m;
    }

    //ottiene un versore a partire da un vettore 3D in coordinate omogenee

    HVECTOR3D normalizeHVector3D(const HVECTOR3D& v) {
        real m = 1 / sqrtr(v.x * v.x + v.y * v.y + v.z * v.z);

        HVECTOR3D r;

        r.x = v.x*m;
        r.y = v.y*m;
        r.z = v.z*m;

        return r;
    }

    void normalizeHVector3DEx(HVECTOR3D_PTR v) {
        real m = 1 / sqrtr(v->x * v->x + v->y * v->y + v->z * v->z);

        v->x *= m;
        v->y *= m;
        v->z *= m;
    }

    //modulo

    real magnitude2(const VECTOR2D& v) {
        return sqrtr(v.x * v.x + v.y * v.y);
    }

    real magniture2Ex(VECTOR2D_PTR v) {
        return sqrtr(v->x * v->x + v->y * v->y);
    }

    real magnitude(const VECTOR3D& v) {
        return sqrtr(v.x * v.x + v.y * v.y + v.z * v.z);
    }

    real magnitudeEx(VECTOR3D_PTR v) {
        return sqrtr(v->x * v->x + v->y * v->y + v->z * v->z);
    }

    real magnitudeh(const HVECTOR3D& v) {
        return sqrtr(v.x * v.x + v.y * v.y + v.z * v.z);
    }

    real magnitudehEx(HVECTOR3D_PTR v) {
        return sqrtr(v->x * v->x + v->y * v->y + v->z * v->z);
    }

    void crossProduct3HEx(HVECTOR3D_PTR result, HVECTOR3D_PTR v1, HVECTOR3D_PTR v2) {
        result->x = (v1->y * v2->z - v1->z * v2->y);
        result->y = (v1->z * v2->x - v1->x * v2->z);
        result->z = (v1->x * v2->y - v1->y * v2->x);
        result->w = 1;
    }

#ifdef _USENAMESPACE_

}

#endif
