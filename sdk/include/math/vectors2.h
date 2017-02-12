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

#ifndef VECTORS2_H_INCLUDED
#define VECTORS2_H_INCLUDED
/*!
 \file vectors2.h
 \brief 2D template-based vector class.
 */

#include "mathdefs.h"
#include <cmath>

#ifdef _USENAMESPACE_
namespace mathengine {
#endif

    using namespace std;

    /**
     *2 Components vector.
     */
    template<class N> class CVector2 {

private:

        void init(N x,N y) {
            this->x=x;
            this->y=y;
        }

public:

        N x,y;

        /**
         * Construct a 2D vector
         * @param x
         * @param y
         */
        CVector2(N x,N y) {
            init(x,y);
        }

        /**
         * Construct a zero vector.
         */
        CVector2() {
            init(0,0);
        }

        /**
         * Get the magnitude (length) of vector
         * @return the magnitude.
         */
        N magnitude() const {
            return sqrt(x*x+y*y);
        }

        /**
         * Negate the vector
         */
        void negate() {
            x=-x;
            y=-y;
        }

        /**
         * Compute the opposite vector
         * @param res destination vector where the opposite vector id copied.
         */
        void get_neg_vector(CVector2<N>& res) const {
            res.x=-x;
            res.y=-y;

        }

        /**
         * Transform this vector into a versor (vector with magnitude 1)
         */
        void normalize() {
            N m=sqrt(x*x+y*y);
            x/=m;
            y/=m;
        }

        /**
         * Computes the versor of this vector
         * @param res return the versor of this vector
         */
        void get_normalized_vector(CVector2<N>& res) const {
            N m=sqrt(x*x+y*y);
            res.x=x/m;
            res.y=y/m;

        }

        /**
         * Dot product between this vector and v
         * @param v second term of the dot product
         * @return the dot product
         */
        N dot_product(const CVector2& v) {
            return x*v.x+y*v.y;
        }

        /**
         * Adds a vector to this vector.
         * @param v the vector to be added to this vector.
         */
        void add(const CVector2& v) {
            x+=v.x;
            y+=v.y;
        }

        /**
         * subtracts a vector from this vector.
         * @param v
         */
        void sub(const CVector2& v) {
            x-=v.x;
            y-=v.y;
        }

        /**
         * Vector by scalar multiplication.
         * @param s
         */
        void mul_scalar(N s) {
            x*=s;
            y*=s;
        }

        /**
         * Sets the components of this vector.
         * @param x the x component
         * @param y the y component
         *
         */
        void set(N x,N y) {
            this->x=x;
            this->y=y;
        }

        /**
         * Rotates this vector counter-clockwise.
         * @param theta angle of rotation in radians.
         */
        void rotate(const float theta) {

            double ct=::cos(theta);
            double st=::sin(theta);

            N x1=x;
            N y1=y;

            x=x1*ct-y1*st;
            y=x1*st+y1*ct;

        }

        /**
         * Gets the angle between this vector and the x axis.
         * @return the angle in radians.
         */
        N getTheta() const {
            return acos(x/sqrt(x*x+y*y));
        }

    };

    /**
     * Vectors addition.
     * @param result v1+v2
     * @param v1 first vector
     * @param v2 second vector
     */
    template<class N> void add_vectors(CVector2<N>* result,const CVector2<N>& v1,const CVector2<N>& v2) {
        result->x=v1.x+v2.x;
        result->y=v1.y+v2.y;

    }

    /**
     * Vector addition,
     * @param v1 first vector.
     * @param v2 secnd vector.
     * @return v1+v1.
     */
    template<class N> CVector2<N> add_vectors(const CVector2<N>& v1,const CVector2<N>& v2) {

        CVector2<N> r;

        r.x=v1.x+v2.x;
        r.y=v1.y+v2.y;

        return r;
    }

    /**
     * Vectors addition.
     * @param v1 first vector
     * @param v2 second vector
     * @return v1+v2
     */
    template<class N> CVector2<N> operator+(const CVector2<N>& v1,const CVector2<N>& v2) {
        CVector2<N> r;

        r.x=v1.x+v2.x;
        r.y=v1.y+v2.y;

        return r;
    }


    /**
     * Vectors subtraction.
     * @param result v1-v2
     * @param v1 first vector
     * @param v2 second vector
     */
    template<class N> void sub_vectors(CVector2<N>* result,const CVector2<N>& v1,const CVector2<N>& v2) {
        result->x=v1.x-v2.x;
        result->y=v1.y-v2.y;

    }

    /**
     * Vectors subtraction.
     * @param v1 first vector
     * @param v2 second vector
     * @return v1-v2
     */
    template<class N> CVector2<N> sub_vectors(const CVector2<N>& v1,const CVector2<N>& v2) {
        CVector2<N> r;

        r.x=v1.x-v2.x;
        r.y=v1.y-v2.y;


        return r;
    }

    /**
     * Vectors subtraction.
     * @param v1 first vector
     * @param v2 second vector
     * @return v1-v2
     */
    template<class N> CVector2<N> operator-(const CVector2<N>& v1,const CVector2<N>& v2) {
        CVector2<N> r;

        r.x=v1.x-v2.x;
        r.y=v1.y-v2.y;


        return r;
    }


    /**
     * Get the versor of a vector.
     * @param v vector.
     * @return the versor of v (the versor has the same direction and orientation of v and unit length)
     */
    template<class N> CVector2<N> norm(const CVector2<N>& v)
    {
        N d=sqrt(v.x*v.x+v.y*v.y);
        CVector2<N> r;
        r.x=v.x/d;
        r.y=v.y/d;


        return r;
    }

    /**
     * Create a vector.
     * @param x x component
     * @param y y component
     * @return the vector (x,y)
     */
    template<class N> CVector2<N> create_vector(N x,N y)
    {
        CVector2<N> v;
        v.x=x;
        v.y=y;
        return v;
    }


#ifdef _USENAMESPACE_
};
#endif

#endif // VECTORS2_H_INCLUDED



