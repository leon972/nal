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
 \file vectors3.h
 \brief 3D template-based homogeneous vector class.
 */

#ifndef VECTORS3_H_INCLUDED
#define VECTORS3_H_INCLUDED

#include "mathdefs.h"
#include <cmath>


namespace mathengine {

    /**
     * 3D Vector in homogeneous coordinates.
     * This is a 3d vector with 4 components
     * made for use with 4x4 transformation matrices.
     */
    template<class N> class CHVector3 {

private:

        void init(N x,N y,N z) {
            this->x=x;
            this->y=y;
            this->z=z;
            this->w=1;
        }

public:

        /**
         *x component.
         */
        N x;
        /**
         *y component.
         */
        N y;
        /**
         *z component.
         */
        N z;
        /**
         *homogeneous term, (usually set to 1)
         *This is a convenience component used to facilitate
         *multiplication between this vector and homogeneous transformation matrices.
         */
        N w; 

        /**
         * Create a vector.
         * The w component is set to 1.
         * @param x x component.
         * @param y y component.
         * @param z z component.
         */
        CHVector3(N x,N y,N z) {
            init(x,y,z); //coord. omogenee
        }
                
        /**
         *Create a zero vector.
         *x,y,z are set to zero, w is set to 1
         */
        CHVector3() {
            init(0,0,0);
        }

        /**
         * Gets the magnitude (length) of the vector.
         * @return the magnitude of this vector.
         */
        N magnitude() const {
            return sqrtl(x*x+y*y+z*z);
        }

        /**
         * Transforms this vector into its opposite vector.
         *
         */
        void negate() {
            x=-x;
            y=-y;
            z=-z;
        }

        /**
         * Returns the opposite vector of this vector.
         * @param res
         */
        void get_neg_vector(CHVector3<N>& res) const {
            res.x=-x;
            res.y=-y;
            res.z=-z;
        }

        /**
         * Transform this vector into its versor.
         * The versor has unit length and same direction.
         */
        void normalize() {
            N m=sqrtl(x*x+y*y+z*z);
            x/=m;
            y/=m;
            z/=m;
            w=1;
        }

        /**
         * Creates the versor of this vector and puts the result in res.
         * @param res the output parameter.
         */
        void get_normalized_vector(CHVector3<N>& res) const {
            N m=sqrt(x*x+y*y+z*z);
            res.x=x/m;
            res.y=y/m;
            res.z=z/m;
        }

        /**
         * Create the versor of this vector.
         * @return a unit length vector with the same direction of this vector.
         */
        CHVector3<N> get_normalized_vector()  const {

            CHVector3<N> res;
            N m=sqrt(x*x+y*y+z*z);
            res.x=x/m;
            res.y=y/m;
            res.z=z/m;
        }

        /**
         * Computes the dot product between this vector and v
         * @param v a vector
         * @return the dot product.
         */
        N dot_product(const CHVector3& v) {
            return x*v.x+y*v.y+z*v.z;
        }

        /**
         * Calculates the cross product between two vectors
         * and places the result into this vector.
         * After the call, this vector contains the cross product of the
         * two vectors v1 x v2
         * @param v1 the first vector
         * @param v2 the second vector
         */
        void cross_product(const CHVector3& v1,const CHVector3& v2) {
            x=(v1.y*v2.z-v1.z*v2.y);
            y=(v1.z*v2.x-v1.x*v2.z);
            z=(v1.x*v2.y-v1.y*v2.x);
            w=1;
        }

        /**
         * Calculates the cross product between this vector and v
         * and places the result into this vector.
         * @param v another vector.
         */
        void cross_product(const CHVector3& v) {
            N x1,y1,z1;

            x1=x;
            y1=y;
            z1=z;

            x=(y1*v.z-z1*v.y);
            y=(z1*v.x-x1*v.z);
            z=(x1*v.y-y1*v.x);
            w=1;
        }

        /**
         * Add to this vector another vector.
         * @param v the vector to be added.
         */
        void add(const CHVector3& v) {
            x+=v.x;
            y+=v.y;
            z+=v.z;
        }

        /**
         * Subtract a vector from this vector.
         * @param v the vector to be subtracted.
         */
        void sub(const CHVector3& v) {
            x-=v.x;
            y-=v.y;
            z-=v.z;
        }

        /**
         * Multiply this vector by a scalar.
         * @param s
         */
        void mul_scalar(N s) {
            x*=s;
            y*=s;
            z*=s;
        }

        /**
         * Set the components of this vector.
         * @param x
         * @param y
         * @param z
         */
        void set(N x,N y,N z) {
            this->x=x;
            this->y=y;
            this->z=z;
            w=1;
        }

        /**
         * Rotate this vector around an axis.
         * This method uses the Rodrigues' formula to accomplish the rotation.
         * @param vers versor of the axis of rotation (must have unit length!).
         * @param theta the angle of rotation in radians.
         */

        void rotate(const CHVector3& vers,const float theta) {

            double ct=::cos(theta);
            double st=::sin(theta);

            //x= prodotto vettoriale
            //v*cos(theta)+vers x v*sin(theta)+(vers,v)*(1-cos(theta)) vers
            CHVector3<N> v1,v2,v3;
            
            v3=*this;
            v1=*this;

            //v*cos(theta)
            v1.mul_scalar(ct);

            //vers x v*sin(theta)
            v3.mul_scalar(st);
            v2=vers;
            v2.cross_product(v3);

            //(vers,v)*(1-cos(theta)) vers
            v3=vers;
            v3.mul_scalar(this->dot_product(vers)*(1-ct));

            x=v1.x;
            y=v1.y;
            z=v1.z;

            add(v2);
            add(v3);

            w=1;
        }

        /**
         * Calculates the angle between the plane containing this vector
         * and normal to the x-y plane and the x-z plane.
         * @return the angle in radians.
         */
        N getLongitudeZ() {
            return acos(x/sqrt(x*x+y*y));
        }

        /**
         * Calculates the angle between this vector and the plane x-y.
         * @return the angle in radians.
         */
        N getAzimutZ() {
            N len=x*x+y*y+z*z;
            return acos(sqrt((x*x+y*y)/len));
        }
    };

    /**
     * Vector addition.
     * @param result the sum v1+v2
     * @param v1 first vector.
     * @param v2 second vector.
     */
    template<class N> void add_vectors(CHVector3<N>* result,const CHVector3<N>& v1,const CHVector3<N>& v2) {

        result->x=v1.x+v2.x;
        result->y=v1.y+v2.y;
        result->z=v1.z+v2.z;
        result->w=1;

    }

    /**
     * Vector addition v1+v2.
     * @param v1 first vector.
     * @param v2 second vector.
     * @return the sum of v1 and v2.
     */
    template<class N> CHVector3<N> add_vectors(const CHVector3<N>& v1,const CHVector3<N>& v2) {

        CHVector3<N> r;

        r.x=v1.x+v2.x;
        r.y=v1.y+v2.y;
        r.z=v1.z+v2.z;
        r.w=1;

        return r;
    }

    /**
     * Vector addition operator.
     * @param v1 first vector.
     * @param v2 second vector.
     * @return v1+v2.
     */
    template<class N> CHVector3<N> operator+(const CHVector3<N>& v1,const CHVector3<N>& v2) {

        CHVector3<N> r;

        r.x=v1.x+v2.x;
        r.y=v1.y+v2.y;
        r.z=v1.z+v2.z;
        r.w=1;

        return r;
    }


    /**
     * Vector subtraction.
     * @param result v1-v2
     * @param v1 first vector.
     * @param v2 second vector.
     */
    template<class N> void sub_vectors(CHVector3<N>* result,const CHVector3<N>& v1,const CHVector3<N>& v2) {

        result->x=v1.x-v2.x;
        result->y=v1.y-v2.y;
        result->z=v1.z-v2.z;
        result->w=1;

    }

    /**
     * Vector subtraction v1-v2.
     * @param v1 first vector.
     * @param v2 second vector.
     * @return v1-v2
     */
    template<class N> CHVector3<N> sub_vectors(const CHVector3<N>& v1,const CHVector3<N>& v2) {

        CHVector3<N> r;

        r.x=v1.x-v2.x;
        r.y=v1.y-v2.y;
        r.z=v1.z-v2.z;
        r.w=1;

        return r;

    }

    /**
     * Vector subtraction operator.
     * @param v1 first vector.
     * @param v2 second vector.
     * @return v1-v2.
     */
    template<class N> CHVector3<N> operator-(const CHVector3<N>& v1,const CHVector3<N>& v2) {

        CHVector3<N> r;

        r.x=v1.x-v2.x;
        r.y=v1.y-v2.y;
        r.z=v1.z-v2.z;
        r.w=1;

        return r;

    }

    /**
     * Calculates the cross production v1 x v2.
     * @param result the cross product between v1 and v2.
     * @param v1 first vector.
     * @param v2 second vector.
     */
    template<class N> void cross_product(CHVector3<N>* result,const CHVector3<N>& v1,const CHVector3<N>& v2) {
        result->x=(v1.y*v2.z-v1.z*v2.y);
        result->y=(v1.z*v2.x-v1.x*v2.z);
        result->z=(v1.x*v2.y-v1.y*v2.x);
        result->w=1;
    }

    /**
     * Calculates the cross product v1 x v2.
     * @param v1
     * @param v2
     * @return v1 x 2.
     */
    template<class N> CHVector3<N> cross_product(const CHVector3<N>& v1,const CHVector3<N>& v2) {
        CHVector3<N> r;

        r.x=(v1.y*v2.z-v1.z*v2.y);
        r.y=(v1.z*v2.x-v1.x*v2.z);
        r.z=(v1.x*v2.y-v1.y*v2.x);
        r.w=1;

        return r;
    }

    /**
     * Cross product
     * @param v1
     * @param v2
     * @return v1 x v2
     */
    template<class N> CHVector3<N> operator * (const CHVector3<N>& v1,const CHVector3<N>& v2) {

        CHVector3<N> r;

        r.x=(v1.y*v2.z-v1.z*v2.y);
        r.y=(v1.z*v2.x-v1.x*v2.z);
        r.z=(v1.x*v2.y-v1.y*v2.x);
        r.w=1;

        return r;
    }

    /**
     * Vector by scalar
     * @param v
     * @param scalar
     * @return the vector v multiplied by the scalar scalar
     */
   template<class N> CHVector3<N> operator * (const CHVector3<N>& v,const N scalar) {

        CHVector3<N> r(v);
        r.mul_scalar(scalar);
        return r;
    }

    /**
     * Scalar by vector
     * @param scalar
     * @param v
     * @return the vector v multiplied by the scalar scalar
     */
    template<class N> CHVector3<N> operator * (const N scalar,const CHVector3<N>& v) {


       
        CHVector3<N> r(v);
        r.mul_scalar(scalar);
        return r;
    }

    /**
     *Dot product.
     */
    template<class N> N dot_product  (const CHVector3<N>& v1,const CHVector3<N>& v2) {

         return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
    }

    
    /**
     * Calculates the versor of a vector.
     * @param v a vector.
     * @return the versor of v.
     */
    template<class N> CHVector3<N> norm(const CHVector3<N>& v)
    {
        N d=sqrtl(v.x*v.x+v.y*v.y+v.z*v.z);
        CHVector3<N> r;
        r.x=v.x/d;
        r.y=v.y/d;
        r.z=v.z/d;
        r.w=1;

        return r;
    }

    /**
     * Creates a vector given its components.
     * @param x
     * @param y
     * @param z
     * @return the vector (x,y,z,1)
     */
    template<class N> CHVector3<N> create_vector(N x,N y,N z)
    {
        CHVector3<N> v;
        v.x=x;
        v.y=y;
        v.z=z;

        return v;
    }

};

#endif 
