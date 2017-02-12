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
 \file quat2.h
 \brief Template-based quaternion class.

 */
#ifndef QUAT2_H_INCLUDED
#define QUAT2_H_INCLUDED

#include <cmath>
#include "vectors3.h"
#include "matrix2.h"

namespace mathengine {

    /**
     * pi / 2
     */
    const double pi_2 = 1.57079632679489661923;

    /**
     * Quaternion class.
     * The quaternion has a real part r and an imaginary part represented
     * by a vector (x,y,z)
     */
    template<class N> class CQuaternion {
    public:

        /**
         Real part
         */
        N r;
        /**
         *Imaginary x
         */
        N x;
        /**
         *Imaginary y
         */
        N y;
        /**
         *Imaginary z
         */
        N z;
        

        /**
         * Sets real and imaginary part
         * @param r real part
         * @param x imaginary component x
         * @param y imaginary component y
         * @param z imaginary component z
         */
        void set(N r, N x, N y, N z) {
            this->r = r;
            this->x = x;
            this->y = y;
            this->z = z;
        }

        /**
         * Create a zero quaternion
         */
        CQuaternion() {
            set(0, 0, 0, 0);
        }

        /**
         * create a quaternion
         * @param r real part
         * @param x imaginary component x
         * @param y imaginary component y
         * @param z imaginary component z
         */
        CQuaternion(N r, N x, N y, N z) {
            set(r, x, y, z);
        }

        /**
         * Assignment operator.
         *
         * @param q
         * @return this quaternion
         */
        CQuaternion<N>& operator=(const CQuaternion<N>& q) {

            x=q.x;
            y=q.y;
            z=q.z;
            r=q.r;

            return *this;

        }

        /**
         * Set all components to v (used to set to zero)
         * @param v
         * @return
         */
        CQuaternion<N>& operator=(const N &v) {

            x=v;
            y=v;
            z=v;
            r=v;

            return *this;
        }

        /**
         * Copies the imaginary part from a vector
         * @param v imaginary part
         */
        void from_vector(const CHVector3<N>& v) {
            x = v.x;
            y = v.y;
            z = v.z;
            r = 0;
        }

        /**
         * Copies the imaginary part into a vector
         * @param v
         */
        void to_vector(CHVector3<N>* v) {
            v->x = x;
            v->y = y;
            v->z = z;
        }

        /**
         * Returns the imaginary part as a vector
         * @return
         */
        CHVector3<N> to_vector() {
            CHVector3<N> v;
            to_vector(&v);
            return v;
        }

        /**
         * Add to this quaternion another quaternion
         * @param q quaternion to be added
         */
        void add(const CQuaternion<N>& q) {
            x += q.x;
            y += q.y;
            z += q.z;
            r += q.r;
        }

        /**
         * Addition operator
         * @param q
         * @return the sum between this quaternion and q
         */
        CQuaternion<N> operator +(const CQuaternion<N> &q) {

            CQuaternion<N> res;

            res=*this;

            res.add(q);

            return res;
        }

        /**
         * Addition and assignment operator
         * @param q
         * @return this quaternion after addition with q
         */
        CQuaternion<N>& operator +=(const CQuaternion<N> &q) {
        
            add(q);
            return *this;
        }

        /**
         * Subtract from this quaternion another quaternion
         * @param q
         */
        void sub(const CQuaternion<N>& q) {
            x -= q.x;
            y -= q.y;
            z -= q.z;
            r -= q.r;
        }

        /**
         * subtraction operator
         * @param q
         * @return
         */
        CQuaternion<N> operator -(const CQuaternion<N> &q) {

            CQuaternion<N> res;

            res.x=x-q.x;
            res.y=y-q.y;
            res.z=z-q.z;
            res.r=r-q.r;

            return res;
        }

        /**
         * subtraction and assignment operator.
         * @return this quaternion after the subtraction with q
         */
        CQuaternion<N>& operator -= (const CQuaternion<N> & q) {

            sub(q);
            return *this;

        }

        /**
         * Normalise this quaternion
         */
        void normalize() {
            N d = 1 / sqrt(r * r + x * x + y * y + z * z);

            r = r*d;
            x = x*d;
            y = y*d;
            z = z*d;
        }

        /**
         * Get the magnitude of this quaternion
         * @return
         */
        N magnitude() {
            return std::sqrt(r * r + x * x + y * y + z * z);
        }

        /**
         * Transform this quaternion to its conjugate
         */
        void conjugate() {
            x = -x;
            y = -y;
            z = -z;
        }

        /**
         * Multiply this quaternion by another quaternion
         * @param q
         */
        void mul(const CQuaternion<N>& q) {

            //riduce il numero di moltiplicazioni ragruppando
            N prd_0 = (z - y) * (q.y - q.z);
            N prd_1 = (r + x) * (q.r + q.x);
            N prd_2 = (r - x) * (q.y + q.z);
            N prd_3 = (y + z) * (q.r - q.x);
            N prd_4 = (z - x) * (q.x - q.y);
            N prd_5 = (z + x) * (q.x + q.y);
            N prd_6 = (r + y) * (q.r - q.z);
            N prd_7 = (r - y) * (q.r + q.z);
            N prd_8 = prd_5 + prd_6 + prd_7;
            N prd_9 = 0.5 * (prd_4 + prd_8);

            r = prd_0 + prd_9 - prd_5;
            x = prd_1 + prd_9 - prd_8;
            y = prd_2 + prd_9 - prd_7;
            z = prd_3 + prd_9 - prd_6;
        }

        /**
         * Multiplication operator
         * @param q another quaternion
         * @return the multiplication between this quaternion and q
         */
        CQuaternion<N> operator * (const CQuaternion<N> &q) {

            CQuaternion<N> qres;
            qres=*this;
            qres.mul(q);
            return qres;
        }

        /**
         * Multiplication and assignment operator
         */
        CQuaternion<N>& operator *= (const CQuaternion<N> &q) {

            mul(q);
            return *this;
        }

        /**
         * Transform this quaternion to a quaternion
         * representing a rotation about and axis.
         * @param versAxis versor of the axis of rotation (has unit length)
         * @param theta angle of rotation in radians
         */
        void set_rotation(const CHVector3<N>& versAxis, N theta) {

            N theta_div_2 = 0.5 * theta;

            N st2 = sin(theta_div_2);

            x = st2 * versAxis.x;
            y = st2 * versAxis.y;
            z = st2 * versAxis.z;
            r = cos(theta_div_2);
        }

        /**
         * Transform this quaternion to the triple product q1*q2*q3
         *  \note The triple product can be used to rotate a vector about any axis.
         * If the quaternion q defines the axis of rotation and rotation angle (@see getRotQuat3D)
         * and v the vector to be rotated, then the imaginary part of q * v * q' represents the rotated vector
         * v' about the axis. q' is the conjugate of q.
         * @param q1 quaternion representing the axis of rotation and rotation
         * @param q2 quaternion representi the vector to be rotated (imaginary part of q2)
         * @param q3 conjugate of q1 (for counter clockwise rotation,for clockwise rotation switch q1 and q3)
         */
        void triple_product(const CQuaternion<N>& q1, const CQuaternion<N>& q2, const CQuaternion<N>& q3) {

            *this = q1;
            mul(q2);
            mul(q3);
        }

        /**
         * If this quaternion is a rotation quaternion, copies in m
         * the rotation matrix corresponding to this quaternion.
         * @param m
         */
        void to_matrix(CHMatrix33<N>* m) {

            N el[16];

            el[3] = el[7] = el[11] = 0;
            el[15] = 1;

            N xy = x*y;
            N xz = x*z;
            N yz = y*z;
            N wx = x*r;
            N wy = y*r;
            N wz = z*r;

            N x2 = x*x;
            N y2 = y*y;
            N z2 = z*z;

            el[0] = 1 - 2 * y2 - 2 * z2;
            el[4] = 2 * xy - 2 * wz;
            el[8] = 2 * xz + 2 * wy;
            el[12] = 0;

            el[1] = 2 * xy + 2 * wz;
            el[5] = 1 - 2 * x2 - 2 * z2;
            el[9] = 2 * yz - 2 * wx;
            el[13] = 0;

            el[2] = 2 * xz - 2 * wy;
            el[6] = 2 * yz + 2 * wx;
            el[10] = 1 - 2 * x2 - 2 * y2;
            el[14] = 0;

            m->set_elements(el);
        }

        /**
         * If this quaternion is a rotation quaternion, copies in m
         * the rotation matrix corresponding to this quaternion.
         * @return
         */
        CHMatrix33<N> to_matrix() {
            CHMatrix33<N> m;
            to_matrix(&m);
            return m;
        }

        /**
         * Transforms this quaternion to the rotation quaternion
         * equivalent to the rotation matrix matrix
         * @param matrix rotation matrix
         */
        void from_matrix(const CHMatrix33<N>& matrix) {

            N m[16];

            matrix.copy_elements(m);

            N T = 1.0 + m[0] + m[5] + m[10]; //1+ diagonale
            N S;

            if (T > EPSILON) { //esegue questo controllo per evitare distorsioni

                S = sqrt(T)*2;
                r = 0.25 * S;
                x = (m[6] - m[9]) / S;
                y = (m[8] - m[2]) / S;
                z = (m[1] - m[4]) / S;

            } else {
                //determina l'elemento diagonale maggiore
                if (m[0] > m[5] && m[0] > m[10]) {

                    S = sqrt(1.0 + m[0] - m[5] - m[10])*2;
                    x = 0.25 * S;
                    y = (m[1] + m[4]) / S;
                    z = (m[8] + m[2]) / S;
                    r = (m[6] - m[9]) / S;

                } else if (m[5] > m[10]) {

                    S = sqrt(1.0 + m[5] - m[0] - m[10]) * 2;
                    x = (m[1] + m[4]) / S;
                    y = 0.25 * S;
                    z = (m[6] + m[9]) / S;
                    r = (m[8] - m[2]) / S;

                } else {
                    S = sqrt(1.0 + m[10] - m[0] - m[5]) * 2;
                    x = (m[8] + m[2]) / S;
                    y = (m[6] + m[9]) / S;
                    z = 0.25 * S;
                    r = (m[1] - m[4]) / S;
                }
            }

        }

        /**
         * Quaternion SLERP (spherical linear interpolation)
         * (for small interpolation uses simple linear interpolation)
         * Interpolate between two quaternions representing rotations.
         * \note This method is useful to achieve smooth interpolation between
         * any two orientations.The first quaternion represents the original orientation
         * of the axes, the second quaternion the final orientation of the axes.
         * @param qa first quaternion representing the initial orientation.
         * @param qb second quaternion representing the final orientation.
         * @param t a scalar between 0.0 and 1.0
         */
        void slerp(const CQuaternion<N>& qa, const CQuaternion<N>& qb, N t) {

            N theta, costheta, sintheta, ratio1, ratio2;

            //calcola il coseno fra gli angoli dei due quaternion
            costheta = qa.x * qb.x + qa.y * qb.y + qa.z * qb.z + qa.r * qb.r;

            //controlla che si stia percorrendo la strada piu' corta
            if ((1.0 + costheta) > EPSILON) {

                // se l'angolo non � troppo piccolo si usa lo slerp
                if ((1.0 - costheta) > EPSILON) {

                    theta = acos(costheta);
                    sintheta = sin(theta);
                    ratio1 = sin((1.0 - t) * theta) / sintheta;
                    ratio2 = sin(t * theta) / sintheta;

                } else {
                    // per piccoli angoli usa il LERP (interpolazione lineare)
                    ratio1 = 1.0 - t;
                    ratio2 = t;
                }

                x = ratio1 * qa.x + ratio2 * qb.x;
                y = ratio1 * qa.y + ratio2 * qb.y;
                z = ratio1 * qa.z + ratio2 * qb.z;
                r = ratio1 * qa.r + ratio2 * qb.r;

            } else {

                x = -qb.y;
                y = qb.x;
                z = -qb.r;
                r = qb.z;

                ratio1 = sin((1.0 - t) * static_cast<N> (pi_2));
                ratio2 = sin(t * static_cast<N> (pi_2));

                x = ratio1 * qa.x + ratio2 * x;
                y = ratio1 * qa.y + ratio2 * y;
                z = ratio1 * qa.z + ratio2 * z;
                r = ratio1 * qa.r + ratio2 * r;
            }
        }

    };
};



#endif // QUAT2_H_INCLUDED
