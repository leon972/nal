/***************************************
 Copyright (c) 2011 by Leonardo Berti
 ***************************************/

/*!
  
 \file matrix2.h
 \brief Template-based geometric homogeneous 3D transformation matrix class.
 The elements of this matrix are stored in column major format, and can be used
 directly as argument of OpenGL functions.
 * 
 */

#ifndef MATRIX2_H_INCLUDED
#define MATRIX2_H_INCLUDED

#include "vectors3.h"
#include "vectors2.h"
#include <cstring>

namespace mathengine {

    using namespace std;
    
    /**
     * Geometric transformation matrix class in homogeneous coordinates.
     * This matrix is used to represent rotations and translations
     * og the axes in a 3 dimensions spaces.
     * The internal elements of the matrix have a column-major arrangement
     * in order to use directly with OpenGL.
     */

    template <class N> class CHMatrix33 {
    public:

        /**
         * Create a zero matrix
         */
        CHMatrix33() {
            reset();
        }

        /**
         * Copy constructor
         * @param m
         */
        CHMatrix33(const CHMatrix33<N>& m) {
            memcpy(a1, m.a1, 16 * sizeof (N));
            cur = a1;
        }

        /**
         * Copy assignement
         * @param m the matrix to assign
         * @return *this
         */
        CHMatrix33<N> & operator=(const CHMatrix33<N>& m) {
            if (&m != this) {
                memcpy(cur, m.cur, 16 * sizeof (N));
            }

            return *this;
        }

        /**
         * Get a matrix component by row and column.
         * \warning for optimization reason this method doesn't permfom bounds check.
         * @param row row (0..3)
         * @param col column (0..3)
         * @return the value of the component.
         */
        N get_value(int row, int col) const {
            return cur[col * 4 + row];
        }

        /**
         * Set a matrix component
         * \warning for optimization reason this method doesn't permfom bounds check.
         * @param row row of the component to set (0..4)
         * @param col column of the component to set (0..3)
         * @param value the value to set.
         */
        void set_value(int row, int col, const N value) {
            cur[col * 4 + row] = value;
        }

        /**
         * Get the internal representation of the matrix in column major form.
         * The destination array may be passed directly to OpenGL functions
         * @param dest the array with the matrix elements.this array must have 16 components at least.
         */
        void copy_elements(N* dest) const {
            memcpy(dest, cur, 16 * sizeof (N));
        }

        /**
         * Get the pointer to the array of the components of the matrix (column major form)
         * @return the pointer to the 16 components of the matrix.
         */
        const N* get_col_major_elements() {
            return cur;
        }

        /**
         * Copy the matrix components from a vector.
         * @param elements components of the matrix in column major form.
         */
        void set_elements(const N* elements) {
            memcpy(cur, elements, 16 * sizeof (N));
        }

        /**
         * Set the identity matrix
         */
        void set_identity() {

            cur[0] = 1;
            cur[5] = 1;
            cur[10] = 1;
            cur[15] = 1;
            cur[1] = cur[2] = cur[3] = cur[4] = 0;
            cur[6] = cur[7] = cur[8] = cur[9] = 0;
            cur[11] = cur[12] = cur[13] = cur[14] = 0;
        }

        /**
         * Multiply this matrix by m
         * @param m
         */
        void mul(const CHMatrix33<N>& m) {

            //1� riga
            next[0] = cur[0] * m.cur[0] + cur[4] * m.cur[1] + cur[8] * m.cur[2] + cur[12] * m.cur[3];
            next[4] = cur[0] * m.cur[4] + cur[4] * m.cur[5] + cur[8] * m.cur[6] + cur[12] * m.cur[7];
            next[8] = cur[0] * m.cur[8] + cur[4] * m.cur[9] + cur[8] * m.cur[10] + cur[12] * m.cur[11];
            next[12] = cur[0] * m.cur[12] + cur[4] * m.cur[13] + cur[8] * m.cur[14] + cur[12] * m.cur[15];

            //2� riga
            next[1] = cur[1] * m.cur[0] + cur[5] * m.cur[1] + cur[9] * m.cur[2] + cur[13] * m.cur[3];
            next[5] = cur[1] * m.cur[4] + cur[5] * m.cur[5] + cur[9] * m.cur[6] + cur[13] * m.cur[7];
            next[9] = cur[1] * m.cur[8] + cur[5] * m.cur[9] + cur[9] * m.cur[10] + cur[13] * m.cur[11];
            next[13] = cur[1] * m.cur[12] + cur[5] * m.cur[13] + cur[9] * m.cur[14] + cur[13] * m.cur[15];

            //3� riga
            next[2] = cur[2] * m.cur[0] + cur[6] * m.cur[1] + cur[10] * m.cur[2] + cur[14] * m.cur[3];
            next[6] = cur[2] * m.cur[4] + cur[6] * m.cur[5] + cur[10] * m.cur[6] + cur[14] * m.cur[7];
            next[10] = cur[2] * m.cur[8] + cur[6] * m.cur[9] + cur[10] * m.cur[10] + cur[14] * m.cur[11];
            next[14] = cur[2] * m.cur[12] + cur[6] * m.cur[13] + cur[10] * m.cur[14] + cur[14] * m.cur[15];

            //4� riga
            next[3] = cur[3] * m.cur[0] + cur[7] * m.cur[1] + cur[11] * m.cur[2] + cur[15] * m.cur[3];
            next[7] = cur[3] * m.cur[4] + cur[7] * m.cur[5] + cur[11] * m.cur[6] + cur[15] * m.cur[7];
            next[11] = cur[3] * m.cur[8] + cur[7] * m.cur[9] + cur[11] * m.cur[10] + cur[15] * m.cur[11];
            next[15] = cur[3] * m.cur[12] + cur[7] * m.cur[13] + cur[11] * m.cur[14] + cur[15] * m.cur[15];

            N* temp = next;
            next = cur;
            cur = temp;
        }

        /**
         * Multiply the rotation part of this matrix by the rotation part of m
         * @param m a transformation matrix
         */

        void mul3x3(const CHMatrix33<N>& m) {
            //1� riga
            next[0] = cur[0] * m.cur[0] + cur[4] * m.cur[1] + cur[8] * m.cur[2];
            next[4] = cur[0] * m.cur[4] + cur[4] * m.cur[5] + cur[8] * m.cur[6];
            next[8] = cur[0] * m.cur[8] + cur[4] * m.cur[9] + cur[8] * m.cur[10];
            next[12] = cur[12];

            //2� riga
            next[1] = cur[1] * m.cur[0] + cur[5] * m.cur[1] + cur[9] * m.cur[2];
            next[5] = cur[1] * m.cur[4] + cur[5] * m.cur[5] + cur[9] * m.cur[6];
            next[9] = cur[1] * m.cur[8] + cur[5] * m.cur[9] + cur[9] * m.cur[10];
            next[13] = cur[13];

            //3� riga
            next[2] = cur[2] * m.cur[0] + cur[6] * m.cur[1] + cur[10] * m.cur[2];
            next[6] = cur[2] * m.cur[4] + cur[6] * m.cur[5] + cur[10] * m.cur[6];
            next[10] = cur[2] * m.cur[8] + cur[6] * m.cur[9] + cur[10] * m.cur[10];
            next[14] = cur[14];

            //4� riga
            next[3] = cur[3] * m.cur[0] + cur[7] * m.cur[1] + cur[11] * m.cur[2];
            next[7] = cur[3] * m.cur[4] + cur[7] * m.cur[5] + cur[11] * m.cur[6];
            next[11] = cur[3] * m.cur[8] + cur[7] * m.cur[9] + cur[11] * m.cur[10];
            next[15] = cur[15];

            N* temp = next;
            next = cur;
            cur = temp;

        }

        /**
         * Add a transformation matrix to this matrix
         * @param m matrix to add
         */
        void add(const CHMatrix33<N>& m) {
            cur[0] += m.cur[0];
            cur[1] += m.cur[1];
            cur[2] += m.cur[2];
            cur[3] += m.cur[3];
            cur[4] += m.cur[4];
            cur[5] += m.cur[5];
            cur[6] += m.cur[6];
            cur[7] += m.cur[7];
            cur[8] += m.cur[8];
            cur[9] += m.cur[9];
            cur[10] += m.cur[10];
            cur[11] += m.cur[11];
            cur[12] += m.cur[12];
            cur[13] += m.cur[13];
            cur[14] += m.cur[14];
            cur[15] += m.cur[15];
        }

        /**
         * Subtract a transformation matrix from this matrix.
         * @param m matrix to subtract.
         */
        void sub(const CHMatrix33<N>& m) {
            cur[0] -= m.cur[0];
            cur[1] -= m.cur[1];
            cur[2] -= m.cur[2];
            cur[3] -= m.cur[3];
            cur[4] -= m.cur[4];
            cur[5] -= m.cur[5];
            cur[6] -= m.cur[6];
            cur[7] -= m.cur[7];
            cur[8] -= m.cur[8];
            cur[9] -= m.cur[9];
            cur[10] -= m.cur[10];
            cur[11] -= m.cur[11];
            cur[12] -= m.cur[12];
            cur[13] -= m.cur[13];
            cur[14] -= m.cur[14];
            cur[15] -= m.cur[15];
        }

        /**
         * Multiply this matrix by the vector v and return the result in result.
         * @param result this matrix multiplied by v
         * @param v a 3 omponents vector in homogeneous coordinates.
         */
        void mul_vect(CHVector3<N>* result, const CHVector3<N>& v) const {
            result->x = cur[0] * v.x + cur[4] * v.y + cur[8] * v.z + cur[12];
            result->y = cur[1] * v.x + cur[5] * v.y + cur[9] * v.z + cur[13];
            result->z = cur[2] * v.x + cur[6] * v.y + cur[10] * v.z + cur[14];
            result->w = 1; //m->a30*v->x+m->a31*v->y+m->a32*v->z+m->a33*v->w;
        }

        /**
         * Multiply this matrix by the vector v.
         * @param v
         * @return the result of the multiplication of this matrix by v.
         */
        CHVector3<N> mul_vect(const CHVector3<N>& v) const {
            CHVector3<N> r;
            mul_vect(&r, v);
            return r;
        }

        /**
         * Multiply this matrix by a 2 components vector (assumes v.z=0)
         * @param result this matrix multiplied by b
         * @param v a 2 components vector.
         */
        void mul_vect(CVector2<N>* result, const CVector2<N>&v) const {

            result->x = cur[0] * v.x + cur[4] * v.y + cur[12];
            result->y = cur[1] * v.x + cur[5] * v.y + cur[13];
        }

        /**
         * Multiply this matrix by a 2 components vector (assumes v.z=0)
         * @return this matrix multiplied by b
         * @param v a 2 components vector.
         */
        CVector2<N> mul_vect(const CVector2<N>& v) const {

            CVector2<N> r;
            mul_vect(&r, v);
            return r;
        }

        /**
         * Transpose this matrix.
         * aij=aji
         */
        void transpose() {
            next[0] = cur[0];
            next[5] = cur[5];
            next[10] = cur[10];
            next[15] = cur[15];

            next[4] = cur[1];
            next[1] = cur[4];
            next[8] = cur[2];
            next[2] = cur[8];
            next[12] = cur[3];
            next[3] = cur[12];
            next[9] = cur[6];
            next[6] = cur[9];
            next[13] = cur[7];
            next[7] = cur[13];
            next[14] = cur[11];
            next[11] = cur[14];

            N* temp = next;
            next = cur;
            cur = temp;
        }

        /**
         * Multiply this matrix by a scalar.
         * @param scalar
         */
        void mul_scalar(N scalar) {
            cur[0] *= scalar;
            cur[1] *= scalar;
            cur[2] *= scalar;
            cur[3] *= scalar;
            cur[4] *= scalar;
            cur[5] *= scalar;
            cur[6] *= scalar;
            cur[7] *= scalar;
            cur[8] *= scalar;
            cur[9] *= scalar;
            cur[10] *= scalar;
            cur[11] *= scalar;
            cur[12] *= scalar;
            cur[13] *= scalar;
            cur[14] *= scalar;
            cur[15] *= scalar;
        }

        /**
         * Get the determinant of the rotation part of this matrix.
         * @return the determinant.
         */
        N det3x3() const {
            return cur[0]*(cur[5] * cur[10] - cur[6] * cur[9]) - cur[4]*(cur[1] * cur[10] - cur[2] * cur[9]) + cur[8]*(cur[1] * cur[6] - cur[2] * cur[5]);
        }

        /**
         * Inverts this matrix.
         * @return true if the matri can be inverted.
         */
        bool inverse3x3() {

            N d = this->det3x3();

            if (d == 0) return false; //non invertibile

            d = 1.0 / d;

            next[0] = d * (cur[5] * cur[10] - cur[6] * cur[9]);
            next[4] = d * (-cur[4] * cur[10] + cur[6] * cur[8]);
            next[8] = d * (cur[4] * cur[9] - cur[5] * cur[8]);

            next[1] = d * (-cur[1] * cur[10] + cur[2] * cur[9]);
            next[5] = d * (cur[0] * cur[10] - cur[2] * cur[8]);
            next[9] = d * (-cur[0] * cur[9] + cur[1] * cur[8]);

            next[2] = d * (cur[1] * cur[6] - cur[2] * cur[5]);
            next[6] = d * (-cur[0] * cur[6] + cur[2] * cur[4]);
            next[10] = d * (cur[0] * cur[5] - cur[4] * cur[1]);

            //coord. omogenee
            next[3] = next[7] = next[11] = next[14] = next[13] = next[12] = 0;
            next[15] = 1;

            N* temp = next;
            next = cur;
            cur = temp;

            return true;
        }

        /**
         * Creates the matrix that performs the rotation of theta radians about an axis.
         * @param vers versor of the axis of rotation (must have magnitude=1)
         * @param theta counter-clockwise rotation in radians.
         */
        void make_rot_matrix(const CHVector3<N>& vers, const N theta) {

            N st = sin(theta);
            N ct = cos(theta);
            N ct1 = 1 - ct;
            N xy = vers.x * vers.y;
            N xz = vers.x * vers.z;
            N yz = vers.y * vers.z;
            N x2 = vers.x * vers.x;
            N y2 = vers.y * vers.y;
            N z2 = vers.z * vers.z;

            cur[3] = cur[7] = cur[11] = cur[14] = cur[13] = cur[12] = 0;
            cur[15] = 1;

            cur[0] = ct + ct1*x2;
            cur[4] = ct1 * xy - st * vers.z;
            cur[8] = ct1 * xz + st * vers.y;

            cur[1] = ct1 * xy + st * vers.z;
            cur[5] = ct + ct1*y2;
            cur[9] = ct1 * yz - st * vers.x;

            cur[2] = ct1 * xz - st * vers.y;
            cur[6] = ct1 * yz + st * vers.x;
            cur[10] = ct + ct1*z2;

        }

        /**
         * Creates the counter-clockwise rotation matrix about the Z axis.
         * @param theta angle of rotation in radians.
         */
        void make_rot_z_matrix(N theta) {

            N ct = cos(theta);
            N st = sin(theta);

            cur[0] = ct;
            cur[4] = st;
            cur[8] = 0;
            cur[12] = 0;
            cur[1] = -st;
            cur[5] = ct;
            cur[9] = 0;
            cur[13] = 0;
            cur[2] = cur[6] = 0;
            cur[10] = 1;
            cur[14] = 0;

            cur[3] = cur[7] = cur[11] = 0;
            cur[15] = 1;
        }

        /**
         * Creates the counter-clockwise rotation matrix about the Y axis.
         * @param theta angle of rotation in radians.
         */
        void make_rot_y_matrix(N theta) {

            N ct = cos(theta);
            N st = sin(theta);
            cur[0] = ct;
            cur[4] = 0;
            cur[8] = -st;
            cur[12] = 0;
            cur[1] = 0;
            cur[5] = 1;
            cur[9] = 0;
            cur[13] = 0;
            cur[2] = st;
            cur[6] = 0;
            cur[10] = ct;
            cur[14] = 0;
            cur[3] = cur[7] = cur[11] = 0;
            cur[15] = 1;
        }

        /**
         * Creates the counter-clockwise rotation matrix about the X axis.
         * @param theta angle of rotation in radians.
         */
        void make_rot_x_matrix(N theta) {

            N ct = cos(theta);
            N st = sin(theta);

            cur[0] = 1;
            cur[4] = 0;
            cur[8] = 0;
            cur[12] = 0;
            cur[1] = 0;
            cur[5] = ct;
            cur[9] = st;
            cur[13] = 0;
            cur[2] = 0;
            cur[6] = -st;
            cur[10] = ct;
            cur[14] = 0;
            cur[3] = cur[7] = cur[11] = 0;
            cur[15] = 1;
        }

        /**
         * Set the rotational part of this matrix to the X axis rotation matrix.
         * Doesn't affect the translation components.
         * @param theta
         */
        void set_rot_x(N theta) {
            N ct = cos(theta);
            N st = sin(theta);

            cur[0] = 1;
            cur[4] = 0;
            cur[8] = 0;

            cur[1] = 0;
            cur[5] = ct;
            cur[9] = st;

            cur[2] = 0;
            cur[6] = -st;
            cur[10] = ct;
        }

        /**
         * Set the rotational part of this matrix to the Z axis rotation matrix.
         * Doesn't affect the translation components.
         * @param theta
         */
        void set_rot_z(N theta) {
            N ct = cos(theta);
            N st = sin(theta);

            cur[0] = ct;
            cur[4] = st;
            cur[8] = 0;

            cur[1] = -st;
            cur[5] = ct;
            cur[9] = 0;

            cur[2] = cur[6] = 0;
            cur[10] = 1;

        }

        /**
         * Set the rotational part of this matrix to the Y axis rotation matrix.
         * Doesn't affect the translation components.
         * @param theta
         */
        void set_rot_y(N theta) {
            N ct = cos(theta);
            N st = sin(theta);
            cur[0] = ct;
            cur[4] = 0;
            cur[8] = -st;

            cur[1] = 0;
            cur[5] = 1;
            cur[9] = 0;

            cur[2] = st;
            cur[6] = 0;
            cur[10] = ct;

        }

        /**
         * Multiply this matrix by the X axis rotation matrix.
         * This corresponds to a rotation about the local axes.
         * @param theta counter-clockwise rotation angle in radians.
         */
        void rotate_x(N theta) {
            N ct, st;

            ct = cos(theta);
            st = sin(theta);

            next[4] = cur[4];
            next[5] = cur[5];
            next[6] = cur[6];
            next[8] = cur[8];
            next[9] = cur[9];
            next[10] = cur[10];

            cur[4] = next[4] * ct - next[8] * st;
            cur[5] = next[5] * ct - next[9] * st;
            cur[6] = next[6] * ct - next[10] * st;
            cur[8] = next[4] * st + next[8] * ct;
            cur[9] = next[5] * st + next[9] * ct;
            cur[10] = next[6] * st + next[10] * ct;
        }

        /**
         * Multiply this matrix by the Y axis rotation matrix.
         * This corresponds to a rotation about the local axes.
         * @param theta counter-clockwise rotation angle in radians.
         */
        void rotate_y(N theta) {
            N ct, st;

            ct = cos(theta);
            st = sin(theta);

            next[0] = cur[0];
            next[1] = cur[1];
            next[2] = cur[2];
            next[8] = cur[8];
            next[9] = cur[9];
            next[10] = cur[10];

            cur[0] = next[0] * ct + next[8] * st;
            cur[1] = next[1] * ct + next[9] * st;
            cur[2] = next[2] * ct + next[10] * st;
            cur[8] = next[8] * ct - next[0] * st;
            cur[9] = next[9] * ct - next[1] * st;
            cur[10] = next[10] * ct - next[2] * st;

        }

        /**
         * Multiply this matrix by the Z axis rotation matrix.
         * This corresponds to a rotation about the local axes.
         * @param theta counter-clockwise rotation angle in radians.
         */
        void rotate_z(N theta) {
            N ct, st;

            ct = cos(theta);
            st = sin(theta);

            next[0] = cur[0];
            next[1] = cur[1];
            next[2] = cur[2];
            next[4] = cur[4];
            next[5] = cur[5];
            next[6] = cur[6];

            cur[0] = next[0] * ct - next[4] * st;
            cur[1] = next[1] * ct - next[5] * st;
            cur[2] = next[2] * ct - next[6] * st;

            cur[4] = next[0] * st + next[4] * ct;
            cur[5] = next[1] * st + next[5] * ct;
            cur[6] = next[2] * st + next[6] * ct;
        }

        /**
         * Create the Euler transformation matrix.
         * @param psi
         * @param theta
         * @param phi
         */
        void make_euler_matrix(const N psi, const N theta, const N phi) {

            N cp = cos(psi);
            N sp = sin(psi);

            N ct = cos(theta);
            N st = sin(theta);

            N cf = cos(phi);
            N sf = sin(phi);

            cur[0] = cf * cp - sf * ct*sp;
            cur[4] = cf * sp + sf * ct*cp;
            cur[8] = sf*st;
            cur[12] = 0;

            cur[1] = -sf * cp - cf * ct*sp;
            cur[5] = -sf * sp + cf * ct*cp;
            cur[9] = cf*st;
            cur[13] = 0;

            cur[2] = st*sp;
            cur[6] = -st*cp;
            cur[10] = ct;
            cur[14] = 0;

            cur[3] = 0;
            cur[7] = 0;
            cur[11] = 0;
            cur[15] = 1;
        }

        /**
         * Create the translation marix.
         *
         * @param tx x translation
         * @param ty y translation
         * @param tz z translation
         */
        void make_translation_matrix(const N tx, const N ty, const N tz) {

            cur[0] = 1;
            cur[4] = 0;
            cur[8] = 0;
            cur[12] = tx;

            cur[1] = 0;
            cur[5] = 1;
            cur[9] = 0;
            cur[13] = ty;

            cur[2] = 0;
            cur[6] = 0;
            cur[10] = 1;
            cur[14] = tz;

            cur[3] = cur[7] = cur[11] = 0;
            cur[15] = 1;
        }

        /**
         * Translation about the local axes
         * @param tx
         * @param ty
         * @param tz
         */
        void translate(const N tx, const N ty, const N tz) {
            next[12] = cur[12];
            next[13] = cur[13];
            next[14] = cur[14];

            cur[12] = cur[0] * tx + cur[4] * ty + cur[8] * tz + next[12];
            cur[13] = cur[1] * tx + cur[5] * ty + cur[9] * tz + next[13];
            cur[14] = cur[2] * tx + cur[6] * ty + cur[10] * tz + next[14];
        }

        /**
         * Set the translation components.
         * @param tx
         * @param ty
         * @param tz
         */
        void set_position(const N tx, const N ty, const N tz) {
            cur[12] = tx;
            cur[13] = ty;
            cur[14] = tz;
        }

        /**
         * Create the scaling matrix.
         *
         * @param sx x scale factor
         * @param sy y scale factor
         * @param sz z scale factor
         */
        void make_scale_matrix(const N sx, const N sy, const N sz) {

            cur[0] = sx;
            cur[4] = 0;
            cur[8] = 0;
            cur[12] = 0;

            cur[1] = 0;
            cur[5] = sy;
            cur[9] = 0;
            cur[13] = 0;

            cur[2] = 0;
            cur[6] = 0;
            cur[10] = sz;
            cur[14] = 0;

            cur[3] = 0;
            cur[7] = 0;
            cur[11] = 0;
            cur[15] = 1;
        }

    private:

        //NB: column major!
        N a1[16];
        N a2[16];

        /**
         *Internal cache
         */
        N* cur;
        N* next;

        void reset() {
            memset(a1, 0, sizeof (N)*16);
            cur = a1;
            next = a2;
        }

    }; //end class

};

#endif // MATRIX2_H_INCLUDED
