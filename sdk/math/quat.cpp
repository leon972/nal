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
 Code by Leonardo Berti (C) 2007
 Portable Math Engine - quaterions

 ****************************************************************/


#include "math/vectors.h"  //header principale (definizione dei vettori e operazioni base)
#include "math/quat.h"
#include "math/matrix.h"   //header per questo modulo

#ifdef _USENAMESPACE_
namespace mathengine {
    #endif

//crea un quaternion
    QUAT3D getQuat3D(const real r,const real x,const real y,const real z) {
        auto QUAT3D q;

        q.r=r;
        q.x=x;
        q.y=y;
        q.z=z;

        return q;
    }

    QUAT3D vectorToQuat3D(HVECTOR3D_PTR v) {
        auto QUAT3D q;

        q.r=0;
        q.x=v->x;
        q.y=v->y;
        q.z=v->z;

        return q;

    }

    void vectorToQuat3DEx(QUAT3D_PTR q,const HVECTOR3D& v) {
        q->r=0;
        q->x=v.x;
        q->y=v.y;
        q->z=v.z;
    }

    HVECTOR3D quatToVector(const QUAT3D& q) {
        auto HVECTOR3D r;

        r.x=q.x;
        r.y=q.y;
        r.z=q.z;
        r.w=1;

        return r;
    }

    void quatToVectorEx(HVECTOR3D_PTR v,const QUAT3D& q) {
        v->x=q.x;
        v->y=q.y;
        v->z=q.z;
        v->w=0;
    }

    //somma di quaternion
    QUAT3D sumQuat3D(const QUAT3D& q1,const QUAT3D& q2) {
        auto QUAT3D r;

        r.r=q1.r+q2.r;
        r.x=q1.x+q2.x;
        r.y=q1.y+q2.y;
        r.z=q1.z+q2.z;

        return r;
    }

    void sumQuat3DEx(QUAT3D_PTR result,const QUAT3D& q1,const QUAT3D& q2) {
        result->r=q1.r+q2.r;
        result->x=q1.x+q2.x;
        result->y=q1.y+q2.y;
        result->z=q1.z+q2.z;
    }

    //sottrazione
    //q1-q2
    QUAT3D subQuat3D(const QUAT3D& q1,const QUAT3D& q2) {
        auto QUAT3D r;

        r.r=q1.r-q2.r;
        r.x=q1.x-q2.x;
        r.y=q1.y-q2.y;
        r.z=q1.z-q2.z;

        return r;

    }

    void subQuat3DEx(QUAT3D_PTR result,const QUAT3D& q1,const QUAT3D& q2) {
        result->r=q1.r-q2.r;
        result->x=q1.x-q2.x;
        result->y=q1.y-q2.y;
        result->z=q1.z-q2.z;
    }

    //normalizzazione
    QUAT3D normalizeQuat3D(const QUAT3D& q) {
        QUAT3D r;

        real d=1/sqrtr(q.r*q.r+q.x*q.x+q.y*q.y+q.z*q.z);

        r.r=q.r*d;
        r.x=q.x*d;
        r.y=q.y*d;
        r.z=q.z*d;

        return r;
    }

    void normalizeQuat3DEx(QUAT3D_PTR q) {
        real d=1/sqrtr(q->r*q->r+q->x*q->x+q->y*q->y+q->z*q->z);

        q->r*=d;
        q->x*=d;
        q->y*=d;
        q->z*=d;
    }

    real magnitudeQuat3D(const QUAT3D& q) {
        return sqrtr(q.r*q.r+q.x*q.x+q.y*q.y+q.z*q.z);
    }

    real magnitudeQuat3DEx(QUAT3D_PTR q) {
        return sqrtr(q->r*q->r+q->x*q->x+q->y*q->y+q->z*q->z);
    }

    QUAT3D quatConjugate3D(const QUAT3D& q) {
        auto QUAT3D r;

        r.r=q.r;
        r.x=-q.x;
        r.y=-q.y;
        r.z=-q.z;

        return r;

    }

    void quatConjugateEx(QUAT3D_PTR result,const QUAT3D& q) {
        result->r=q.r;
        result->x=-q.x;
        result->y=-q.y;
        result->z=-q.z;
    }

    void printQuat3D(FILE* out,const QUAT3D& q) {
        fprintf(out,"(%f+%fi+%fj+%fk)",q.r,q.x,q.y,q.z);
    }

    //moltiplicazione fra quaternion q1*q2
    QUAT3D mulQuat3D(const QUAT3D& q1,const QUAT3D& q2) {
        auto QUAT3D r;

        //riduce il numero di moltiplicazione ragruppando
        real prd_0 = (q1.z - q1.y) * (q2.y - q2.z);
        real prd_1 = (q1.r + q1.x) * (q2.r + q2.x);
        real prd_2 = (q1.r - q1.x) * (q2.y + q2.z);
        real prd_3 = (q1.y + q1.z) * (q2.r - q2.x);
        real prd_4 = (q1.z - q1.x) * (q2.x - q2.y);
        real prd_5 = (q1.z + q1.x) * (q2.x + q2.y);
        real prd_6 = (q1.r + q1.y) * (q2.r - q2.z);
        real prd_7 = (q1.r - q1.y) * (q2.r + q2.z);
        real prd_8 = prd_5 + prd_6 + prd_7;
        real prd_9 = 0.5 * (prd_4 + prd_8);

        r.r = prd_0 + prd_9 - prd_5;
        r.x = prd_1 + prd_9 - prd_8;
        r.y = prd_2 + prd_9 - prd_7;
        r.z = prd_3 + prd_9 - prd_6;

        return r;
    }

    void mulQuat3DEx(QUAT3D_PTR result,const QUAT3D& q1,const QUAT3D& q2) {

        //riduce il numero di moltiplicazioni ragruppando
        real prd_0 = (q1.z - q1.y) * (q2.y - q2.z);
        real prd_1 = (q1.r + q1.x) * (q2.r + q2.x);
        real prd_2 = (q1.r - q1.x) * (q2.y + q2.z);
        real prd_3 = (q1.y + q1.z) * (q2.r - q2.x);
        real prd_4 = (q1.z - q1.x) * (q2.x - q2.y);
        real prd_5 = (q1.z + q1.x) * (q2.x + q2.y);
        real prd_6 = (q1.r + q1.y) * (q2.r - q2.z);
        real prd_7 = (q1.r - q1.y) * (q2.r + q2.z);
        real prd_8 = prd_5 + prd_6 + prd_7;
        real prd_9 = 0.5 * (prd_4 + prd_8);

        result->r=prd_0 + prd_9 - prd_5;
        result->x=prd_1 + prd_9 - prd_8;
        result->y=prd_2 + prd_9 - prd_7;
        result->z=prd_3 + prd_9 - prd_6;

    }

    //ottiene un quatenion a partire dal versore della direzione di rotazione e dall'angolo in gradi
    QUAT3D getRotQuat3D(const HVECTOR3D& versAxis,real theta) {
        auto QUAT3D q;

        real theta_div_2 = (0.5)*theta;

        real st2 = sinr(theta_div_2);

        q.x = st2 * versAxis.x;
        q.y = st2 * versAxis.y;
        q.z = st2 * versAxis.z;
        q.r = cosr(theta_div_2);

        return q;

    }

    //ottiene un quaternion che esegue la rotazione di theta radianti intorno al versore versAxis
    //nota: se versAxis � un versore , allora il quatenion ottenuto ha modulo 1
    void getRotQuat3DEx(QUAT3D_PTR q,const HVECTOR3D& versAxis,real theta) {
        real theta_div_2 = (0.5)*theta;

        real st2 = sinr(theta_div_2);

        q->x = st2 * versAxis.x;
        q->y = st2 * versAxis.y;
        q->z = st2 * versAxis.z;
        q->r = cosr(theta_div_2);
    }

    //q1*q2*q3
    //Il triploprodotto serve ad eseguire la rotazione del vettore rappresentato dal quaternion q2
    void quatTripleProduct(QUAT3D_PTR result,const QUAT3D& q1,const QUAT3D& q2,const QUAT3D& q3) {
        QUAT3D qtemp;

        mulQuat3DEx(&qtemp,q1,q2);
        mulQuat3DEx(result,qtemp,q3);
    }

    //ottiene la marice di rotazione a partire da un quaternion di rotazione
    /*

    Matrix =  [ 1 - 2y2 - 2z2    2xy - 2wz      2xz + 2wy
              2xy + 2wz    1 - 2x2 - 2z2    2yz - 2wx
              2xz - 2wy      2yz + 2wx    1 - 2x2 - 2y2 ]
    */
    void quatToRotMatrixEx(HMATRIX33_PTR m,const QUAT3D& q) {

        m->a30=m->a31=m->a32=0;
        m->a33=1; //coord-> omogenea

        real xy=q.x*q.y;
        real xz=q.x*q.z;
        real yz=q.y*q.z;
        real wx=q.x*q.r;
        real wy=q.y*q.r;
        real wz=q.z*q.r;

        real x2=q.x*q.x;
        real y2=q.y*q.y;
        real z2=q.z*q.z;

        m->a00=1-2*y2-2*z2;
        m->a01=2*xy-2*wz;
        m->a02=2*xz+2*wy;
        m->a03=0;

        m->a10=2*xy+2*wz;
        m->a11=1-2*x2-2*z2;
        m->a12=2*yz-2*wx;
        m->a13=0;

        m->a20=2*xz-2*wy;
        m->a21=2*yz+2*wx;
        m->a22=1-2*x2-2*y2;
        m->a23=0;

    }

    HMATRIX33 quatToRotMatrix(const QUAT3D& q) {
        auto HMATRIX33 m;

        quatToRotMatrixEx(&m,q);

        return m;

    }

    /**converte da matrice di rotazione a quternion
     NB:la matrice deve essere una matrice di rotazione!*/
    void rotMatToQuatEx(QUAT3D_PTR q,const HMATRIX33& m) {

        real T=1.0+m.a00+m.a11+m.a22;
        real S;

        if (T>EPSILON) { //esegue questo controllo per evitare distorsioni

            S=sqrt(T)*2;
            q->r=0.25*S;
            q->x=(m.a21-m.a12)/S;
            q->y=(m.a02-m.a20)/S;
            q->z=(m.a10-m.a01)/S;

        } else {
            //determina l'elemento diagonale maggiore
            if (m.a00>m.a11 && m.a00>m.a22) {

                S=sqrt(1.0+m.a00-m.a11-m.a22)*2;
                q->x = 0.25 * S;
                q->y = (m.a10 + m.a01) / S;
                q->z = (m.a02 + m.a20 ) / S;
                q->r = (m.a21 - m.a12) / S;

            } else if (m.a11>m.a22) {
                S  = sqrt( 1.0 + m.a11 - m.a00 - m.a22 ) * 2;
                q->x = (m.a10 + m.a01) / S;
                q->y = 0.25 * S;
                q->z = (m.a21 + m.a12) / S;
                q->r = (m.a02 - m.a20 ) / S;
            } else {
                S  = sqrt( 1.0 + m.a22 - m.a00 - m.a11 ) * 2;
                q->x = (m.a02 + m.a20 ) / S;
                q->y = (m.a21 + m.a12 ) / S;
                q->z = 0.25 * S;
                q->r = (m.a10 - m.a01 ) / S;
            }
        }
    }

    QUAT3D rotMatToQuat(const HMATRIX33& m) {
        auto QUAT3D q;

        rotMatToQuatEx(&q,m);

        return q;
    }

    /**SLERP:(spherical interpolation)
    restituisce un quaternion compreso fra il quaternion iniziale qa e quello finale qb
    basato sul parametro t compreso fra 0 e 1.t=0 rende qa t=1 rende qb
    */
    void quatSlerp(QUAT3D_PTR r,const QUAT3D& qa,const QUAT3D& qb,real t) {
        real theta,costheta,sintheta,ratio1,ratio2;

        //calcola il coseno fra gli angoli dei due quaternion
        costheta=qa.x*qb.x+qa.y*qb.y+qa.z*qb.z+qa.r*qb.r;

        //controlla che si stia percorrendo la strada piu' corta
        if ((1.0+costheta)>EPSILON) {

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

            r->x = ratio1 * qa.x + ratio2 * qb.x;
            r->y = ratio1 * qa.y + ratio2 * qb.y;
            r->z = ratio1 * qa.z + ratio2 * qb.z;
            r->r = ratio1 * qa.r + ratio2 * qb.r;

        } else {

            r->x = -qb.y;
            r->y = qb.x;
            r->z = -qb.r;
            r->r = qb.z;

            ratio1 = sin((1.0 - t) * (real)M_PI_2);
            ratio2 = sin(t * (real)M_PI_2);

            r->x = ratio1 * qa.x + ratio2 * r->x;
            r->y = ratio1 * qa.y + ratio2 * r->y;
            r->z = ratio1 * qa.z + ratio2 * r->z;
            r->r = ratio1 * qa.r + ratio2 * r->r;
        }
    }

    #ifdef _USENAMESPACE_
};
#endif
