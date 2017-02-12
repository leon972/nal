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
 Portable Math Engine - matrici 2x2 3x3 e 4x4 (coordinate omogenee)
 definizioni e operazioni principali

 ****************************************************************/


#include "math/vectors.h"  //header principale (definizione dei vettori e operazioni base)
#include "math/matrix.h"   //header per questo modulo

#ifdef _USENAMESPACE_
namespace mathengine
{
#endif

void getIdentityMatrix22Ex(MATRIX22_PTR m)
{
    m->a00=m->a11=1;
    m->a01=m->a10=0;
}

void getIdentityMatrix33Ex(MATRIX33_PTR m)
{
    m->a00=m->a11=m->a22=1;
    m->a01=m->a02=0;
    m->a10=m->a12=0;
    m->a20=m->a21=0;
}

void getIdentityHMatrix33Ex(HMATRIX33_PTR m)
{
    m->a00=m->a11=m->a22=m->a33=1;
    m->a01=m->a02=m->a03=0;
    m->a10=m->a12=m->a13=0;
    m->a20=m->a21=m->a23=0;
    m->a30=m->a31=m->a32=0;
}

MATRIX22 getIdentityMatrix22(void)
{
    auto MATRIX22 m;

    m.a00=m.a11=1;
    m.a01=m.a10=0;

    return m;
}

MATRIX33 getIdentityMatrix33(void)
{
    auto MATRIX33 m;

    m.a00=m.a11=m.a22=1;
    m.a01=m.a02=0;
    m.a10=m.a12=0;
    m.a20=m.a21=0;

    return m;
}

HMATRIX33 getIdentityHMatrix33(void)
{
    auto HMATRIX33 m;

    m.a00=m.a11=m.a22=m.a33=1;
    m.a01=m.a02=m.a03=0;
    m.a10=m.a12=m.a13=0;
    m.a20=m.a21=m.a23=0;
    m.a30=m.a31=m.a32=0;

    return m;
}

//moltiplicazioni fra matrici (m1 x m2)

void mulMatrix33(MATRIX33_PTR result,const MATRIX33& m1,const MATRIX33& m2)
{
    //1� riga
    result->a00=m1.a00*m2.a00+m1.a01*m2.a10+m1.a02*m2.a20;
    result->a01=m1.a00*m2.a01+m1.a01*m2.a11+m1.a02*m2.a21;
    result->a02=m1.a00*m2.a02+m1.a01*m2.a12+m1.a02*m2.a22;
    //2� riga
    result->a10=m1.a10*m2.a00+m1.a11*m2.a10+m1.a12*m2.a20;
    result->a11=m1.a10*m2.a01+m1.a11*m2.a11+m1.a12*m2.a21;
    result->a12=m1.a10*m2.a02+m1.a11*m2.a12+m1.a12*m2.a22;
    //3� riga
    result->a20=m1.a20*m2.a00+m1.a21*m2.a10+m1.a22*m2.a20;
    result->a21=m1.a20*m2.a01+m1.a21*m2.a11+m1.a22*m2.a21;
    result->a22=m1.a20*m2.a02+m1.a21*m2.a12+m1.a22*m2.a22;

}

void mulHMatrix33(HMATRIX33_PTR result,const HMATRIX33& m1,const HMATRIX33& m2)
{
    //1� riga
    result->a00=m1.a00*m2.a00+m1.a01*m2.a10+m1.a02*m2.a20+m1.a03*m2.a30;
    result->a01=m1.a00*m2.a01+m1.a01*m2.a11+m1.a02*m2.a21+m1.a03*m2.a31;
    result->a02=m1.a00*m2.a02+m1.a01*m2.a12+m1.a02*m2.a22+m1.a03*m2.a32;
    result->a03=m1.a00*m2.a03+m1.a01*m2.a13+m1.a02*m2.a23+m1.a03*m2.a33;

    //2� riga
    result->a10=m1.a10*m2.a00+m1.a11*m2.a10+m1.a12*m2.a20+m1.a13*m2.a30;
    result->a11=m1.a10*m2.a01+m1.a11*m2.a11+m1.a12*m2.a21+m1.a13*m2.a31;
    result->a12=m1.a10*m2.a02+m1.a11*m2.a12+m1.a12*m2.a22+m1.a13*m2.a32;
    result->a13=m1.a10*m2.a03+m1.a11*m2.a13+m1.a12*m2.a23+m1.a13*m2.a33;

    //3� riga
    result->a20=m1.a20*m2.a00+m1.a21*m2.a10+m1.a22*m2.a20+m1.a23*m2.a30;
    result->a21=m1.a20*m2.a01+m1.a21*m2.a11+m1.a22*m2.a21+m1.a23*m2.a31;
    result->a22=m1.a20*m2.a02+m1.a21*m2.a12+m1.a22*m2.a22+m1.a23*m2.a32;
    result->a23=m1.a20*m2.a03+m1.a21*m2.a13+m1.a22*m2.a23+m1.a23*m2.a33;

    //4� riga
    result->a30=m1.a30*m2.a00+m1.a31*m2.a10+m1.a32*m2.a20+m1.a33*m2.a30;
    result->a31=m1.a30*m2.a01+m1.a31*m2.a11+m1.a32*m2.a21+m1.a33*m2.a31;
    result->a32=m1.a30*m2.a02+m1.a31*m2.a12+m1.a32*m2.a22+m1.a33*m2.a32;
    result->a33=m1.a30*m2.a03+m1.a31*m2.a13+m1.a32*m2.a23+m1.a33*m2.a33;

}

//prodotto fra matrice 2x2 e vettore
void mulMatrixVect2(VECTOR2D_PTR result,const MATRIX22& m,const VECTOR2D& v)
{
    result->x=m.a00*v.x+m.a01*v.y;
    result->y=m.a10*v.x+m.a11*v.y;
}

//prodotto fra matrice 3X3 e vettore
void mulMatrixVect3(VECTOR3D_PTR result,const MATRIX33& m,const VECTOR3D& v)
{
    result->x=m.a00*v.x+m.a01*v.y+m.a02*v.z;
    result->y=m.a10*v.x+m.a11*v.y+m.a12*v.z;
    result->z=m.a20*v.x+m.a21*v.y+m.a22*v.z;
}

void mulMatrixHVect3Ex(HVECTOR3D_PTR result,const HMATRIX33& m,const HVECTOR3D& v)
{
    result->x=m.a00*v.x+m.a01*v.y+m.a02*v.z+m.a03*v.w;
    result->y=m.a10*v.x+m.a11*v.y+m.a12*v.z+m.a13*v.w;
    result->z=m.a20*v.x+m.a21*v.y+m.a22*v.z+m.a23*v.w;
    result->w=1; //m->a30*v->x+m->a31*v->y+m->a32*v->z+m->a33*v->w;
}

HVECTOR3D mulMatrixHVect3(const HMATRIX33& m,const HVECTOR3D& v)
{
    auto HVECTOR3D result;

    result.x=m.a00*v.x+m.a01*v.y+m.a02*v.z+m.a03*v.w;
    result.y=m.a10*v.x+m.a11*v.y+m.a12*v.z+m.a13*v.w;
    result.z=m.a20*v.x+m.a21*v.y+m.a22*v.z+m.a23*v.w;
    result.w=1;

    return result;
}

//calcola la matrice trasposta
void transposeMatrixEx(HMATRIX33_PTR result,const HMATRIX33& m)
{
    result->a00=m.a00;
    result->a01=m.a10;
    result->a02=m.a20;
    result->a03=m.a30;


    result->a10=m.a01;
    result->a11=m.a11;
    result->a12=m.a21;
    result->a13=m.a31;

    result->a20=m.a02;
    result->a21=m.a12;
    result->a22=m.a22;
    result->a23=m.a32;

    result->a30=m.a03;
    result->a31=m.a13;
    result->a32=m.a23;
    result->a33=m.a33;
}

HMATRIX33 transposeMatrix(const HMATRIX33& m)
{
    auto HMATRIX33 result;

    result.a00=m.a00;
    result.a01=m.a10;
    result.a02=m.a20;
    result.a03=m.a30;


    result.a10=m.a01;
    result.a11=m.a11;
    result.a12=m.a21;
    result.a13=m.a31;

    result.a20=m.a02;
    result.a21=m.a12;
    result.a22=m.a22;
    result.a23=m.a32;

    result.a30=m.a03;
    result.a31=m.a13;
    result.a32=m.a23;
    result.a33=m.a33;

    return result;
}

//somma fra matrici
void addMatrix33(MATRIX33_PTR result,const MATRIX33& m1,const MATRIX33&m2)
{
    result->a00=m1.a00+m2.a00;
    result->a01=m1.a01+m2.a01;
    result->a02=m1.a02+m2.a02;

    result->a10=m1.a10+m2.a10;
    result->a11=m1.a11+m2.a11;
    result->a12=m1.a12+m2.a12;

    result->a20=m1.a20+m2.a20;
    result->a21=m1.a21+m2.a21;
    result->a22=m1.a22+m2.a22;
}

//aggiunge alla matrice dest la matrice m
void addMatrix33Ex(MATRIX33_PTR dest,const MATRIX33& m)
{
    dest->a00+=m.a00;
    dest->a01+=m.a01;
    dest->a02+=m.a02;
    dest->a10+=m.a10;
    dest->a11+=m.a11;
    dest->a12+=m.a12;
    dest->a20+=m.a20;
    dest->a21+=m.a21;
    dest->a22+=m.a22;
}

//somma fra matrici in coordinate omogenee
void addHMatrix33(HMATRIX33_PTR result,const HMATRIX33& m1,const HMATRIX33&m2)
{
    result->a00=m1.a00+m2.a00;
    result->a01=m1.a01+m2.a01;
    result->a02=m1.a02+m2.a02;
    result->a03=m1.a03+m2.a03;

    result->a10=m1.a10+m2.a10;
    result->a11=m1.a11+m2.a11;
    result->a12=m1.a12+m2.a12;
    result->a13=m1.a13+m2.a13;

    result->a20=m1.a20+m2.a20;
    result->a21=m1.a21+m2.a21;
    result->a22=m1.a22+m2.a22;
    result->a23=m1.a23+m2.a23;

    result->a30=m1.a30+m2.a30;
    result->a31=m1.a31+m2.a31;
    result->a32=m1.a32+m2.a32;
    result->a33=m1.a33+m2.a33;

}

//aggiunge alla matrice dest la matrice m
void addHMatrix33Ex(HMATRIX33_PTR dest,const HMATRIX33& m)
{
    dest->a00+=m.a00;
    dest->a01+=m.a01;
    dest->a02+=m.a02;
    dest->a03+=m.a03;

    dest->a10+=m.a10;
    dest->a11+=m.a11;
    dest->a12+=m.a12;
    dest->a13+=m.a13;

    dest->a20+=m.a20;
    dest->a21+=m.a21;
    dest->a22+=m.a22;
    dest->a23+=m.a23;

    dest->a30+=m.a30;
    dest->a31+=m.a31;
    dest->a32+=m.a32;
    dest->a33+=m.a33;
}

void subMatrix33(MATRIX33_PTR result,const MATRIX33& m1,const MATRIX33&m2)
{
    result->a00=m1.a00-m2.a00;
    result->a01=m1.a01-m2.a01;
    result->a02=m1.a02-m2.a02;

    result->a10=m1.a10-m2.a10;
    result->a11=m1.a11-m2.a11;
    result->a12=m1.a12-m2.a12;

    result->a20=m1.a20-m2.a20;
    result->a21=m1.a21-m2.a21;
    result->a22=m1.a22-m2.a22;

}
//aggiunge alla matrice dest la matrice m
void subMatrix33Ex(MATRIX33_PTR dest,const MATRIX33& m)
{
    dest->a00-=m.a00;
    dest->a01-=m.a01;
    dest->a02-=m.a02;
    dest->a10-=m.a10;
    dest->a11-=m.a11;
    dest->a12-=m.a12;
    dest->a20-=m.a20;
    dest->a21-=m.a21;
    dest->a22-=m.a22;

}
//somma fra matrici in coordinate omogenee
void subHMatrix33(HMATRIX33_PTR result,const HMATRIX33& m1,const HMATRIX33&m2)
{
    result->a00=m1.a00-m2.a00;
    result->a01=m1.a01-m2.a01;
    result->a02=m1.a02-m2.a02;
    result->a03=m1.a03-m2.a03;

    result->a10=m1.a10-m2.a10;
    result->a11=m1.a11-m2.a11;
    result->a12=m1.a12-m2.a12;
    result->a13=m1.a13-m2.a13;

    result->a20=m1.a20-m2.a20;
    result->a21=m1.a21-m2.a21;
    result->a22=m1.a22-m2.a22;
    result->a23=m1.a23-m2.a23;

    result->a30=m1.a30-m2.a30;
    result->a31=m1.a31-m2.a31;
    result->a32=m1.a32-m2.a32;
    result->a33=m1.a33-m2.a33;
}

//aggiunge alla matrice dest la matrice m
void subHMatrix33Ex(HMATRIX33_PTR dest,const HMATRIX33& m)
{
    dest->a00-=m.a00;
    dest->a01-=m.a01;
    dest->a02-=m.a02;
    dest->a03-=m.a03;

    dest->a10-=m.a10;
    dest->a11-=m.a11;
    dest->a12-=m.a12;
    dest->a13-=m.a13;

    dest->a20-=m.a20;
    dest->a21-=m.a21;
    dest->a22-=m.a22;
    dest->a23-=m.a23;

    dest->a30-=m.a30;
    dest->a31-=m.a31;
    dest->a32-=m.a32;
    dest->a33-=m.a33;

}

void mulScalMatrix33(MATRIX33_PTR result,real a,const MATRIX33& m)
{
    result->a00=a*m.a00;
    result->a01=a*m.a01;
    result->a02=a*m.a02;

    result->a10=a*m.a10;
    result->a11=a*m.a11;
    result->a12=a*m.a12;

    result->a20=a*m.a20;
    result->a21=a*m.a21;
    result->a22=a*m.a22;
}

void mulScalMatrix33Ex(MATRIX33_PTR dest,real a)
{
    dest->a00*=a;
    dest->a01*=a;
    dest->a02*=a;

    dest->a10*=a;
    dest->a11*=a;
    dest->a12*=a;

    dest->a20*=a;
    dest->a21*=a;
    dest->a22*=a;
}

void mulScalHMatrix33(HMATRIX33_PTR result,real a,const HMATRIX33& m)
{
    result->a00=a*m.a00;
    result->a01=a*m.a01;
    result->a02=a*m.a02;
    result->a03=a*m.a03;

    result->a10=a*m.a10;
    result->a11=a*m.a11;
    result->a12=a*m.a12;
    result->a13=a*m.a13;

    result->a20=a*m.a20;
    result->a21=a*m.a21;
    result->a22=a*m.a22;
    result->a23=a*m.a23;

    result->a30=a*m.a30;
    result->a31=a*m.a31;
    result->a32=a*m.a32;
    result->a33=a*m.a33;
}

void mulScalHMatrix33Ex(HMATRIX33_PTR dest,real a)
{
    dest->a00*=a;
    dest->a01*=a;
    dest->a02*=a;
    dest->a03*=a;

    dest->a10*=a;
    dest->a11*=a;
    dest->a12*=a;
    dest->a13*=a;

    dest->a20*=a;
    dest->a21*=a;
    dest->a22*=a;
    dest->a23*=a;

    dest->a30*=a;
    dest->a31*=a;
    dest->a32*=a;
    dest->a33*=a;
}

/**Calcola il determinante della matrice 3x3 usando i cofattori*/
real det(const MATRIX33& m)
{
    return m.a00*(m.a11*m.a22-m.a21*m.a12)-m.a01*(m.a10*m.a22-m.a20*m.a12)+m.a02*(m.a10*m.a21-m.a20*m.a11);
}

real det(const MATRIX22& m) {

    return m.a00*m.a11-m.a01*m.a10;
}

/**Calcola la matrice inversa di una matrice 3x3*/
bool inverse(MATRIX33_PTR r,const MATRIX33& m)
{
    real d=det(m);

    if (d==0) return false; //non invertibile
    d=1.0/d;

    r->a00=d*(m.a11*m.a22-m.a21*m.a12);
    r->a01=d*(-m.a01*m.a22+m.a21*m.a02);
    r->a02=d*(m.a01*m.a12-m.a11*m.a02);

    r->a10=d*(-m.a10*m.a22+m.a20*m.a12);
    r->a11=d*(m.a00*m.a22-m.a20*m.a02);
    r->a12=d*(-m.a00*m.a12+m.a10*m.a02);

    r->a20=d*(m.a10*m.a21-m.a20*m.a11);
    r->a21=d*(-m.a00*m.a21+m.a20*m.a01);
    r->a22=d*(m.a00*m.a11-m.a01*m.a10);

    return true;
}

/**trsforma una matrice omogenee in una semplice 3x3*/
void toMatrix33(MATRIX33_PTR result,const HMATRIX33& m)
{
    result->a00=m.a00;
    result->a01=m.a01;
    result->a02=m.a02;

    result->a10=m.a10;
    result->a11=m.a11;
    result->a12=m.a12;

    result->a20=m.a20;
    result->a21=m.a21;
    result->a22=m.a22;
}

/**trasforma una 3x3 in una matrice in coordinate omogenee*/
void toHMatrix33(HMATRIX33_PTR result,const MATRIX33& m)
{
    result->a00=m.a00;
    result->a01=m.a01;
    result->a02=m.a02;
    result->a03=0;

    result->a10=m.a10;
    result->a11=m.a11;
    result->a12=m.a12;
    result->a13=0;

    result->a20=m.a20;
    result->a21=m.a21;
    result->a22=m.a22;
    result->a23=0;

    result->a30=result->a31=result->a32=0;
    result->a33=1;

}

/**
 * Risolve il sistema lineare 3x3 con la regola di Cramer
 * @param x vettore nel quale vengono restituiti i risultati
 * @param aij matrice dei coefficienti
 * @param d vettore termini noti
 * @return true se il sistema è risolvibile, false se ha infinite soluzioni (det(aij)=0)
 *
 */
bool solveEquationsSet3(VECTOR3D_PTR x,const MATRIX33& aij,const VECTOR3D& d) {

    double dt=det(aij);

    if (dt==0) return false; //infinite soluzioni

    if (x) {

        MATRIX33 tm;

        tm=aij;

        //sostituisce la colonna i-esima con la colonna dei termini noti
        tm.a00=d.x;
        tm.a10=d.y;
        tm.a20=d.z;

        x->x=det(tm)/dt;

        tm=aij;

        tm.a01=d.x;
        tm.a11=d.y;
        tm.a21=d.z;

        x->y=det(tm)/dt;

        tm=aij;

        tm.a02=d.x;
        tm.a12=d.y;
        tm.a22=d.z;

        x->z=det(tm)/dt;
    }

    return true;

}

/**
 * Risolve il sistema lineare 2x2 con la regola di Cramer
 * @param x vettore nel quale vengono restituiti i risultati
 * @param aij matrice dei coefficienti
 * @param d vettore termini noti
 * @return true se il sistema è risolvibile, false se ha infinite soluzioni (det(aij)=0)
 *
 */
bool solveEquationsSet2(VECTOR2D_PTR x,const MATRIX22& aij,const VECTOR2D& d) {

    double dt=det(aij);

    if (dt==0) return false; //infinite soluzioni

    if (x) {

        MATRIX22 tm;

        tm=aij;

        tm.a00=d.x;
        tm.a10=d.y;

        x->x=det(tm)/dt;

        tm=aij;

        tm.a01=d.x;
        tm.a11=d.y;       

        x->y=det(tm)/dt;
        
    }

    return true;


}

#ifdef _USENAMESPACE_
};
#endif
