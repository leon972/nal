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
 Math engine - trasformazioni 3D , rotazioni,proiezioni ecc...

 ****************************************************************/

#include "math/vectors.h"
#include "math/matrix.h"
#include "math/3dtransf.h"

#ifdef _USENAMESPACE_
namespace mathengine
{
#endif

//Ruota il vettore v intorno ad un versore (formula di Rodrigues')
//(nota: usare la matrice di rotazione quando si richiede maggiore velocit�)
HVECTOR3D rotateVector(const HVECTOR3D& v,const HVECTOR3D& vers,real theta)
{
    real ct=cosr(theta);
    real st=sinr(theta);
    HVECTOR3D t;

    //v*cos(theta)+vers x v*sin(theta)+vers.v*(1-cos(theta)) vers
    t=sumHVector3D(sumHVector3D(mulHVector3D(v,ct),crossProduct3H(vers,mulHVector3D(v,st))),mulHVector3D(vers,dotProduct3H(vers,v)*(1-ct)));
    t.w=1;

    return t;
}

//TRASFORMAZIONI
//restituisce l'angolo in radianti della longitudine cio� l'angolo formato dal piano che contiene il versore
//ed � perpendicolare al piano XY e il piano XZ
real getLongitudeZ(const HVECTOR3D& vers)
{
    return acosr(vers.x/sqrtr(vers.x*vers.x+vers.y*vers.y));
}

//restituisce l'angolo in radianti fra il versore vers e il piano XY
//NB vers ha modulo 1 (� un versore)
real getAzimutZ(const HVECTOR3D& vers)
{
    real len=vers.x*vers.x+vers.y*vers.y+vers.z*vers.z;
    return acosr(sqrtr((vers.x*vers.x+vers.y*vers.y)/len));
}

//Ottiene la matrice che esegue la rotazione di theta radianti intorno al versore vers
//(formula di Rodrigues)
HMATRIX33 getRotMatrix(const HVECTOR3D& vers,real theta)
{
    auto HMATRIX33 r;

    real st=sinr(theta);
    real ct=cosr(theta);
    real ct1=1-ct;

    real xy=vers.x*vers.y;
    real xz=vers.x*vers.z;
    real yz=vers.y*vers.z;

    real x2=vers.x*vers.x;
    real y2=vers.y*vers.y;
    real z2=vers.z*vers.z;

    r.a30=r.a31=r.a32=0;
    r.a33=1; //coord-> omogenea

    r.a00=ct+ct1*x2;
    r.a01=ct1*xy-st*vers.z;
    r.a02=ct1*xz+st*vers.y;
    r.a03=0;

    r.a10=ct1*xy+st*vers.z;
    r.a11=ct+ct1*y2;
    r.a12=ct1*yz-st*vers.x;
    r.a13=0;

    r.a20=ct1*xz-st*vers.y;
    r.a21=ct1*yz+st*vers.x;
    r.a22=ct+ct1*z2;
    r.a23=0;

    return r;
}

//Ottiene la matrice che esegue la rotazione di theta radianti intorno al versore vers
//(formula di Rodrigues)
void getRotMatrixEx(HMATRIX33_PTR r,const HVECTOR3D& vers,real theta)
{
    r->a30=r->a31=r->a32=0;
    r->a33=1; //coord-> omogenea

    real st=sinr(theta);
    real ct=cosr(theta);
    real ct1=1-ct;

    real xy=vers.x*vers.y;
    real xz=vers.x*vers.z;
    real yz=vers.y*vers.z;

    real x2=vers.x*vers.x;
    real y2=vers.y*vers.y;
    real z2=vers.z*vers.z;

    r->a00=ct+ct1*x2;
    r->a01=ct1*xy-st*vers.z;
    r->a02=ct1*xz+st*vers.y;
    r->a03=0;

    r->a10=ct1*xy+st*vers.z;
    r->a11=ct+ct1*y2;
    r->a12=ct1*yz-st*vers.x;
    r->a13=0;

    r->a20=ct1*xz-st*vers.y;
    r->a21=ct1*yz+st*vers.x;
    r->a22=ct+ct1*z2;
    r->a23=0;
}

//matrici ortonormali di rotazione
HMATRIX33 getRotZMatrix(real theta)
{
    auto HMATRIX33 r;

    real ct=cosr(theta);
    real st=sinr(theta);

    r.a00=ct;r.a01=st;r.a02=0;r.a03=0;
    r.a10=-st;r.a11=ct;r.a12=0;r.a13=0;
    r.a20=r.a21=0;r.a22=1;r.a23=0;
    r.a30=r.a31=r.a32=0;r.a33=1;

    return r;
}

HMATRIX33 getRotYMatrix(real theta)
{
    auto HMATRIX33 r;

    real ct=cosr(theta);
    real st=sinr(theta);

    r.a00=ct;r.a01=0;r.a02=-st;r.a03=0;
    r.a10=0;r.a11=1;r.a12=0;r.a13=0;
    r.a20=st;r.a21=0;r.a22=ct;r.a23=0;
    r.a30=r.a31=r.a32=0;r.a33=1;

    return r;
}

HMATRIX33 getRotXMatrix(real theta)
{
    auto HMATRIX33 r;

    real ct=cosr(theta);
    real st=sinr(theta);

    r.a00=1;r.a01=0;r.a02=0;r.a03=0;
    r.a10=0;r.a11=ct;r.a12=st;r.a13=0;
    r.a20=0;r.a21=-st;r.a22=ct;r.a23=0;
    r.a30=r.a31=r.a32=0;r.a33=1;

    return r;
}

//da verificare
HMATRIX33 getEulerMatrix(const real psi,const real theta,const real phi)
{
    auto HMATRIX33 r;

    real cp=cosr(psi);
    real sp=sinr(psi);

    real ct=cosr(theta);
    real st=sinr(theta);

    real cf=cosr(phi);
    real sf=sinr(phi);

    r.a00=cf*cp-sf*ct*sp;
    r.a01=cf*sp+sf*ct*cp;
    r.a02=sf*st;
    r.a03=0;

    r.a10=-sf*cp-cf*ct*sp;
    r.a11=-sf*sp+cf*ct*cp;
    r.a12=cf*st;
    r.a13=0;

    r.a20=st*sp;
    r.a21=-st*cp;
    r.a22=ct;
    r.a23=0;

    r.a30=0;
    r.a31=0;
    r.a32=0;
    r.a33=1;

    return r;
}

void getEulerMatrixEx(HMATRIX33* m,const real psi,const real theta,const real phi)
{
    real cp=cosr(psi);
    real sp=sinr(psi);

    real ct=cosr(theta);
    real st=sinr(theta);

    real cf=cosr(phi);
    real sf=sinr(phi);

    m->a00=cf*cp-sf*ct*sp;
    m->a01=cf*sp+sf*ct*cp;
    m->a02=sf*st;
    m->a03=0;

    m->a10=-sf*cp-cf*ct*sp;
    m->a11=-sf*sp+cf*ct*cp;
    m->a12=cf*st;
    m->a13=0;

    m->a20=st*sp;
    m->a21=-st*cp;
    m->a22=ct;
    m->a23=0;

    m->a30=0;
    m->a31=0;
    m->a32=0;
    m->a33=1;
}

//Ottiene la matrice di traslazione (per traslare un vettore occorre moltiplicare Mxv)
HMATRIX33 getTranslationMatrix(const real tx,const real ty,const real tz)
{
    auto HMATRIX33 r;

    r.a00=1;
    r.a01=0;
    r.a02=0;
    r.a03=tx;

    r.a10=0;
    r.a11=1;
    r.a12=0;
    r.a13=ty;

    r.a20=0;
    r.a21=0;
    r.a22=1;
    r.a23=tz;

    r.a30=r.a31=r.a32=0;
    r.a33=1;

    return r;

}

void getTranslationMatrixEx(HMATRIX33* m,const real tx,const real ty,const real tz)
{
    m->a00=1;
    m->a01=0;
    m->a02=0;
    m->a03=tx;

    m->a10=0;
    m->a11=1;
    m->a12=0;
    m->a13=ty;

    m->a20=0;
    m->a21=0;
    m->a22=1;
    m->a23=tz;

    m->a30=m->a31=m->a32=0;
    m->a33=1;
}

HMATRIX33 getScaleMatrix(const real sx,const real sy,const real sz)
{
    auto HMATRIX33 r;

    r.a00=sx;
    r.a01=0;
    r.a02=0;
    r.a03=0;

    r.a10=0;
    r.a11=sy;
    r.a12=0;
    r.a13=0;

    r.a20=0;
    r.a21=0;
    r.a22=sz;
    r.a23=0;

    r.a30=0;
    r.a31=0;
    r.a32=0;
    r.a33=1;

    return r;

}

void getScaleMatrixEx(HMATRIX33* m,const real sx,const real sy,const real sz)
{
    m->a00=sx;
    m->a01=0;
    m->a02=0;
    m->a03=0;

    m->a10=0;
    m->a11=sy;
    m->a12=0;
    m->a13=0;

    m->a20=0;
    m->a21=0;
    m->a22=sz;
    m->a23=0;

    m->a30=0;
    m->a31=0;
    m->a32=0;
    m->a33=1;
}

#ifdef _USENAMESPACE_
};
#endif


