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
 * \file mathdefs.h
 * \brief Math constants and global definitions.
 
 */

#ifndef MATHDEFS_H_INCLUDED
#define MATHDEFS_H_INCLUDED

#define _USEDEFAULT_

/**Default settings*/
#ifdef _USEDEFAULT_
#undef _USEDOUBLE_
#define _USENAMESPACE_
#endif

/**
 *Deg to Rad conversion factor
 */
#define DEG_TO_RAD_CONV 0.01745329251994329577

/**
 *Rad to Deg conversion factor
 */
#define RAD_TO_DEG_CONV 57.2957795130823208768

/**
 *Degrees to radians conversion
 */
#define DEG_TO_RAD(x) x*DEG_TO_RAD_CONV

/**
 *Radians to degrees conversion
 */
#define RAD_TO_DEG(x) x*RAD_TO_DEG_CONV

#define _USENAMESPACE_

#ifdef _USENAMESPACE_
namespace mathengine
{
#endif

/**
 * Real number
 */
#define REAL double

const double EPSILON=0.0000000001;

/**
 *Single variable function pointer  y=f(x)
 */
typedef REAL (*FUNCT_X)(REAL);

/**
 *Double variable function pointer  z=f(x,y)
 */
typedef REAL (*FUNCT_XY)(REAL,REAL);

/**
 *Triple variable function pointer v=f(x,y,z)
 */
typedef REAL (*FUNCT_XYZ)(REAL,REAL,REAL);

/**
 *\brief Division by zero exception
 *Thrown when a division by zero occurs. 
 */
class DivideByZeroVideoException {};

/**
 *\brief Generic invalid argument exception
 *Thrown when the arguments are not valid.
 */
class InvalidArgumentException {};

/**
 *\brief No solution exception.
 * Thrown when a system or equation cannot be solved.
 * 
 */
class NoSolutionException {};

/**
 *Functor interface y=f(x)
 */
class IFunctionX
{
    public:
    virtual REAL operator()(REAL x)=0;
};

/**
 *Functor interface z=f(x,y)
 */
class IFunctionXY
{
    public:
    virtual REAL operator()(REAL x,REAL y)=0;
};

/**
 *Functor interface v=f(x,y,z)
 */
class IFunctionXYZ
{
    public:
    virtual REAL operator()(REAL x,REAL y,REAL z)=0;
};


/**
 *Template  based functor interface y=f(x)
 */
template <class T> class TFunctionX
{
    public:
    virtual T operator()(T x)=0;
};

/**
 *Template  based functor z=f(x,y)
 */
template <class T> class TFunctionXY
{
    public:
    virtual T operator()(T x,T y)=0;
};

/**
 *Template  based functor v=f(x,y,z)
 */
template <class T> class TFunctionXYZ
{
    public:
    virtual T operator()(T x,T y,T z)=0;
};


#ifdef _USENAMESPACE_
};
#endif

#endif // MATHDEFS_H_INCLUDED
