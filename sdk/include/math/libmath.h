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

/**
 \file libmath.h
 *\brief pre processor definitions to compile the project as a DLL library
 */

#ifndef LIBMATH_H_INCLUDED
#define LIBMATH_H_INCLUDED
//----------------------------------------------------------------------

//definire questo flag quando si vuol esportare su DLL
#define CREATE_MATH_DLL
//#define CREATE_MATH_VB_EXPORT //funzioni esportate per Visual Basic

#ifdef CREATE_MATH_DLL
#include <windows.h>
#define MATH_EXPORT __declspec(dllexport)
#else
//compilato come eseguibile
#define MATH_EXPORT
#endif

#ifdef CREATE_MATH_VB_EXPORT
#undef MATH_EXPORT
#define MATH_EXPORT __stdcall __declspec(dllexport) //le funzioni esportabili per visual basic devono avere __stdcall
#endif

#endif // LIBMATH_H_INCLUDED
