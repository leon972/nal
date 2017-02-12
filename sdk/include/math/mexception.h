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
 \file mexception.h
 \brief Generic math exception class.
 */

#ifndef MATHEXCEPTION_H_INCLUDED
#define MATHEXCEPTION_H_INCLUDED

#include "mathdefs.h"
#include "libmath.h"
#include <list>
#include <string>

#ifdef _USENAMESPACE_
namespace mathengine {
#endif

using namespace std;

/**
 Generic math exception 
 */
class MATH_EXPORT MathException
{
    private:

    list<string> err_stack;
    int id;
    TCHAR* tmp_message;

    public:

    MathException(const string &msg);
    MathException(const TCHAR* formatText ...);
    MathException(MathException* parent,const TCHAR* formatText ...);
   // VideoException(string msg,VideoException* parent);//ambigua
    int SetErrId(int err_id);
    MathException();
    virtual ~MathException();
    const TCHAR* GetMessage() const;
    int GetId() const;
    const TCHAR* GetStackTrace() const;

};


#ifdef _USENAMESPACE_
};
#endif

#endif // MVideoException_H_INCLUDED
