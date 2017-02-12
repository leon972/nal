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


#include "math/mexception.h"
#include <cstdarg>
#include <cstdio>

#ifdef _USENAMESPACE_
namespace mathengine {
#endif

    using namespace std;

    //costruttore protetto usato dalla classi derivate
    MathException::MathException():err_stack()
    {
        tmp_message=0;
        id=0;
    }

    MathException::MathException(const string &msg):err_stack()
    {
        id=0;
        err_stack.push_back(msg);
    }

    MathException::MathException(const TCHAR* formatText ...):err_stack() {

        const unsigned msg_size=1024;

        tmp_message=new TCHAR[msg_size];
        memset(tmp_message,0,msg_size*sizeof(TCHAR));

        va_list		argumentPtr; //punta alla lista degli argomenti opzionali

        if (formatText == NULL) return;

        va_start(argumentPtr, formatText);

        vsprintf(tmp_message, formatText, argumentPtr);

        va_end(argumentPtr);

        err_stack.push_back(tmp_message);

        id=0;

    }

    MathException::MathException(MathException* parent,const TCHAR* formatText ...):err_stack()
    {
        const unsigned msg_size=1024;

        if (parent)
        {
            list<string>::iterator it=err_stack.begin();
            err_stack.insert(it,parent->err_stack.begin(),parent->err_stack.end());
        }

        tmp_message=new TCHAR[msg_size];

        memset(tmp_message,0,msg_size*sizeof(TCHAR));

        va_list		argumentPtr; //punta alla lista degli argomenti opzionali

        if (formatText == NULL) return;

        va_start(argumentPtr, formatText);

        vsprintf(tmp_message, formatText, argumentPtr);

        va_end(argumentPtr);

        err_stack.push_back(tmp_message);

        id=0;
    }

    int MathException::SetErrId(int err_id)
    {
        id=err_id;
        return id;
    }

    MathException::~MathException() {

        if (tmp_message)
        {
            delete[] tmp_message;
            tmp_message=0;
        }
    }

    const TCHAR* MathException::GetMessage() const {
        return err_stack.back().c_str();

    }

    int MathException::GetId() const {
        return id;
    }

    const TCHAR* MathException::GetStackTrace() const
    {
        string res="";

        if (err_stack.size()) res.reserve(30*err_stack.size());

        for (list<string>::const_reverse_iterator it=err_stack.rbegin();it!=err_stack.rend();it++)
        {
            res+=(*it);
            res+="\n";
        }

        return res.c_str();
    }

#ifdef _USENAMESPACE_
};
#endif
