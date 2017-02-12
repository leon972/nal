/**************************************************
Implementa un vettore di REAL a dimensione fissa
code by L.Berti (c) 2009
***************************************************/

#include "math/nvector.h"
#include <cmath>

#ifdef _USENAMESPACE_
namespace mathengine {
    #endif

    NVector::NVector(REAL* values,size_t size)
    {
        this->values=new REAL[size];
        memcpy(this->values,values,size);
        n_components=size; //numero di componenti
        mem_size=size*sizeof(REAL); //memoria allocata

    }

    NVector::NVector(size_t size)
    {
        values=new REAL[size];
        n_components=size;
        mem_size=size*sizeof(REAL);
        memset(values,0,mem_size);
    }

    NVector::~NVector()
    {
        if (values)
        {
            delete[] values;
            values=0;
            n_components=0;
            mem_size=0;
        }
    }

    /**Resetta la matrice e ne cambia le dimensioni*/
    void NVector::Reset(size_t size,bool set_to_zero)
    {
        if (values)
        {
            delete[] values;
            values=0;
        }

        values=new REAL[size];

        n_components=size;
        mem_size=size*sizeof(REAL);

        if (set_to_zero) memset(values,0,mem_size);

    }

    REAL NVector::get_value(size_t index) const
    {
        return values[index];
    }

    const REAL* NVector::get_values() const
    {
        return values;
    }

    REAL& NVector::operator [] (size_t index)
    {
        return values[index];
    }

    REAL NVector::set_value(size_t index,REAL value)
    {
        values[index]=value;
        return value;
    }

    void NVector::set_values(const REAL* val) throw (MathException)
    {
        if (!val) throw MathException("null values pointer");

        memcpy(values,val,mem_size);
    }

    void NVector::ZeroVector()
    {
        memset(values,0,mem_size);
    }

    NVector& NVector::operator = (const NVector& v) throw (MathException)
    {
        if (&v != this)
        {
            if (v.GetSize() != n_components) throw MathException("cannot assign to vector of different size.");
            memcpy(values,v.values,mem_size);
        }

        return *this;
    }

    void NVector::add(NVector* result,const NVector& v) const
    {
        for (size_t i=0;i<n_components;++i)
        {
            values[i] += v.values[i];
        }
    }

    NVector NVector::operator + (const NVector& v) const
    {
        NVector vd(n_components);
        add(&vd,v);
        return vd;
    }

    void NVector::sub(NVector* result,const NVector& v) const
    {
        for (size_t i=0;i<n_components;++i)
        {
            values[i] -= v.values[i];
        }
    }

    NVector NVector::operator - (const NVector& v) const
    {
        NVector vd(n_components);
        add(&vd,v);
        return vd;
    }

    /**restituisce il numero di elementi*/
    size_t NVector::GetSize() const
    {
        return n_components;
    }

    REAL NVector::dot_product(const NVector& v) const
    {
        REAL r=0;
        for (size_t i=0;i<n_components;++i) r+=values[i]*v.values[i];
        return r;
    }

    REAL NVector::operator * (const NVector& v) const
    {
        REAL r=0;
        for (size_t i=0;i<n_components;++i) r+=values[i]*v.values[i];
        return r;
    }

    /**nega ogni componente*/
    void NVector::negate()
    {
         for (size_t i=0;i<n_components;++i) values[i] = -values[i];
    }

    /**normalizza (lo trasforma in versore)*/
    void NVector::normalize() throw (MathException)
    {
        REAL m=magnitude();
        m =1/m;

        if (fabs(m)<EPSILON) throw MathException("cannot normalize vector:magnitude is zero");

        for (size_t i=0;i<n_components;++i) values[i] *= m;
    }

    /**acquisisce il modulo*/
    REAL NVector::magnitude() const
    {
        REAL r;
        REAL c;
        for (size_t i=0;i<n_components;++i)
        {
            c=values[i];
            r+= c*c;
        }

        return sqrtf(r);
    }

    /**moltiplica per uno scalare*/
    void NVector::mul_scalar(REAL a)
    {
        for (size_t i=0;i<n_components;++i) values[i] *= a;
    }

     //operatori & utility
    ostream& operator<<(ostream& s,const NVector& v)
    {
        size_t sz=v.GetSize()-1;

        for (size_t i=0;i<sz;++i)
        {
            s<<v.get_value(i)<<',';
        }

        s<<v.get_value(sz);

        return s;
    }

#ifdef _USENAMESPACE_
};
#endif


