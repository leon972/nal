/**************************************************
Implementa un vettore di REAL a dimensione fissa
code by L.Berti (c) 2009
***************************************************/
/*!
 \file nvector.h
 \brief N components vector class. 
 
 */


#ifndef NVECTOR_H_INCLUDED
#define NVECTOR_H_INCLUDED

#include "mathdefs.h"
#include "mexception.h"
#include "libmath.h"
#include <stddef.h>
#include <iostream>

namespace mathengine {

    using namespace std;

    /**
     * Numeric vector class with n real components.
     *
     */
    class MATH_EXPORT NVector
    {
        protected:

            /**
             *The data array used to store the vector components.
             */
        REAL* values;

        /**
         *Number of components of the vector.
         */
        size_t n_components;

        /**
         * Size in bytes of the data used to store the vector.
         */
        size_t mem_size; 

        public:

        /**
         * Constructs a vector specifing its components.
         * @param values an array with the components of the vector.
         * @param size number of components.
         */
        NVector(REAL* values,size_t size);

        /**
         * Constructs a vector with size components.
         * @param size number of components.
         */
        NVector(size_t size);

        /**
         * Memory deallocation.
         */
        virtual ~NVector();

        /**
         * Recreate the vector. All existing components are erased.
         * @param size new vector size (number of components)
         * @param set_to_zero if set, all components are set to zero.
         */
        void Reset(size_t size,bool set_to_zero);

        /**
         * Get a component
         * @param index index of component, the first component has index zero.
         * @return the value of the component.
         */
        REAL get_value(size_t index) const;

        /**
         * Access a component.
         * @param index index of component, the first component has index zero.
         * @return the value of the component.
         * \warning this method doesn't perform bounds check form optimization reasons.
         */
        REAL& operator [] (size_t index);

        /**
         * Set a component.
         * @param index index of component, the first component has index zero.
         * @param value the value of the component.
         * @return the value of the component.
         * \warning this method doesn't perform bounds check form optimization reasons.
         */
        REAL set_value(size_t index,REAL value);

        /**
         * Set all components to zero.
         */
        void ZeroVector();

        /**
         * Set all components using an array.
         * @param val an array with the components' values. Must have at aleast the same size of this vector.
         */
        void set_values(const REAL* val) throw (MathException);

        /**
         * Get the values of the components.
         * @return an unmodifiable array with the values of this vector.
         */
        const REAL* get_values() const;

        /**
         * Assignment operator.
         * @param v
         * @return
         */
        NVector& operator = (const NVector& v) throw (MathException);

        /**
         * Vector addition. Adds to this vector the vector v and places the
         * result in result.
         * @param result the sum of this vector with v.
         * @param v the vector to be added. Must have at least the same size of this vector.
         */
        void add(NVector* result,const NVector& v) const;

        /**
         * Vector addition.
         * @param v the vector to be added. Must have at least the same size of this vector.
         * @return  the sum of this vector and v.
         */
        NVector operator + (const NVector& v) const;

        /**
         * Vector subctraction.
         * @param result the difference of this vector and v.
         * @param v the vector to be subctracted.
         */
        void sub(NVector* result,const NVector& v) const;

        /**
         * Vector subctraction.
         * @param v the vector to be subctracted.
         * @return the difference of this vector and v.
         */
        NVector operator - (const NVector& v) const;

        /**
         * Get the number of components.
         * @return the number of components.
         */
        size_t GetSize() const;

        /**
         * Dot (scalar) product.
         * @param v a vector with the same size of this vector.
         * @return the dot product of this vector with v.
         */
        REAL dot_product(const NVector& v) const;

        /**
         * Dot (scalar) product.
         * @param v a vector with the same size of this vector.
         * @return
         */
        REAL operator * (const NVector& v) const;

        /**nega ogni componente*/
        void negate();

        /**normalizza (lo trasforma in versore)*/
        void normalize() throw (MathException);

        /**acquisisce il modulo*/
        REAL magnitude() const;

        /**moltiplica per uno scalare*/
        void mul_scalar(REAL a);

    };

    //operatori & utility
    ostream& operator<<(ostream& s,const NVector& v);

};

#endif // NVECTOR_H_INCLUDED
