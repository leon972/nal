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
 \file nalcomplex.h
 \brief Complex number class and functions.Some operators are optimized to reduce the risk of overflow.
 * This class provides additional functions for complex numbers not found in the std::complex class.
 */

#ifndef NALCOMPLEX_H
#define	NALCOMPLEX_H

#include <cmath>
#include <complex>
#include <iostream>
#include "mathdefs.h"

#ifdef _USENAMESPACE_
namespace mathengine {
#endif

    /**
     * Complex number class.
     * Some operations in this complex class have been optimized in order
     * to avoid overflow.
     */
    template<typename T> class Complex {
    private:

        T re, im; //real and imaginary part

    public:

        /**
         * Create a complex number.
         * Some operations are optimized to avoid overflow during intermediate
         * computations.
         * @param real_part real part
         * @param imaginary_part imaginary part
         */
        Complex(const T real_part = T(), const T imaginary_part = T()) {

            re = real_part;
            im = imaginary_part;
        };

        /**
         * Set the real part
         * @param re
         */
        void real(const T& re) {
            this->re = re;
        }

        /**
         * set the imaginary part
         * @param im
         */
        void imag(const T& im) {
            this->im = im;
        }

        /**
         * Set the real and imaginary part
         * @param re
         * @param im
         */
        void set(const T& re, const T& im) {
            this->re = re;
            this->im = im;
        }

        /**
         * Set this complex in polar notation
         * @modulus modulus
         * @angle argument in radians
         * @return
         */
        void setPolar(const T& modulus, const T& angle) {

            re = modulus * cos(angle);
            im = modulus * sin(angle);

        }

        /**
         * Get the real part (const version)
         * @return the real part
         */
        const T& real() const {
            return re;
        }

        /**
         * Get the real part
         * @return
         */
        T real() {
            return re;
        }

        /**
         * Get the imaginary part (const version)
         * @return
         */
        const T& imag() const {
            return im;
        }

        /**
         * Get the imaginary part
         * @return
         */
        T imag() {
            return im;
        }

        /**
         *Assign this complex number to scalar scalar
         */
        Complex<T> & operator =(const T& scalar) {

            re = scalar;
            im = T();
            return *this;
        }

        /**
         * Assign this complex number to a complex number z
         * @param z
         * @return
         */
        Complex<T> & operator =(const Complex<T>& z) {

            re = z.re;
            im = z.im;
            return *this;
        }

        /**
         * Add a scalar to this complex number
         * @param scalar
         * @return this complex number after addition
         */
        Complex<T> & operator +=(const T& scalar) {

            re += scalar;
            return *this;
        }

        /**
         * Add the complex number z to this complex number
         * @param z
         * @return this complex number after addition
         */
        Complex<T> & operator +=(const Complex<T>& z) {

            re += z.re;
            im += z.im;
            return *this;
        }

        /**
         * Subtract a scalar from this complez number
         * @param scalar
         * @return this complex number after subtraction
         */
        Complex<T> & operator -=(const T& scalar) {

            re -= scalar;
            return *this;
        }

        /**
         * Subtract the complex number z from this complex number
         * @param z
         * @return this number after subtraction
         */
        Complex<T> & operator -=(const Complex<T>& z) {

            re -= z.re;
            im -= z.im;
            return *this;
        }

        /*Multiply this complex number by a scalar
         @return this complex number after the multiplication
         */
        Complex<T> & operator *=(const T& scalar) {

            re *= scalar;
            im *= scalar;

        }

        /**
         * Multiply this complex number by a complex number z.
         * The multiplication is made using only 3 scalar multiplication and 2 additions
         * and 3 subtraction.
         * @param z
         * @return this complex number after the multiplication
         */
        Complex<T> & operator *=(const Complex<T>& z) {

            T ac = z.re*re;
            T bd = z.im*im;

            T re1 = ac - bd;

            im = (re + im)*(z.re + z.im) - ac - bd;
            re = re1;

            return *this;

        }

        /**
         *Compute the modulus of this complex number.
         *The calculation is made in a way to prevent intermediate overflow.
         */
        T modulus() const {

            if (re == 0) {
                return fabs(im);
            } else if (im == 0) {
                return fabs(re);
            } else if (re >= im) {

                T q = im / re;
                q *= q; //q al quadrato

                return fabs(re) * sqrt(1 + q);

            } else {

                T q = re / im;
                q *= q; //q al quadrato

                return fabs(im) * sqrt(1 + q);
            }
        }

        /**
         * Returns the squared modulus
         * @return
         */
        T squaredNorm() const {
            return re * re + im*im;
        }

        /**
         * Compute the argument of this complex number
         * @param z
         * @return the argument in radians.The result is in the range [0;pi] if im>=0 and in the range (0;-pi] if im<0
         */
        T arg() const {

            T m = modulus();

            T a1 = re / m;

            if (a1 == static_cast<T> (1)) {

                return T(); //0

            } else if (a1 == static_cast<T> (-1)) {

                return M_PI;

            } else {

                T a = acos(re / m); //restituisce un angolo fra 0 e M_PI

                if (im >= 0) {
                    return a;
                } else {
                    return -a;
                }

            }
        }

        /**
         * Compute the argument of this complex number in the range [0;2*pi)
         * @return the argument in randians in the range [0;2*pi)
         */
        T arg2() const {

            T m = modulus();

            T a = acos(re / m);

            if (im >= 0) return a;
            else {
                return 2 * M_PI - a;
            }
        }

        /**
         * Same as modulus.Introduced for compatibility with std::complex
         * @return the modulus of this complex number
         */
        T abs() const {
            return modulus();
        }

        /**
         * Divide this complex number by a complex number z.
         * The calculation is made in a way to prevent intermediate overflow.
         * @param z
         * @return
         */
        Complex<T> & operator /=(const Complex<T>& z) {

            if (z.re >= z.im) {

                T dc = z.im / z.re;
                T q = z.re + z.im*dc;

                T a = re;
                T b = im;

                re = (a + b * dc) / q;
                im = (b - a * dc) / q;

            } else {

                T cd = z.re / z.im;

                T a = re;
                T b = im;

                T q = z.re * cd + z.im;

                re = (a * cd + b) / q;
                im = (b * cd - a) / q;
            }

            return *this;
        }

        /**
         *Divide this complex number by a scalar
         */
        Complex<T> operator /=(const T& scalar) {

            re /= scalar;
            im /= scalar;

        }

        /**
         * Change this complex into its conjugate
         * @return this complex
         */
        Complex<T> conjugate() {

            im = -im;
            return *this;
        }

    };

    //---------------------- operators -----------------------------------------

    /**
     * Complex number addition
     * @param z1
     * @param z2
     * @return z1+z2
     */
    template <typename T> inline Complex<T> operator +(const Complex<T>& z1, const Complex<T>& z2) {

        Complex<T> temp = z1;
        temp += z2;
        return temp;
    }

    /**
     * Add a scalar to a complex number z1
     * @param z1
     * @param scalar
     * @return z1+scalar
     */
    template<typename T> inline Complex<T> operator +(const Complex<T>& z1, const T& scalar) {
        Complex<T> temp = z1;
        temp += scalar;
        return temp;
    }

    /**
     * Add a scalar to a complex number z1
     * @param scalar
     * @param z1
     * @return z1+scalar
     */
    template<typename T> inline Complex<T> operator +(const T& scalar, const Complex<T>& z1) {
        Complex<T> temp = z1;
        temp += scalar;
        return temp;
    }

    /**
     *Unary minus
     * @return -z
     */
    template <typename T> inline Complex<T> operator -(const Complex<T>&z) {
        Complex<T> c(-z.real(), -z.imag());
        return c;
    }

    /**
     * Complex number subtraction
     * @param z1
     * @param z2
     * @return z1-z2
     */
    template <typename T> inline Complex<T> operator -(const Complex<T>& z1, const Complex<T>& z2) {

        Complex<T> temp = z1;
        temp -= z2;
        return temp;
    }

    /**
     * Subtract a scalar from a complex number.
     * @param z1
     * @param scalar
     * @return z1-scalar
     */
    template<typename T> inline Complex<T> operator -(const Complex<T>& z1, const T& scalar) {
        Complex<T> temp = z1;
        temp -= scalar;
        return temp;
    }

    /**
     * Subtract a complex from a scalar
     * @param scalar
     * @param z
     * @return scalar - z
     */
    template<typename T> inline Complex<T> operator -(const T& scalar, const Complex<T>& z) {

        Complex<T> temp(scalar - z.real(), -z.imag());

        return temp;

    }

    /**
     * Multiply two complex numbers
     * @param z1
     * @param z2
     * @return  z1*z2
     */
    template <typename T> inline Complex<T> operator *(const Complex<T>& z1, const Complex<T>& z2) {

        Complex<T> temp = z1;
        temp *= z2;
        return temp;
    }

    /**
     * Multiply a complex number by  a scalar
     * @param z1
     * @param scalar
     * @return z1*scalar
     */
    template <typename T> inline Complex<T> operator *(const Complex<T>& z1, const T& scalar) {

        Complex<T> temp = z1;
        temp *= scalar;
        return temp;

    }

    /**
     * Multiply a scalar by a complex
     * @param scalar
     * @param z
     * @return  scalar * z
     */
    template <typename T> inline Complex<T> operator *(const T& scalar, const Complex<T>& z) {

        Complex<T> temp = z;
        temp *= scalar;
        return temp;

    }

    /**
     * Complex number division
     * @param z1
     * @param z2
     * @return z1/z2
     */
    template <typename T> inline Complex<T> operator /(const Complex<T>& z1, const Complex<T>& z2) {

        Complex<T> temp = z1;
        temp /= z2;
        return temp;
    }

    /**
     * Divide a complex number by a scalar
     * @param z1
     * @param scalar
     * @return z1/scalar
     */
    template <typename T> inline Complex<T> operator /(const Complex<T>& z1, const T& scalar) {
        Complex<T> temp = z1;
        temp /= scalar;
        return temp;
    }

    /**
     * Divide a scalar by a complex
     * @param scalar
     * @param z
     * @return scalar / z
     */
    template <typename T> inline Complex<T> operator /(const T& scalar, const Complex<T>& z) {

        Complex<T> zn(scalar, T());

        return zn / z;
    }

    /**
     * Test whether z1 and z2 are equals
     * @param z1
     * @param z2
     * @return true if z1 equals z2
     */
    template <typename T> inline bool operator ==(const Complex<T>& z1, const Complex<T>& z2) {

        return z1.real() == z2.real() && z1.imag() == z2.imag();
    }

    /**
     *
     * @param z1
     * @param z2
     * @return true if z1 != z2
     */
    template <typename T> inline bool operator !=(const Complex<T>& z1, const Complex<T>& z2) {

        return z1.real() != z2.real() && z1.imag() != z2.imag();
    }

    //--------------------- functions ------------------------------------------

    //all functions for the mathengine::Complex use the prefix c_ to avoid
    //using scope (::) to access cmath corresponding functions

    /**
     * Argument of a complex number
     * @param z
     * @return the argument in radians in the range [0;pi] if im>=0 and in the range (0;-pi] if im<0
     */
    template<typename T> inline T c_arg(const Complex<T>& z) {
        return z.arg();

    }

    /**
     * Argument of a complex number
     * @param z
     * @return the argument in radians in the interval [0,2*pi)
     */
    template<typename T> inline T c_arg2(const Complex<T>& z) {
        return z.arg2();
    }

    /**
     * Cosine of a complex number z
     * @param z
     * @return cos(z)
     */
    template<typename T> inline Complex<T> c_cos(const Complex<T>& z) {

        T re = z.real();
        T im = z.imag();

        Complex<T> c(::cos(re) * ::cosh(im), -::sin(re) * ::sinh(im));
        return c;
    }

    /**
     * Arc cosine (inverse cosine) of a complex number
     * @param z
     * @return acos(z)
     */
    template<typename T> inline Complex<T> c_acos(const Complex<T>& z) {

        Complex<T> i(0, 1);

        return -i * c_log(z + c_sqrt(z * z - static_cast<T> (1)));
    }

    /**
     * Hyperbolic cosine of a complex number
     * @param z a complex number
     * @return cosh(z)
     */
    template<typename T> inline Complex<T> c_cosh(const Complex<T>& z) {

        T re = z.real();
        T im = z.imag();

        Complex<T> c(::cosh(re) * ::cos(im), ::sinh(re) * ::sin(im));

        return c;
    }

    /**
     *Inverse huperbolic cosine
     */
    template<typename T> inline Complex<T> c_acosh(const Complex<T>& z) {
        return c_log(z + c_sqrt(z - 1) * c_sqrt(z + 1));
    }

    /**
     * Sine of a complex number z
     * @param z a complex number
     *
     * @return sin(z)
     */

    template<typename T> inline Complex<T> c_sin(const Complex<T>& z) {

        T re = z.real();
        T im = z.imag();

        Complex<T> c(::sin(re) * ::cosh(im), ::cos(re) * ::sinh(im));

        return c;
    }

    /**
     * Arcsin of a complex number
     * @param z a complex number
     * @return the complex arcsin of z
     */

    template<typename T> inline Complex<T> c_asin(const Complex<T>& z) {

        Complex<T> i(0, 1);
        Complex<T> i1(0, -1);

        return i1 * c_log(i * z + c_sqrt(static_cast<T> (1) - z * z));
    }

    /**
     * Hyperbolic sine of a complex number
     * @param z
     * @return sinh(z)
     */
    template<typename T> inline Complex<T> c_sinh(const Complex<T>& z) {

        T re = z.real();
        T im = z.imag();

        Complex<T> c(::sinh(re) * ::cos(im), ::cosh(re) * ::sin(im));

        return c;
    }

    /**
     * Inverse hyperbolic sine
     * @param z
     * @return
     */
    template<typename T> inline Complex<T> c_asinh(const Complex<T>& z) {

        return c_log(z + c_sqrt(z * z + 1));
    }

    /**
     * Tangent of a complex number
     * @param z
     * @return tan(z)
     */
    template<typename T> inline Complex<T> c_tan(const Complex<T>& z) {

        return c_sin(z) / c_cos(z);
    }

    /**
     * Arctangent (inverse tangent) of a complex number
     * @param z
     * @return atan(z)
     */
    template<typename T> inline Complex<T> c_atan(const Complex<T>& z) {

        Complex<T> i(0, 1);

        Complex<T> iz = i*z;

        return (i / static_cast<T> (2))*(c_log(static_cast<T> (1) - iz) - c_log(static_cast<T> (1) + iz));
    }

    /**
     * Hyperbolic tangent of a complex number
     * @param z
     * @return tanh(z)
     */
    template<typename T> inline Complex<T> c_tanh(const Complex<T>& z) {

        return c_sinh(z) / c_cosh(z);
    }

    /**
     * Inverse hyperbolic tangent
     * @param z
     * @return
     */
    template<typename T>inline Complex<T> c_atanh(const Complex<T>& z) {

        return 0.5 * (c_log(1 + z) - c_log(1 - z));
    }

    /**
     * Cotangent of a complex number
     * @param z
     * @return cot(z)
     */
    template<typename T> inline Complex<T> c_cot(const Complex<T>& z) {

        return c_cos(z) / c_sin(z);

    }

    /**
     * Arccotangent (inverse cotangent) of a complex number
     * @param z
     * @return acot(z)
     */
    template<typename T> inline Complex<T> c_acot(const Complex<T>& z) {

        Complex<T> i(0, 1);

        Complex<T> iz = i / z;

        return (i / static_cast<T> (2))*(c_log(static_cast<T> (1) - iz) - c_log(static_cast<T> (1) + iz));
    }

    /**
     * Hyperbolic cotangent of a complex number
     * @param z
     * @return coth(z)
     */
    template<typename T> inline Complex<T> c_coth(const Complex<T>& z) {
        return c_cosh(z) / c_sinh(z);
    }

    /**
     *Inverse hyperbolic cotangent
     */
    template<typename T> inline Complex<T> c_acoth(const Complex<T>& z) {
        Complex<T> z1 = 1 / z;
        return 0.5 * (c_log(1 + z1) - c_log(1 - z1));
    }

    /**
     * Secant of a complex number
     * @param z
     * @return sec(z)
     */
    template<typename T> inline Complex<T> c_sec(const Complex<T>& z) {

        return static_cast<T> (1) / c_cos(z);
    }

    /**
     * Arcsecant (inverse secant) of a complex number
     * @param z
     * @return asec(z)
     */
    template<typename T> inline Complex<T> c_asec(const Complex<T>& z) {

        Complex<T> z1 = static_cast<T> (1) / z;

        Complex<T> i(0, 1);

        return -i * c_log(c_sqrt(z1 * z1 - static_cast<T> (1)) + z1);
    }

    /**
     *Hyperbolic secant of a complex number
     * @return sech(z)
     */
    template<typename T> inline Complex<T> c_sech(const Complex<T>& z) {
        return static_cast<T> (1) / c_cosh(z);
    }

    /**
     *Inverse hyperbolic secant
     */
    template<typename T> inline Complex<T> c_asech(const Complex<T>& z) {

        Complex<T> z1 = 1 / z;

        return c_log(c_sqrt(z1 - 1) * c_sqrt(1 + z1) + z1);
    }

    /**
     *Cosecant of a complex number
     * @return cosec(z)
     */
    template<typename T> inline Complex<T> c_cosec(const Complex<T>& z) {
        return static_cast<T> (1) / c_sin(z);
    }

    /**
     * Arccosecant (inverse cosecant) of a complex number
     * @return acosec(z)
     */
    template<typename T> inline Complex<T> c_acosec(const Complex<T>& z) {

        Complex<T> i(0, 1);

        return -i * log(c_sqrt(static_cast<T> (1) - static_cast<T> (1) / (z * z)) + i / z);
    }

    /**
     *Hyperbolic cosecant of a complex number
     * @return cosech(z)
     */
    template<typename T> inline Complex<T> c_cosech(const Complex<T>& z) {
        return static_cast<T> (1) / c_sinh(z);
    }

    /**
     * Inverse hyperbolic cosecant
     * @param z
     * @return
     */
    template<typename T>inline Complex<T> c_acosech(const Complex<T>& z) {

        Complex<T> z1 = 1 / z;

        return c_log(c_sqrt(1 + 1 / (z * z)) + z1);

    }

    /**
     *Exponentiation of a complex number
     * @return e raised to complex z
     */
    template<typename T> inline Complex<T> c_exp(const Complex<T>& z) {

        T im = z.imag();

        Complex<T> c(::cos(im), ::sin(im));

        c *= ::exp(z.real());

        return c;
    }

    /**
     * Natural logarithm of a complex number
     * @param z
     * @return log(z)
     */
    template<typename T> inline Complex<T> c_log(const Complex<T>& z) {

        Complex<T> c(::log(z.modulus()), z.arg());

        return c;
    }

    /**
     * Exponentiation
     * @param z complex base
     * @param a complex exponent
     * @return z raised to a
     */
    template <typename T> inline Complex<T> c_pow(const Complex<T>& z, const Complex<T>& a) {

        T A = z.squaredNorm();

        T c = a.real();
        T d = a.imag();

        T argab = z.arg();

        T r = ::exp(-d * argab) * ::pow(A, c * 0.5);

        T lA = ::log(A);

        Complex<T> r2(::cos(c * argab + 0.5 * d * lA), ::sin(c * argab + 0.5 * d * lA));

        r2 *= r;

        return r2;
    }

    /**
     * A real base raised to a complex exponent
     * @return base^z
     */
    template<typename T> inline Complex<T> c_pow(const T& base, const Complex<T>& z) {

        Complex<T> b(base, 0);

        return c_pow(b, z);
    }

    /**
     * z raised to a real exponent
     * @param z
     * @param exponent
     * @return z^exponent
     */
    template<typename T> inline Complex<T> c_pow(const Complex<T>& z, const T& exponent) {

        Complex<T> e1(exponent, 0);

        return c_pow(z, e1);
    }

    /**
     * Square root of a complex number
     * (there are two roots)
     * @return one of the roots of sqrt(z).(The other root is the opposite
     * of the complex returned)
     */
    template <typename T> inline Complex<T> c_sqrt(const Complex<T>& z) {

        T r = z.modulus();

        if (r == 0) {

            Complex<T> c(0, 0);

            return c;

        } else if (r == z.real()) {

            Complex<T> c(::sqrt(z.real()), 0);

            return c;

        } else {

            Complex<T> c(0, -::sqrt((r - z.real()) / static_cast<T> (2)));

            c.real(z.imag() / (2 * c.imag()));

            return c;
        }

    }

    /**
     * Calculate the 2 roots of the square root of a complex number
     * @param z the complex number
     * @param root1 out parameter, is set to the 1st root
     * @param root2 out parameter, is set to the 2nd root
     * root1=-root2
     */
    template <typename T> void c_sqrt(const Complex<T>& z, Complex<T>& root1, Complex<T>& root2) {

        root1 = sqrt(z);
        root2 = -root1;

    }

    //-------------------- utility ---------------------------------------------

    /**
     * Create a complex number using polar notation
     * @param modulus
     * @param angle argument angle in radians
     * @return the complex number with specified modulus and argument.
     */
    template<typename T> inline Complex<T> polar(const T& modulus, const T& angle) {

        Complex<T> c(modulus * cos(angle), modulus * sin(angle));
        return c;
    }

    /**
     * Get the conjugate of a complex number z
     * @param z a complex number
     * @return the conjugate of z
     */
    template<typename T> inline Complex<T> conjugate(const Complex<T>& z) {

        Complex<T> c(z.real(), -z.imag());
        return c;
    }

    /**
     * Display a xomplex number on tha output stream ostrm
     * @param ostrm
     * @param z
     * @return the output stream ostrm
     */
    template<typename T> ostream & operator <<(ostream& ostrm, const Complex<T>& z) {

        if (z.imag() < 0) {

            //parte immaginaria negativa
            ostrm << z.real() << z.imag() << 'i';

        } else {

            ostrm << z.real() << '+' << z.imag() << 'i';

        }

        return ostrm;
    }

    //istream& operator>>(istream& s,CHVector3<REAL>& vi)

    /**
     * Read a complex from an input stream.
     * The format must be x+yi where x and y are floating point numbers.
     * @param s input stream (eg. cin)
     * @param c the complex
     * @return the input stream
     */
    template<typename T> istream & operator>>(istream& s, Complex<T>& c) {

        T r, i;

        c.set(T(), T());

        if (!(s >> r)) return s;

        char ch;

        if (!(s >> ch)) return s;

        if (!(ch == '+' || ch == '-')) {
            s.setstate(std::ios_base::failbit);
            return s;
        }

        if (!(s >> i)) return s;

        if (!(s >> ch)) return s;

        if (ch != 'i') {

            s.setstate(std::ios_base::failbit);
            return s;
        }

        c.set(r, i);

        return s;
    }

    //functions pointers definitions

    /**
     *Single variable complex function
     */
    typedef Complex<double> (*COMPLEX_FUNCTION)(const Complex<double>&);

    /**
     *Two variable function
     */
    typedef Complex<double> (*COMPLEX_FUNCTION_XY)(const Complex<double>&,const Complex<double>&);



#ifdef _USENAMESPACE_
};
#endif

#endif	/* NALCOMPLEX_H */

