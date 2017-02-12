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
 * \file ftable.h
 * \brief Lookup table for single variable function.
 * This lookup table is used for function interpolation.
 * 
 */



#ifndef _FTABLE_INCLUDE_
#define _FTABLE_INCLUDE_

#include "mathdefs.h"
#include "libmath.h"


namespace mathengine {

    /**
     *Lookup table of single variable function values (scissas - ordinates pairs) not necessarily equally spaced.
     */
    template <class T> class InterpLookupTable {
    private:

        REAL *x; /*!<-abscissas*/
        T *y; /*!<-ordinates*/
        int npoints; /*!<-number of pairs*/

    public:

        /**
         * Construct a lookup table with npoints values.
         * @param npoints numer of x-y pairs
         */
        InterpLookupTable(int npoints) {

            this->npoints = npoints;
            x = new REAL[npoints];
            y = new T[npoints];

        }

        /**
         * Release allocated memory
         */
        virtual ~InterpLookupTable() {

            if (x) {
                delete[] x;
                x = 0;
            }
            if (y) {
                delete[] y;
                y = 0;
            }

        }

        /**
         * Set a pair abscissa - ordinate
         * @param index index of point (from 0 to npoints)
         * @param x abscissa
         * @param y ordinate
         */
        void setValue(int index, REAL x, T y) {

            if (index >= 0 && index < npoints) {
                this->x[index] = x;
                this->y[index] = y;
            }

        }

        /**
         * Get the ordinate (value) of point with index index
         * @param index index of point
         * @return  the ordinate
         */
        T getY(int index) const {

            if (index >= 0 && index < npoints) {
                return y[index];
            }

        }

        /**
         * Get the abscissa of point with index index
         * @param index index of point
         * @return the abscissa of point with index index (from 0 to npoints)
         */
        REAL getX(int index) const {

            if (index >= 0 && index < npoints) {
                return x[index];
            }

        }

        /**
         * Set the ordinate of point index
         * @param index
         * @param valuey the ordinate (the value of the function of point index)
         */
        void setY(int index, T valuey) {

            if (index >= 0 && index < npoints) {
                y[index] = valuey;
            }

        }

        /**
         * Set the abscissa of point index.
         * The abscissas must be in ascent order.
         * @param index index of point (from 0 to npoints)
         * @param valuex value of the abscissa
         * @throw InvalidArgumentException if the index is out of range or the absissa in invalid.
         */
        void setX(int index, REAL valuex) throw (InvalidArgumentException) {
            if (index >= 0 && index < npoints) {

                if (index > 0 && !(x[index - 1] < valuex)) {
                    throw InvalidArgumentException();
                }
                x[index] = valuex;

            } else throw InvalidArgumentException();

        }

        /**
         * Sets the abscissa and ordinate for point with index index
         * @param index
         * @param x abscissa
         * @param y known value of function for x y=f(x)
         */
        void setXY(int index, REAL x, T y) throw (InvalidArgumentException) {

            if (index >= 0 && index < npoints) {

                if (index > 0 && !(this->x[index - 1] < x)) {
                    throw InvalidArgumentException();
                }

                this->x[index] = x;
                this->y[index] = y;

            } else throw InvalidArgumentException();

        }

        /**
         * Get the abscissas vector.
         * @return
         */
        const REAL* getXs() const {
            return x;

        }

        /**
         * Get the ordinates.
         * @return
         */
        const T* getYs() const {
            return y;
        }

        /**
         * Get the number of x-y pairs
         * @return the number of points in the table.
         */
        int getNumEntries() const {

            return npoints;

        }

    };

};


#endif



