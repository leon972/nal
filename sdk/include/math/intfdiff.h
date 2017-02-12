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
 \file intfdiff.h
 \brief Interpolation of non equally spaced tabulated values of a function
 with the Divided Differences method.

 */


#ifndef _INTFDIFF_INCLUDE_
#define _INTFDIFF_INCLUDE_

#include "mathdefs.h"
#include "libmath.h"
#include "ftable.h" 


namespace mathengine {

    const int _one = 1; //(used for doxygen bug)

    /**
     *Interpolation of single variable functions with the divided differeces method.
     *The tamplete parameter T in this class is the type of the function values.
     */
    template <class T> class DividedDiffInterpolation : public TFunctionX<T> {
    public:

        /**
         * Consruct the divided differences interpolation engine.
         * @param lookupTable lookup table with the known values of function. (x-y pairs)
         * @param interpDegree degree of interpolation (must be between 1 and the number of lookup tables' entries-1)
         */
        DividedDiffInterpolation(const InterpLookupTable<T>& lookupTable, int interpDegree = 2) throw (InvalidArgumentException) {

            init();
            setLookupTable(lookupTable);
            setInterpDegree(interpDegree);

        }

        /**
         * Construct the divided differences interpolation engine.
         * @param interpDegree degree of interpolation, in the degree of the polynomials used to interpolate
         * the given x-y pairs
         */
        DividedDiffInterpolation(int interpDegree = 2) {

            init();
            setInterpDegree(interpDegree);

        }

        virtual ~DividedDiffInterpolation() {

            releaseTempTable();

        }

        /**
         * Set the degree of the approximating polynomial (must be greater than zero and less than the
         * number of entries in the lookup table.
         * @param degree degree of the interpolant polynomial.
         */
        void setInterpDegree(int degree) throw (InvalidArgumentException) {

            if (degree <= 0 || (lookupTable != 0 && degree >= lookupTable->getNumEntries())) throw InvalidArgumentException();
            interpDegree = degree;
            releaseTempTable();
            //tabella delle ascisse usata per la interpolazione
            //il numero di elementi dipende dal grado di interpolazione
            py = new REAL[interpDegree + 1];
            px = new REAL[interpDegree + 1];
            //diff=new REAL[interpDegree];

        }

        /**
         * Set the lookup table
         * @param lookupTable the table with the known values of function to be interpolated
         */
        void setLookupTable(const InterpLookupTable<T>& lookupTable) {

            this->lookupTable = &lookupTable;

            //vettore delle ascisse ordinato in ordine crescente (xs[i]<xs[i+1])
            xs = lookupTable.getXs();
            ys = lookupTable.getYs();
            n = lookupTable.getNumEntries();

        }

        /**
         * Get the degree of the approximating polynomial.
         * @return 
         */
        int getInterpDegree() {
            return interpDegree;
        }

        /**
         * Compute the interpolated value of then function.
         * @param x abscissa of the point.
         * @return the interpolated value of f(x)
         */
        T getInterpY(REAL x) throw () {

            //riarrangia la tabella delle ascisse immettendo nella tabella
            //ptrx le ascisse in ordine di distanza crescente da x (ptrx[0]
            //è l'ascissa della tabella da interpolare più vicina a x

            arrangeTable(x);

            //la tabella px[i] py[i] contiene interpDegree+1 valori
            //estratti dalla tabella dei valori originali


            for (int i = 1; i <= interpDegree; i++) {

                for (int j = i; j <= interpDegree; j++) {

                    double x2 = px[j];
                    double x1 = px[i - 1];

                    py[j] = (py[i - 1]*(x2 - x) - py[j]*(x1 - x)) / (x2 - x1);
                }
            }

            return py[interpDegree];
        }

        /**
         * Computes the linear interpolation of f(x)
         * @param x abscissa of point.
         * @return the linear approximation of f(x)
         */
        T getLinearInterpY(REAL x) throw () {

            int i1 = 0;
            int i2 = n - 1;
            int i = (i2 + i1) / 2;
            REAL vx;

            if (x >= xs[i2]) {
                i1 = i2 - 1;
            } else if (x <= xs[i1]) {
                i2 = i1 + 1;
            } else {

                //ricerca binaria : determina le ascisse che definiscono
                //l'intervallo in cui si trova x
                while (i2 - i1 > 1) {

                    vx = xs[i];

                    if (x >= vx) {
                        i1 = i;
                    } else {
                        i2 = i;
                    }

                    i = (i2 + i1) / 2;
                }
            }

            //-------------------------
            //differenza divisa del primo ordine
            double x2 = xs[i2];
            double x1 = xs[i1];

            return (ys[i1]*(x2 - x) - ys[i2]*(x1 - x)) / (x2 - x1);

        }

        /**
         * Interpolated value of f(x)
         * Use this method to compute the interpolated values of the function.
         * @param x
         * @return
         */
        T operator()(REAL x) {
            return getInterpY(x);
        }

        /**
         * Utility method to compute the linear interpolated value
         * @param f1 f(x1)
         * @param f0 f(x0)
         * @param x1 abscissa of the second point
         * @param x0 abscissa of the first point
         * @param x abscissa of the point
         * @return the linear interpolated value.
         */
        static T firstDegInterp(T f1, T f0, REAL x1, REAL x0, REAL x);

    protected:

        const InterpLookupTable<T> *lookupTable;

        /**
         * Compute the first degree divided difference (f1-f0)/(x1-x0)
         * @param f1
         * @param f0
         * @param x1
         * @param x0
         * @return
         */
        T finiteDiff(T f1, T f0, REAL x1, REAL x0) {
            return (f1 - f0) / (x1 - x0);
        }

        /**
         * Arrange the values in the lookup table so that x0 is the nearest to x
         * @param x point to be interpolated
         */
        void arrangeTable(REAL x) {

            //riarrangia la tabella usando gli indici in modo da avere come x0
            //l'ascissa più vicina al punto da interpolare
            int x0 = findNearestX(x, xs, n);

            //punta alla ascissa più vicina
            py[0] = ys[x0];
            px[0] = xs[x0];

            REAL vx;
            int ifwd = x0 + 1, ifback = x0 - 1;
            int i = 1;

            while (i <= interpDegree) {

                if (ifback >= 0 && ifwd < n) {

                    if (xs[ifwd] - x <= x - xs[ifback]) {
                        px[i] = xs[ifwd];
                        py[i] = ys[ifwd];
                        ++ifwd;
                    } else {
                        px[i] = xs[ifback];
                        py[i] = xs[ifback];
                        --ifback;
                    }
                } else if (ifback >= 0) {
                    px[i] = xs[ifback];
                    py[i] = xs[ifback];
                    --ifback;
                } else if (ifwd < n) {
                    px[i] = xs[ifwd];
                    py[i] = ys[ifwd];
                    ++ifwd;
                }

                ++i;

            } //while

        }

        /**
         * Find the nearest abscissa in the vector xs to x
         * @param x
         * @param xs vector of abscissas of known points
         * @param n size of xs
         * @return the index of the nearest xs to x
         */
        int findNearestX(const REAL x, const REAL* xs, int n) {

            int i1 = 0;
            int i2 = n - 1;
            int i = (i2 + i1) / 2;
            REAL vx;

            if (x >= xs[i2]) return i2;

            else if (x <= xs[i1]) return i1;

            //ricerca binaria
            while (i2 - i1 > 1) {

                vx = xs[i];

                if (x >= vx) {
                    i1 = i;
                } else {
                    i2 = i;
                }

                i = (i2 + i1) / 2;
            }

            //i1 e i2 a questo punto contengono gli indici nella tabella di interpolazione
            //xs delle ascisse che definiscono un intervallo che contiene x
            //xs[i1]<=x<=xs[i2]
            return x - xs[i1] < xs[i2] - x ? i1 : i2;

        }

        void releaseTempTable() {

            if (py) {
                delete[] py;
                py = 0;
            }

            if (px) {
                delete[] px;
                px = 0;
            }

        }

        int interpDegree; /*!<-interpolation degree*/

        T* py; /*!<-temporary data*/
        REAL* px;

        //REAL* diff; /**tabella temporanea per il calcolo dei polinomi di interpolazione*/

        const REAL* xs; /*!<-abscissas vector*/
        const T* ys; /*!<-ordinates*/

        int n; /*!<-values count*/

    private:

        /**
         * Initialise
         */
        void init() {
            py = 0;
            px = 0;
            n = 0;
            xs = 0;
            ys = 0;
            //diff = 0;
        }
    };

};


#endif


