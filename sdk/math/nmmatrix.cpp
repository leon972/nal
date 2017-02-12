/***************************************
Matrice rettangolare nxm
Code by L.Berti (c) 2009
 ****************************************/




#include "math/nmmatrix.h"
#include "math/nvector.h"
#include <cmath>
#include <mem.h>
#include <iostream>

#ifdef _USENAMESPACE_
namespace mathengine {
#endif

    ///////////////////////////// NMMatrix ////////////////////////////////

    //(costruttore usato internamente dalle classi derivate)

    NMMatrix::NMMatrix() : data(0) {
    }

    NMMatrix::NMMatrix(int nrows, int ncolumns, bool set_to_zero) : data(0) {
        Reset(nrows, ncolumns, set_to_zero);
    }

    /**assegnamento*/
    NMMatrix& NMMatrix::operator =(const NMMatrix& mat) {

        if (this != &mat && mat.GetColumns() == m && mat.GetRows() == n) {

            memcpy(data, mat.data, matrix_size);
        }

        return *this;
    }

    /**Resetta la matrice e ne cambia le dimensioni*/
    void NMMatrix::Reset(int nrows, int ncolumns, bool set_to_zero = true) {
        if (data) delete[] data;
        data = 0;

        n = nrows;
        m = ncolumns;

        matrix_size = n * m * sizeof (REAL);

        data = new REAL[n * m];

        if (set_to_zero) ZeroMatrix();
    }

    NMMatrix::~NMMatrix() {
        if (data) delete[] data;
        data = 0;
    }

    void NMMatrix::ZeroMatrix() {
        memset(data, 0, matrix_size);
    }

    int NMMatrix::GetRows() const {
        return n;
    }

    int NMMatrix::GetColumns() const {
        return m;
    }

    /**Restituisce la dimensione in byte dei dati della mattice*/
    size_t NMMatrix::GetDataSize() const {
        return matrix_size;
    }

    REAL NMMatrix::get_value(int row, int col) const {
        return data[row * m + col];
    }

    /**Imposta un valore della matrice*/
    void NMMatrix::set_value(int row, int col, REAL value) {
        data[row * m + col] = value;
    }

    /**Imposta una intera riga*/
    void NMMatrix::set_row(int row, const REAL* values) {
        memcpy(&data[row * m], values, sizeof (REAL) * m);
    }

    const REAL* NMMatrix::get_row(int row) const {
        return &data[row * m];
    }

    /**Imposta una parte di riga*/
    void NMMatrix::set_row(int row, const REAL* values, int row_size) throw (MathException) {
        if (row_size <= 0 || row_size > m) throw MathException("invalid row size");
        memcpy(&data[row * m], values, sizeof (REAL) * row_size);
    }

    /**Imposta una colonna*/
    void NMMatrix::set_column(int col, const REAL* values) {
        for (int i = 0; i < n; ++i) data[i * m + col] = values[i];
    }

    void NMMatrix::set_column(int col, const NVector& col_data) {
        for (int i = 0; i < n; ++i) data[i * m + col] = col_data.get_value(i);
    }

    /**scambia due righe*/
    void NMMatrix::SwapRows(int row1, int row2) {
        if (row1 == row2) return;

        size_t rsize = m * sizeof (REAL);

        REAL* t = new REAL[m];

        memcpy(t, &data[row1 * m], rsize);

        memcpy(&data[row1 * m], &data[row2 * m], rsize);

        memcpy(&data[row2 * m], t, rsize);

        delete[] t;

    }

    /**scambia due colonne*/
    void NMMatrix::SwapColumns(int col1, int col2) {
        if (col1 == col2) return;

        REAL t;

        for (int i = 0; i < n; ++i) {
            t = data[i * m + col1]; //n=numero righe
            data[i * m + col1] = data[i * m + col2];
            data[i * m + col2] = t;
        }
    }

    REAL* NMMatrix::operator [] (size_t row) {
        return &data[row * m];
    }

    /**moltiplica la matrice per il vettore v e restituisce un vettore*/
    void NMMatrix::mul(NVector* res, const NVector& v) const throw (MathException) {
        if (static_cast<size_t> (m) != v.GetSize()) throw MathException("wrong vector size");
        if (static_cast<size_t> (n) != res->GetSize()) throw MathException("wrong vector size");

        const REAL* val = v.get_values();

        int i, j;
        REAL r;

        for (i = 0; i < n; ++i) {
            int im = i*m;

            r = 0;
            for (j = 0; j < m; ++j) {
                r += (data[im + j] * val[j]);
            }

            res->set_value(i, r);
        }
    }

    NVector NMMatrix::operator *(const NVector& v) const throw (MathException) {
        NVector vet(n);

        mul(&vet, v);

        return vet;
    }

    /**Moltiplica la matrice corrente per m e mette il risultato in result*/
    void NMMatrix::mul(NMMatrix* result, const NMMatrix& mat) const {
        //controlla se � possibile effettuare la moltiplicazione

        //mxn per nxp = mxp

        if (!result || mat.GetRows() != GetColumns() || result->GetRows() != GetRows() || result->GetColumns() != mat.GetColumns()) return;

        int p = mat.GetColumns();

        REAL* res = result->data;

        REAL* b = mat.data;

        REAL s;

        int k;

        for (int row = 0; row < n; ++row) {
            for (int col = 0; col < p; ++col) {
                s = 0;
                //prodotto scalare fra la riga row dalla matrice corrente e la colonna col della matrice m
                for (k = 0; k < m; ++k) s += data[row * m + k] * b[k * p + col];

                res[row * p + col] = s;
            }
        }
    }

    /**Somma la matrice corrente a m e mette il risultato in result*/
    void NMMatrix::add(NMMatrix* result, const NMMatrix& mat) const {
        if (!result || result->GetRows() != n || result->GetColumns() != m || mat.GetRows() != m || mat.GetColumns() != n) return;

        REAL* res = result->data;
        REAL* b = mat.data;
        size_t sz=m*n;

        for (size_t i=0;i<sz;i++) {
            res[i]=data[i]+b[i];
        }
       
    }

    /**Sottrae m dalla matrice corrente e mette il risultato in result*/
    void NMMatrix::sub(NMMatrix* result, const NMMatrix& mat) const {
        if (!result || result->GetRows() != n || result->GetColumns() != m || mat.GetRows() != m || mat.GetColumns() != n) return;

        REAL* res = result->data;
        REAL* b = mat.data;
        size_t offs;

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; ++j) {
                offs = i * m + j;
                res[offs] = data[offs] - b[offs];
            }
        }
    }

    /**Moltiplica la matrice corrente per m e mette il risultato in result*/
    //non calola l'elemento 0,0 da controllare

    NMMatrix NMMatrix::operator *(const NMMatrix& mat) const {

        NMMatrix res(n, mat.GetColumns(), true);

        mul(&res, mat);

        return res;
    }

    NMMatrix operator +(const NMMatrix& m1, const NMMatrix& m2) {

        NMMatrix res(m1.GetRows(), m1.GetColumns(), false);

        m1.add(&res, m2);

        return res;

    }

    NMMatrix NMMatrix::operator -(const NMMatrix& mat) const {

        NMMatrix res(n, m, false);

        sub(&res, mat);

        return res;

    }

    /**
     * Logical equal operator
     * @param m
     * @return true if m equals this matrix
     */
    bool NMMatrix::operator ==(const NMMatrix& m) {

        if (m.GetColumns() != GetColumns() || m.GetRows() != GetRows()) return false;

        size_t sz = GetRows() * GetColumns();

        const REAL* pdata = m.data;

        for (size_t i = 0; i < sz; i++) {
            if (data[i] != pdata[i]) return false;
        }

        return true;

    }

    /**
     * Inequality operator
     * @param m
     * @return true if m is different from this matrix
     */
    bool NMMatrix::operator !=(const NMMatrix& m) {

        return !(*this == m);

    }


    NMMatrix NMMatrix::operator *(const REAL a) const {

        size_t sz=GetRows() * GetColumns();

        NMMatrix temp(GetRows(),GetColumns(),false);

        REAL *pdata=temp.data;

        for (size_t i=0;i<sz;i++) {
            pdata[i]=data[i]*a;
        }

        return temp;
    }


    NMMatrix& NMMatrix::operator*=(const REAL a) {

        size_t sz=GetRows() * GetColumns();

        for (size_t i=0;i<sz;i++) {
            data[i]*=a;
        }

        return *this;
    }

    ////////////////////////////////// Matrix ////////////////////////////////

    Matrix::Matrix(int n, bool set_to_zero) {

        data = 0;

        Reset(n, n, set_to_zero);

    }

    Matrix::~Matrix() {
        if (data) delete[] data;
        data = 0;
    }

    void Matrix::LoadIdentity() {
        ZeroMatrix();

        for (int i = 0; i < n; ++i) data[n * i + i] = IDENTITY_VALUE;

    }

    void Matrix::GetDiagonal(REAL* diag) {
        if (diag) {
            for (int i = 0; i < n; ++i) diag[i] = data[i * n + i];
        }
    }

    Matrix& Matrix::operator = (const Matrix& m) {

        if (&m != this && m.GetRows()==GetRows()) {
            memcpy(data,m.data,matrix_size);
        }
    }

    Matrix Matrix::operator *(const NMMatrix& m) const {

        Matrix res(n, false);

        mul(&res, m);

        return res;
    }

    Matrix& Matrix::operator *=(const NMMatrix& m) {

        Matrix res(n, false);

        mul(&res, m);

        memcpy(data,res.data,matrix_size);

        return *this;

    }

    NVector Matrix::operator *(const NVector& v) const throw (MathException) {
        NVector vet(n);

        mul(&vet, v);

        return vet;

    }

    Matrix Matrix::operator * (const REAL a) const {

        size_t sz=GetRows() * GetColumns();

        Matrix temp(GetRows(),false);

        REAL *pdata=temp.data;

        for (size_t i=0;i<sz;i++) {
            pdata[i]=data[i]*a;
        }

        return temp;
    }

    Matrix& Matrix::operator *= (const REAL a) {
        
        size_t sz=GetRows() * GetColumns();

        for (size_t i=0;i<sz;i++) {
            data[i]*=a;
        }

        return *this;
    }

    Matrix Matrix::operator +(const NMMatrix& m) const {

        auto Matrix res(n, true);

        add(&res, m);

        return res;

    }

    Matrix& Matrix::operator +=(const NMMatrix& m) {

        add(this,m);

        return *this;
    }

    Matrix Matrix::operator -(const NMMatrix& m) const {
        Matrix res(n, true);

        sub(&res, m);

        return res;
    }

    Matrix& Matrix::operator -=(const NMMatrix& m) {

        Matrix res(n, false);

        sub(&res, m);

        memcpy(data,res.data,matrix_size);

        return *this;
    }

    /**
     * Logical equal operator
     * @param m
     * @return true if m equals this matrix
     */
    bool Matrix::operator ==(const NMMatrix& m) {
        return NMMatrix::operator ==(m);

    }

    /**
     * Inequality operator
     * @param m
     * @return true if m is different from this matrix
     */
    bool Matrix::operator!=(const NMMatrix& m) {
        return NMMatrix::operator!=(m);
    }


    //operatori & utility

    ostream & operator<<(ostream& s, const NMMatrix& m) {

        int cols = m.GetColumns();
        int rows = m.GetRows();

        REAL v;

        const int width = 12; //width (numero di carateri usati per ogni elemento)
        const int prec = 6; //precision
        //salva lo stato di formattazione dell'output corrente
        ios_base::fmtflags oldf = s.flags();

        s.setf(ios_base::left, ios_base::adjustfield);

        streamsize oldp = s.precision(prec); //usa 6 cifre

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                s.width(width); //usa sempre 8 caratteri
                s.fill(' '); //padding con spazio
                v = m.get_value(i, j);

                if (fabs(v) < EPSILON) v = 0;

                s << v;

            }

            s << endl;
        }

        //ripristina la precisione
        s.precision(oldp);
        //ripristina lo stato della formattazione
        s.flags(oldf);

        return s;
    }

    /**Rende la matrice simmetrica*/
    void Matrix::MakeSymmetric(MATRIX_ELEMS elems) {
        int i, j;


        if (elems == ABOVE_DIAGONAL) {
            for (i = 0; i < n; ++i) {
                for (j = 0; j < i; ++j) {
                    data[n * i + j] = data[n * j + i];
                }
            }
        } else {
            for (i = 0; i < n; ++i) {
                for (j = 0; j < i; ++j) {
                    data[n * j + i] = data[n * i + j];
                }
            }
        }
    }

    /**Rende la matrice triangolare*/
    void Matrix::MakeTriangular(MATRIX_ELEMS elems) {

        int i, j;

        if (elems == ABOVE_DIAGONAL) {
            //mantiene gli elementi sopra la diagonale
            for (i = 0; i < n; ++i) {
                for (j = 0; j < i; ++j) {
                    data[n * i + j] = 0;
                }
            }
        } else {
            //mantiene quelli sotto la diagonale
            for (i = 0; i < n; ++i) {
                for (j = 0; j < i; ++j) {
                    data[n * j + i] = 0;
                }
            }
        }
    }

    /**Rende la matrice a banda
       width=numero di bande/2 se width=0 la matrice � diagonale se n=1 tridiagonale ecc...
     */

    void Matrix::MakeDiagonal(int width) {
        if (width < 0) return;



        int i, j;

        for (i = 0; i < n; ++i) {
            //azzera sotto banda
            for (j = 0; j < i - width; j++) {
                data[n * i + j] = 0;
            }

            //azzera sopra banda
            for (j = i + 1 + width; j < n; j++) {
                data[n * i + j] = 0;
            }
        }
    }

#ifdef _USENAMESPACE_
};
#endif


