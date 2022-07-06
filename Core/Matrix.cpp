/*!
 * Contains all the code for handling matrices
 * which supports basic error checking / handling 
 * through the use of exceptions
 */
#include "stdafx.h"
#include "matrix.h"

using namespace std;

namespace Common {

Matrix::Matrix()
{ 
    initialize(1,1);
    Fill(0.0);
}

Matrix::Matrix(size_t rows, size_t cols)
{
    initialize(rows, cols);
    Fill(0.0);
}

Matrix::Matrix(const Matrix& src)
{
    copyHelper(src);
}

Matrix::~Matrix()
{
    deleteHelper();
}

Matrix Matrix::Unit(double v)
{
    Matrix m(1,1);
    m(1,1) = v;
    return m;
}

Matrix Matrix::Identity(size_t size)
{
    Matrix m(size,size);
    for (size_t i=0; i<size; i++) { m(i,i) = 1.0; }
    return m;
}

Matrix Matrix::FromVector(const vector<double>& source, bool vertical)
{
    if (vertical) 
    {
        Matrix m(source.size(), 1);
        for (size_t i=0; i<m.Rows(); i++) { m(i,0) = source[i]; }
        return m;
    }
    else
    {
        Matrix m(1, source.size());
        for (size_t i=0; i<m.Cols(); i++) { m(0,i) = source[i]; }
        return m;
    }
}

vector<double> Matrix::ToVector(size_t index, bool vertical) const
{
    if (vertical) 
    {
        vector<double> v(_rows);
        for (size_t i=0; i<_rows; i++) { v[i] = (*this)(i,index); }
        return v;
    }
    else 
    {
        vector<double> v(_cols);
        for (size_t i=0; i<_cols; i++) { v[i] = (*this)(index,i); }
        return v;
    }
}

Matrix& Matrix::operator=(const Matrix& rhs) 
{
    if (this == &rhs)
    {
        return (*this);
    }
    deleteHelper();
    copyHelper(rhs);
    return *this;
}

const Matrix operator+(const Matrix& lhs, const Matrix& rhs) 
{
    if ((lhs._rows != rhs._rows) || (lhs._cols != rhs._cols)) {throw domain_error("Matrix dimensions do not agree");}
    
    Matrix R(lhs._rows, lhs._cols);
    for (size_t i=0; i<lhs._rows; i++) {
        for (size_t j=0; j<lhs._cols; j++) {
            R(i,j) = lhs(i,j) + rhs(i,j);
        }
    }
    return R;
}

const Matrix operator-(const Matrix& lhs, const Matrix& rhs)
{
    if ((lhs._rows != rhs._rows) || (lhs._cols != rhs._cols)) {throw domain_error("Matrix dimensions do not agree");}
    
    Matrix R(lhs._rows, lhs._cols);
    for (size_t i=0; i<lhs._rows; i++) {
        for (size_t j=0; j<lhs._cols; j++) {
            R(i,j) = lhs(i,j) - rhs(i,j);
        }
    }
    return R;
}

const Matrix operator*(const Matrix& lhs, const Matrix& rhs) 
{
    if (lhs._cols != rhs._rows) {throw domain_error("Internal matrix dimensions do not agree");}
    
    Matrix R(lhs._rows, rhs._cols);

    for (size_t i=0; i<lhs._rows; i++) {
        for (size_t j=0; j<rhs._cols; j++) {
            for (size_t k=0; k<lhs._cols; k++)
            {
                R(i,j)+=lhs(i,k) * rhs(k,j);
            }
        }
    }
    return R;
}

const Matrix operator*(const Matrix& lhs, const double scalar)
{
    Matrix R(lhs._rows, lhs._cols);
    for (size_t i=0; i<lhs._rows; i++) {
        for (size_t j=0; j<lhs._cols; j++) {
            R(i,j) = lhs(i,j) * scalar;
        }
    }
    return R;
}
const Matrix operator*(const double scalar,const Matrix& rhs)
{ 
     return rhs*scalar; 
};

const Matrix operator/(const Matrix& lhs, const double scalar)
{
    if (scalar == 0.0) { throw invalid_argument("Divide by zero"); }
    return lhs * (1/scalar);
}

double Matrix::operator() (size_t row, size_t col) const
{
    if ((row>=_rows) || (col>=_cols)) {
        throw out_of_range("Requested matrix element out of range"); 
    }
    return _M[row][col];
}

double& Matrix::operator() (size_t row, size_t col)
{
    if ((row>=_rows) || (col>=_cols)) 
    {
        throw out_of_range("Requested matrix element out of range"); 
    }
    return _M[row][col];
}

size_t Matrix::Rows() const
{
    return _rows;
}

size_t Matrix::Cols() const
{
    return _cols;
}

Matrix Matrix::Transpose() const
{
    Matrix T(_cols, _rows);
    
    for (size_t i=0; i<_cols; i++) {
        for (size_t j=0; j<_rows; j++) {
            T(i,j) = (*this)(j,i);
        }
    }
    return T;
}

Matrix Matrix::Normals() const
{
    return Transpose()*(*this);
}

Matrix Matrix::removeRowCol(size_t row, size_t col) const
{
    int k;
    int n;

    if ((row>_rows) || (col>_cols)) { throw out_of_range(""); }
    Matrix R(_rows-1, _cols-1);
    for (size_t i = 0; i<_rows; i++) {
        if (i==row) { continue;}
        for (size_t j = 0; j<_cols; j++) {
            if (j==col) { continue;}
            if (i<=row) {k=i;} else {k=i-1;}
            if (j<=col) {n=j;} else {n=j-1;}
            R(k,n) = (*this)(i,j);
        }
    }
    return R;
}

Matrix Matrix::AugmentR(const Matrix &M) const
{
    // augment the matrix
    // make a new one
    // with the max of the rows
    // and the sum of the columns
    Matrix R(max(_rows,M._rows),_cols+M._cols);
    for (size_t i=0; i<_rows; i++) {
        for (size_t j=0; j<_cols; j++) {
            R(i,j) = (*this)(i,j);
        }
    }
    for (size_t i=0; i<M._rows; i++) {
        for (size_t j=0; j<M._cols; j++) {
            R(i,j+_cols) = M(i,j);
        }
    }
    return R;
}

Matrix Matrix::Row(size_t r) const
{
    if (r>_rows) { throw out_of_range(""); }
    Matrix S(1,_cols);
    for(size_t i=0; i<_cols; i++) {
        S(0,i) = (*this)(r,i);
    }
    return S;
}

Matrix Matrix::Column(size_t c) const
{
    if (c>_cols) { throw out_of_range(""); }
    Matrix S(_rows,1);
    for(size_t i=0; i<_rows; i++) {
        S(i,0) = (*this)(i,c);
    }
    return S;
}

Matrix Matrix::Extract(size_t top_r, size_t top_c, size_t bot_r, size_t bot_c) const
{
    if ((bot_r>_rows) || (bot_c>_cols)) { throw out_of_range(""); }
    Matrix R(bot_r - top_r+1,bot_c - top_c+1);
    for (size_t i = top_r; i<=bot_r; i++) {
        for (size_t j = top_c; j<=bot_c; j++) {
            R(i-top_r,j-top_c) = (*this)(i,j);
        }
    }
    return R;
}

double Matrix::Det() const
{
    // cf http://mathworld.wolfram.com/Condensation.html

    if (_rows!=_cols) {throw logic_error("Matrix must be square");}
    auto R = (*this);
    auto K1 = (*this);
    auto K2 = (*this);

    try {
        for (size_t k=_rows-1; k>=1; k--) {
            if (k<=_rows-2) {
                K1 = K1.Extract(2,2,2+k-1,2+k-1); 
            }
            K2 = R;

            for (size_t i = 0; i<k; i++) { // rows
                for (size_t j = 0; j<k; j++) { // cols
                    R(i,j) = R(i,j)*R(i+1,j+1) - R(i+1,j)*R(i,j+1);
                }
            }
            // can shave off the R now
            R = R.removeRowCol(R._rows,R._cols);

            if (k<=_rows-2) {
                // do the necessary divisions
                // get the divison matrix, it's k x k
                for (size_t s=0; s<k; s++) {
                    for (size_t t=0; t<k; t++) {
                        R(s,t)/=K1(s,t);
                    }
                }    
            }
            K1=K2;
        }
        return R(0,0);
    }
    catch (...) 
    {
        throw logic_error("Condensation method failed. Use DetSlow()");
    }
}

double Matrix::DetSlow() const
{
    // use expansion by first row
    // (basic algorithm!)
    if (_rows!=_cols) {throw logic_error("Matrix must be square");}
    auto adder = 1; // toggles
    auto det = 0.0;
    if (_rows == 2) {
        det = (*this)(0,0)*(*this)(1,1) - (*this)(0,1)*(*this)(1,0);
    }
    else {
        for (size_t i=0; i<_cols; i++) {
            det=det+(0,i)*adder*(removeRowCol(0,i).DetSlow());
            adder*=-1;
        }
    }
    if (det==0) {throw logic_error("Singular matrix");}
    return(det);
}

double Matrix::Trace() const
{
    if (_rows!=_cols) {throw logic_error("Matrix must be square");}
    auto sum = 0.0;
    for (size_t i = 0; i<_rows; i++) {
        sum+=(*this)(i,i);
    }
    return sum;
}

double Matrix::l2Norm(size_t col) const
{
    auto sum = 0.0;
    for (size_t i=0; i<_rows; ++i)
    {
        sum += (*this)(i,col) * (*this)(i,col);
    }
    return sqrt(sum);

}
double Matrix::infinityNorm(size_t col) const
{
    double maximum = 0.0;
    for (size_t i=0; i<_rows; ++i)
    {
        maximum = max (abs((*this)(i,col)), maximum);
    }
    return maximum;
}

void Matrix::Fill(double d)
{
    for (size_t i=0; i<_rows; i++) {
        for (size_t j=0; j<_cols; j++) {
            (*this)(i,j) = d;
        }
    }
}

PCAResult Matrix::PrincipalComponentAnalysis (double coverageFactor) const 
{
    // Currently using Power method
    auto trace = Trace();
    Matrix source((*this));
    vector<double> eigenvalues;
    vector<Matrix> eigenvectors;
    vector<double> norms;
    Matrix eigenvector1(source._cols,1);
    Matrix eigenvector2(source._cols,1);
    auto current = &eigenvector1;
    auto next = &eigenvector2;
    auto tol = 1e-11;
    auto explained = 0.0;
    auto factorId = 0;
    // stop looking once we have breached the limit
    // NOTE: numerical limits mean checking explained == 1 fails, hence using of ? operator
    while (coverageFactor >= 1 ? factorId!=source._cols : explained <= coverageFactor)
    {
        current->Fill(0.01); 
        next->Fill(0.0);
        do {
            *next = source * *current;
            *next = *next / next->infinityNorm(0);    
            swap(current, next);
        } 
        while ((*next-*current).l2Norm(0) > tol);
        auto eigenvalue = (source * *current).infinityNorm(0); 
        // normalise for removal
        auto norm = current->l2Norm(0);
        auto eigenvector(*current / norm);
        // remove this from the covariance matrix
        source = source - eigenvalue * (eigenvector * eigenvector.Transpose());
        // store the answer, after putting the eigenvalue in the first column
        eigenvalues.push_back(eigenvalue);
        eigenvectors.push_back(eigenvector);
        norms.push_back(norm);
        // update coverage and tracker
        explained += eigenvalue / trace;
        factorId++;
    }

    PCAResult pcar;
    pcar.Coverage = explained;
    pcar.Eigenvalues = eigenvalues;
    pcar.Norms = norms;
    pcar.Eigenvectors = eigenvectors;
    return pcar;

}

Matrix Matrix::fit(const Matrix& observed) const
{
    auto L = Normals().Cholesky();
    auto b = Transpose() * observed; 

    auto c = L.LSolve(b);
    return L.Transpose().USolve(c);
}

Matrix Matrix::Cholesky() const
{
    // uses the lower half of the matrix to produce L, s.t L*L^T = original matrix
    try {
        auto A = (*this);
        for (size_t i=0; i<_rows; i++) {
            auto s=0.0;
            for (size_t j=0; j<=i-1; j++) {
                double tm=0;
                for (size_t k=0; k<=j-1; k++) {tm+=A(i,k)*A(j,k);}
                auto t = A(i,j)-tm;
                t /= A(j,j);
                A(i,j) = t;
                s+=t*t;
            }
            s=A(i,i)-s;
            A(i,i)=sqrt(s);
        }
        // erase those above diagonal
        
        for (size_t i=0; i<_rows; i++) {
            for (size_t j=0; j<_rows; j++) {
                if (i<j) {A(i,j)=0;}
            }
        }
        return A;
    }
    catch (...) 
    {
        throw logic_error("Cholesky factorisation failed. Matrix might not be spd.");
    }
}

Matrix Matrix::LSolve(const Matrix& b) const
{
    Matrix x(_rows,1);
    for (size_t i=0; i<_rows; i++)
    {
        auto g=0.0;
        for (size_t j=0; j<i; j++) { g+=(*this)(i,j)*x(j,0); }
        x(i,0) = (b(i,0)-g) / (*this)(i,i);
    }
    return x;
}

Matrix Matrix::USolve(const Matrix &b) const
{
    Matrix x(_rows,1);
    // warning: i = 0, then i-- -> (big int) as i is unsigned long long
    // so we work with index+1 instead.
    for (size_t i=_rows; i>0; i--) 
    {
        auto g=0.0;
        for (size_t j=i-1; j<_rows; j++) { 
            g+=(*this)(i-1,j)*x(j,0); 
        }
        x(i-1,0) = (b(i-1,0)-g) / (*this)(i-1,i-1);
    }
    return x;
}

Matrix Matrix::Solve(const Matrix& b, const double omega) const
{
    Matrix phi(b.Rows(), 1);
    auto tol = 1e-7;
    auto conv = tol * 2;
    while (conv >= tol)
    {
        for (size_t i = 0; i<_rows; i++)
        {
            auto sigma = 0.0;
            for (size_t j = 0; j<_cols; j++)
            {
                if (j != i) {
                    sigma += (*this)(i,j)*phi(j,0);
                }
            }
            phi(i,0) += omega * ((b(i,0)-sigma)/(*this)(i,i) - phi(i,0));
        }

        conv = ((*this)*phi - b).infinityNorm(0);
    }

    return phi;
}

Matrix Matrix::TriSolve(const Matrix& d) const
{
    auto n = d.Rows();
    Matrix x(n, 2);
    auto& m = *this;
    

    x(0,0) = m(0,1) / m(0,0);
    x(0,1) = d(0,0) / m(0,0);

    for (size_t i = 1; i<=n-1; i++)
    {
        auto q = m(i,i) - m(i,i-1) * x(i-1,0);
        auto w = (i == n-1) ? 0 : m(i,i+1) / q;
        auto e = (d(i,0) - x(i-1,1)*m(i,i-1)) / q;
        x(i,0) = w;
        x(i,1) = e;
    }

    x(n-1,0) = x(n-1,1);

    for (size_t i = n-1; i-- > 0;)
    {
        size_t k = x(i,1);
        size_t j = x(i,0);
        auto h = x(i+1,0);
        x(i,0) = k - j * h; // x(i,1)-x(i,0)*x(i+1,0);
    }


    return x.Column(0);
}

void Matrix::initialize(size_t rows, size_t cols)
{
    _M = new double* [rows];
    for (size_t i=0; i<rows; i++) {
        _M[i] = new double [cols];
    }

    _rows = rows;
    _cols = cols;
}

void Matrix::copyHelper(const Matrix& src) 
{
    _rows = src._rows;
    _cols = src._cols;

    _M = new double* [_rows];
    for (size_t i=0; i<_rows; i++) {
        _M[i] = new double [_cols];
    }

    for (size_t i=0; i<_rows; i++) {
        for (size_t j=0; j<_cols; j++) {
            (*this)(i,j) = src(i,j);
        }
    }
}

void Matrix::deleteHelper()
{
    for (size_t i=0; i<_rows; i++) {
        delete[] _M[i];
    }
    delete[] _M;
}

};