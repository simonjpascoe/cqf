/*!
 * Matrix class header
 */
#pragma once
#include <vector>
#include "stdafx.h"

using namespace std;

namespace Common {

class Matrix;
class PCAResult
{
public:
	double Coverage;
	vector<double> Eigenvalues;
	vector<double> Norms;
	vector<Matrix> Eigenvectors;
};

class Matrix 
{
public:
	// construction and destruction
	Matrix();
	Matrix(size_t inSizeR, size_t inSizeC);
	Matrix(const Matrix& src);
	virtual ~Matrix();

	static Matrix Identity(size_t size);
	static Matrix Unit(double v);
	static Matrix FromVector(const vector<double>& source, bool vertical = true);
	vector<double> ToVector(size_t index, bool vertical = true) const;

	// operators
	Matrix& operator=(const Matrix& rhs);
	friend const Matrix operator+(const Matrix& lhs, const Matrix& rhs);
	friend const Matrix operator-(const Matrix& lhs, const Matrix& rhs);
	friend const Matrix operator*(const Matrix& lhs, const Matrix& rhs);
	friend const Matrix operator*(const Matrix& lhs, const double scalar);
	friend const Matrix operator*(const double scalar,const Matrix& rhs);
	friend const Matrix operator/(const Matrix& lhs, const double scalar);

	// accessing and basic information
	double  operator() (size_t row, size_t col) const;
	double& operator() (size_t row, size_t col);
	size_t Rows() const;
	size_t Cols() const;

	// const methods
	Matrix Transpose() const;							// returns matrix transpose
	Matrix Normals() const;								// form M'MX;
	Matrix AugmentR(const Matrix& M) const;				// augments the matrix by M;
	Matrix Row(size_t r) const;							// returns a row matrix
	Matrix Column(size_t c) const;						// returns a column matrix
	Matrix Extract(size_t top_r, size_t top_c, size_t bot_r, size_t bot_c) const;

	double Det() const;			// returns determinant
	double DetSlow() const;		// returns determinant (by first row, very slow)
	double Trace() const;
	double l2Norm(size_t col) const;
	double infinityNorm(size_t col) const;

	// mutable methods
	void Fill(double d);

	// matrix utility methods (also const)
	PCAResult PrincipalComponentAnalysis (double coverageFactor) const;
	Matrix fit(const Matrix& observed) const;
	Matrix Cholesky() const; // for spd matrices. returns lower triangular
	Matrix LSolve(const Matrix& b) const; // solve L*x = b, L = lower triangular
	Matrix USolve(const Matrix& b) const; // solve U*x = b, U = upper triangular
	Matrix Solve(const Matrix& b, const double omega = 1.0) const; 
	Matrix TriSolve(const Matrix& d) const;

private:

	Matrix removeRowCol(size_t row, size_t col) const;
	

	// methods
	void copyHelper(const Matrix& src);
	void deleteHelper(void);

	void initialize(size_t rows, size_t cols);

	// members
	size_t _rows;
	size_t _cols;
	double** _M;		
};

};