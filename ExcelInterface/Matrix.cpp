#include <xlw/xlw.h>
#include <xlw/XlfServices.h>

#include "ExcelHelpers.h"
#include "Matrix.h"
#include "ObjectCache.h"
#include "ObjectOwner.h"

using namespace xlw;
using namespace std;
using namespace Common;

extern "C" {
	LPXLFOPER EXCEL_EXPORT xlMatrixCreate(LPXLFOPER source)
	{
		EXCEL_BEGIN;
		NO_FUNC_WIZARD;

		XlfOper xSource(source);
		auto src(xSource.AsMatrix("source"));
		Matrix inter(src.size1(), src.size2());
		for (int i = 0; i<src.size1(); i++)
		{
			for (int j =0; j<src.size2(); j++)
			{
				inter(i,j) = src(i,j);
			}
		}

		// ownership management
		auto key = cache_item(inter);

		return XlfOper(key);
		EXCEL_END;
	}

	LPXLFOPER EXCEL_EXPORT xlMatrixRows(LPXLFOPER source)
	{
		EXCEL_BEGIN;
		NO_FUNC_WIZARD;

		XlfOper xSource(source);
		auto inter = ObjectCache<Matrix>::Instance().Lookup(xSource.AsString());
		unsigned long rowcount = inter.Rows();
		return XlfOper(rowcount);
		EXCEL_END;
	}

	LPXLFOPER EXCEL_EXPORT xlMatrixColumns(LPXLFOPER source)
	{
		EXCEL_BEGIN;
		NO_FUNC_WIZARD;

		XlfOper xSource(source);
		auto inter = ObjectCache<Matrix>::Instance().Lookup(xSource.AsString());
		unsigned long colcount = inter.Cols();
		return XlfOper(colcount);
		EXCEL_END;
	}

	LPXLFOPER EXCEL_EXPORT xlMatrixData(LPXLFOPER source)
	{
		EXCEL_BEGIN;
		NO_FUNC_WIZARD;

		XlfOper xSource(source);
		auto inter = ObjectCache<Matrix>::Instance().Lookup(xSource.AsString());
		MyMatrix dest(inter.Rows(), inter.Cols());
		for (int i = 0; i<dest.size1(); i++)
		{
			for (int j =0; j<dest.size2(); j++)
			{
				dest(i,j) = inter(i,j);
			}
		}
		return XlfOper(dest);
		EXCEL_END;
	}

	LPXLFOPER EXCEL_EXPORT xlMatrixRow(LPXLFOPER source, LPXLFOPER row)
	{
		EXCEL_BEGIN;
		NO_FUNC_WIZARD;

		XlfOper xSource(source);
		XlfOper xRow(row);
		auto inter = ObjectCache<Matrix>::Instance().Lookup(xSource.AsString());
		auto dest = inter.Row(xRow.AsULong());

		// ownership management
		auto key = cache_item(dest);

		return XlfOper(key);
		EXCEL_END;
	}

	LPXLFOPER EXCEL_EXPORT xlMatrixColumn(LPXLFOPER source, LPXLFOPER column)
	{
		EXCEL_BEGIN;
		NO_FUNC_WIZARD;

		XlfOper xSource(source);
		XlfOper xColumn(column);
		auto inter = ObjectCache<Matrix>::Instance().Lookup(xSource.AsString());
		auto dest = inter.Column(xColumn.AsULong());

		// ownership management
		auto key = cache_item(dest);

		return XlfOper(key);
		EXCEL_END;
	}

	LPXLFOPER EXCEL_EXPORT xlMatrixValue(LPXLFOPER source, LPXLFOPER row, LPXLFOPER col)
	{
		EXCEL_BEGIN;
		NO_FUNC_WIZARD;

		XlfOper xSource(source);
		auto M = ObjectCache<Matrix>::Instance().Lookup(xSource.AsString());
		XlfOper xr(row);
		XlfOper xc(col);
		return XlfOper(M(xr.AsULong(), xc.AsULong()));

		EXCEL_END;
	}

	LPXLFOPER EXCEL_EXPORT xlMatrixTranspose(LPXLFOPER source)
	{
		EXCEL_BEGIN;
		NO_FUNC_WIZARD;

		XlfOper xSource(source);
		auto M = ObjectCache<Matrix>::Instance().Lookup(xSource.AsString());
		auto T = M.Transpose();
		// ownership management
		auto key = cache_item(T);

		return XlfOper(key);
		EXCEL_END;
	}

	LPXLFOPER EXCEL_EXPORT xlMatrixNormals(LPXLFOPER source)
	{
		EXCEL_BEGIN;
		NO_FUNC_WIZARD;

		XlfOper xSource(source);
		auto M = ObjectCache<Matrix>::Instance().Lookup(xSource.AsString());
		auto N = M.Normals();
		// ownership management
		auto key = cache_item(N);

		return XlfOper(key);
		EXCEL_END;
	}

	LPXLFOPER EXCEL_EXPORT xlMatrixCholesky(LPXLFOPER source)
	{
		EXCEL_BEGIN;
		NO_FUNC_WIZARD;

		XlfOper xSource(source);
		auto M = ObjectCache<Matrix>::Instance().Lookup(xSource.AsString());
		auto C = M.Cholesky();
		// ownership management
		auto key = cache_item(C);

		return XlfOper(key);
		EXCEL_END;
	}

	LPXLFOPER EXCEL_EXPORT xlMatrixExtract(LPXLFOPER source, LPXLFOPER x1, LPXLFOPER y1, LPXLFOPER x2, LPXLFOPER y2)
	{
		EXCEL_BEGIN;
		NO_FUNC_WIZARD;

		XlfOper xSource(source);
		XlfOper xx1(x1);
		XlfOper xx2(x2);
		XlfOper xy1(y1);
		XlfOper xy2(y2);

		auto M = ObjectCache<Matrix>::Instance().Lookup(xSource.AsString());
		auto E = M.Extract(xx1.AsULong(),xy1.AsULong(),xx2.AsULong(),xy2.AsULong());
		// ownership management
		auto key = cache_item(E);

		return XlfOper(key);
		EXCEL_END;
	}

	LPXLFOPER EXCEL_EXPORT xlMatrixMultiply(LPXLFOPER first, LPXLFOPER second)
	{
		EXCEL_BEGIN;
		NO_FUNC_WIZARD;

		XlfOper xSourceL(first);
		XlfOper xSourceR(second);

		auto L = ObjectCache<Matrix>::Instance().Lookup(xSourceL.AsString());
		auto R = ObjectCache<Matrix>::Instance().Lookup(xSourceR.AsString());
		auto M = L*R;
		// ownership management
		auto key = cache_item(M);

		return XlfOper(key);
		EXCEL_END;
	}

	LPXLFOPER EXCEL_EXPORT xlMatrixTSolve(LPXLFOPER A, LPXLFOPER b, LPXLFOPER isLT)
	{
		EXCEL_BEGIN;
		NO_FUNC_WIZARD;

		XlfOper xSource(A);
		XlfOper xVar(b);
		XlfOper op(isLT);

		auto M = ObjectCache<Matrix>::Instance().Lookup(xSource.AsString());
		auto vb = ObjectCache<Matrix>::Instance().Lookup(xVar.AsString());
		auto result = op.AsBool() ? M.LSolve(vb) : M.USolve(vb);

		// ownership management
		auto key = cache_item(result);

		return XlfOper(key);
		EXCEL_END;
	}

	LPXLFOPER EXCEL_EXPORT xlMatrixSolve(LPXLFOPER A, LPXLFOPER b)
	{
		EXCEL_BEGIN;
		NO_FUNC_WIZARD;

		XlfOper xSource(A);
		XlfOper xVar(b);

		auto M = ObjectCache<Matrix>::Instance().Lookup(xSource.AsString());
		auto vb = ObjectCache<Matrix>::Instance().Lookup(xVar.AsString());
		auto result = M.Solve(vb);

		// ownership management
		auto key = cache_item(result);

		return XlfOper(key);
		EXCEL_END;
	}

	LPXLFOPER EXCEL_EXPORT xlMatrixTriSolve(LPXLFOPER A, LPXLFOPER b)
	{
		EXCEL_BEGIN;
		NO_FUNC_WIZARD;

		XlfOper xSource(A);
		XlfOper xVar(b);

		auto M = ObjectCache<Matrix>::Instance().Lookup(xSource.AsString());
		auto vb = ObjectCache<Matrix>::Instance().Lookup(xVar.AsString());
		auto result = M.TriSolve(vb);

		// ownership management
		auto key = cache_item(result);

		return XlfOper(key);
		EXCEL_END;
	}

	/*LPXLFOPER EXCEL_EXPORT xlLinearRegression(LPXLFOPER projx, LPXLFOPER y)
	{
		EXCEL_BEGIN;
		NO_FUNC_WIZARD;

		XlfOper xXs(projx);
		XlfOper xY(y);
		
		auto Xs = ObjectCache<Matrix>::Instance().Lookup(xXs.AsString());
		auto Y  = ObjectCache<Matrix>::Instance().Lookup(xY.AsString());
		auto result = Xs.fit(Y);

		// ownership management
		auto key = cache_item(result);

		return XlfOper(key);
		EXCEL_END;
	}*/
}

namespace {
	// create matrix
	XLRegistration::Arg argsXlCreateMatrix[] = {
        { "source", "Input cells", "XLF_OPER" }
    };

	XLRegistration::XLFunctionRegistrationHelper register__matrix_create(
        "xlMatrixCreate", "matrix_create", "Make a matrix",
        "CQF", argsXlCreateMatrix, 1);

	// one matrix args funcs
	XLRegistration::Arg argsXlOneMatrixFunc[] = {
        { "source", "Handle", "XLF_OPER" }
    };

	XLRegistration::XLFunctionRegistrationHelper register__matrix_data(
        "xlMatrixData", "matrix_data", "Return matrix to cells",
        "CQF", argsXlOneMatrixFunc, 1);

	XLRegistration::XLFunctionRegistrationHelper register__matrix_rows(
        "xlMatrixRows", "matrix_rows", "Return number of rows",
        "CQF", argsXlOneMatrixFunc, 1);

	XLRegistration::XLFunctionRegistrationHelper register__matrix_columns(
        "xlMatrixColumns", "matrix_columns", "Return number of columns",
        "CQF", argsXlOneMatrixFunc, 1);

	XLRegistration::XLFunctionRegistrationHelper register__matrix_transpose(
        "xlMatrixTranspose", "matrix_transpose", "Return matrix tranpose",
        "CQF", argsXlOneMatrixFunc, 1);

	XLRegistration::XLFunctionRegistrationHelper register__matrix_normals(
        "xlMatrixNormals", "matrix_normals", "Return matrix normal M'M",
        "CQF", argsXlOneMatrixFunc, 1);

	XLRegistration::XLFunctionRegistrationHelper register__matrix_cholesky(
        "xlMatrixCholesky", "matrix_cholesky", "Calculate Cholseky decomposition",
        "CQF", argsXlOneMatrixFunc, 1);

	// more complicated funcs
	
	XLRegistration::Arg argsXlMatrixExtract[] = {
		{ "source", "Handle", "XLF_OPER" },
		{ "row1", "Row", "XLF_OPER" },
		{ "col1", "Col", "XLF_OPER" },
		{ "row2", "Row", "XLF_OPER" },
		{ "col2", "Col", "XLF_OPER" }
    };

	XLRegistration::XLFunctionRegistrationHelper register__matrix_extract(
        "xlMatrixExtract", "matrix_extract", "Extract a sub matrix",
        "CQF", argsXlMatrixExtract, 5);

	XLRegistration::Arg argsXlMatrixValue[] = {
		{ "source", "Handle", "XLF_OPER" },
		{ "row", "Row", "XLF_OPER" },
		{ "col", "Col", "XLF_OPER" }
    };

	XLRegistration::XLFunctionRegistrationHelper register__matrix_value(
        "xlMatrixValue", "matrix_value", "Return single matrix value",
        "CQF", argsXlMatrixValue, 3);

	XLRegistration::Arg argsXlMatrixRow[] = {
		{ "source", "Handle", "XLF_OPER" },
		{ "row", "Row", "XLF_OPER" }
	};

	XLRegistration::XLFunctionRegistrationHelper register__matrix_row(
        "xlMatrixRow", "matrix_row", "Return row from matrix",
        "CQF", argsXlMatrixRow, 2);
	
	XLRegistration::Arg argsXlMatrixColumn[] = {
		{ "source", "Handle", "XLF_OPER" },
		{ "column", "Column", "XLF_OPER" }
	};

	XLRegistration::XLFunctionRegistrationHelper register__matrix_column(
        "xlMatrixColumn", "matrix_column", "Return column from matirx",
        "CQF", argsXlMatrixColumn, 2);


	XLRegistration::Arg argsXlMatrixMultiply[] = {
        { "first", "Handle", "XLF_OPER" },
		{ "second", "Handle", "XLF_OPER" },
    };

	XLRegistration::XLFunctionRegistrationHelper register__matrix_multiply(
        "xlMatrixMultiply", "matrix_multiply", "Multiply two matrices",
        "CQF", argsXlMatrixMultiply,2);

	XLRegistration::Arg argsXlMatrixTSolve[] = {
        { "matrix", "Handle", "XLF_OPER" },
		{ "vector", "Handle", "XLF_OPER" },
		{ "isLower", "Bool", "XLF_OPER" }
    };

	XLRegistration::XLFunctionRegistrationHelper register__matrix_tsolve(
        "xlMatrixTSolve", "matrix_tsolve", "Solve lower or upper triangular system",
        "CQF", argsXlMatrixTSolve, 3);

	XLRegistration::Arg argsXlMatrixSolve[] = {
        { "matrix", "Handle", "XLF_OPER" },
		{ "vector", "Handle", "XLF_OPER" }
    };

	XLRegistration::XLFunctionRegistrationHelper register__matrix_solve(
        "xlMatrixSolve", "matrix_solve", "Solve square linear system with SOR",
        "CQF", argsXlMatrixSolve, 2);


	XLRegistration::XLFunctionRegistrationHelper register__matrix_trisolve(
        "xlMatrixTriSolve", "matrix_trisolve", "Solve tri-diagonal system",
        "CQF", argsXlMatrixSolve, 2);

	
/*	XLRegistration::Arg argsXlLinearRegression[] = {
		{ "projx", "Matrix projections x", "XLF_OPER" },
        { "y", "Input matrix", "XLF_OPER" },
	};

	XLRegistration::XLFunctionRegistrationHelper register__linear_regression(
        "xlLinearRegression", "linear_regression", "Perform linear regression",
        "CQF", argsXlLinearRegression, 2);*/
}