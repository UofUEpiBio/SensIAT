// [[Rcpp::plugins(cpp20)]]

#include "common.hpp"



[[nodiscard]] String pcoriaccel_hello()
{
	return String("Hello from the PCORI Acceleration C++ sub-library!");
}

[[nodiscard]] NumericMatrix mmul( NumericMatrix matrA  , NumericMatrix matrB   )
{
	if ( matrA.ncol() != matrB.nrow() ) [[unlikely]]
	{
		stop(std::format(
			"Matrix dimension mismatch {}x{} by {}x{}!\n",
			matrA.nrow(),matrA.ncol(), matrB.nrow(),matrB.ncol()
		));
	}

	NumericMatrix ret(Dimension( matrA.nrow(), matrB.ncol() ));
	for ( int irow=0; irow<matrA.nrow(); ++irow )
	for ( int icol=0; icol<matrB.ncol(); ++icol )
	{
		double dp = 0.0;
		for ( int k=0; k<matrA.ncol(); ++k )
		{
			dp += matrA(irow,k) * matrB(k,icol);
		}
		ret(irow,icol) = dp;
	}
	return ret;
}
[[nodiscard]] NumericVector mmul( NumericVector row_vec, NumericMatrix matr    )
{
	if ( row_vec.length() != matr.nrow() ) [[unlikely]]
	{
		stop(std::format(
			"Matrix dimension mismatch {}x{} by {}x{}!\n",
			1,row_vec.length(), matr.nrow(),matr.ncol()
		));
	}

	NumericVector ret( matr.ncol() );
	for ( int icol=0; icol<matr.ncol(); ++icol )
	{
		double dp = 0.0;
		for ( int k=0; k<row_vec.length(); ++k )
		{
			dp += row_vec[k] * matr(k,icol);
		}
		ret[icol] = dp;
	}
	return ret;
}
[[nodiscard]] NumericVector mmul( NumericMatrix matr   , NumericVector col_vec )
{
	if ( matr.ncol() != col_vec.length() ) [[unlikely]]
	{
		stop(std::format(
			"Matrix dimension mismatch {}x{} by {}x{}!\n",
			matr.nrow(),matr.ncol(), col_vec.length(),1
		));
	}

	NumericVector ret( matr.nrow() );
	for ( int irow=0; irow<matr.nrow(); ++irow )
	{
		double dp = 0.0;
		for ( int k=0; k<matr.ncol(); ++k )
		{
			dp += matr(irow,k) * col_vec[k];
		}
		ret[irow] = dp;
	}
	return ret;
}
[[nodiscard]] SEXP pcoriaccel_mmul( SEXP matrA, SEXP matrB )
{
	if      ( Rf_isVector(matrA) && Rf_isMatrix(matrB) )
	{
		return mmul( as<NumericVector>(matrA), as<NumericMatrix>(matrB) );
	}
	else if ( Rf_isMatrix(matrA) && Rf_isVector(matrB) )
	{
		return mmul( as<NumericMatrix>(matrA), as<NumericVector>(matrB) );
	}
	else if ( Rf_isMatrix(matrA) && Rf_isMatrix(matrB) )
	{
		return mmul( as<NumericMatrix>(matrA), as<NumericMatrix>(matrB) );
	}
	else
	{
		stop("Unknown types for matrix multiplication.");
	}
}

[[nodiscard]] double pcoriaccel_inner_prod( NumericVector vecA, NumericVector vecB )
{
	if ( vecA.length() != vecB.length() ) [[unlikely]]
	{
		stop(std::format(
			"Matrix dimension mismatch {}x{} by {}x{}!\n",
			1,vecA.length(), vecB.length(),1
		));
	}

	double dp = 0.0;
	for ( int k=0; k<vecA.length(); ++k )
	{
		dp += vecA[k] * vecB[k];
	}
	return dp;
}
[[nodiscard]] NumericMatrix pcoriaccel_outer_sum( NumericVector vecA, NumericVector vecB )
{
	NumericMatrix ret(Dimension( vecA.length(), vecB.length() ));
	for ( int irow=0; irow<vecA.length(); ++irow )
	for ( int icol=0; icol<vecB.length(); ++icol )
	{
		ret(irow,icol) = vecA[irow] + vecB[icol];
	}
	return ret;
}
[[nodiscard]] NumericMatrix pcoriaccel_outer_prod( NumericVector vecA, NumericVector vecB )
{
	NumericMatrix ret(Dimension( vecA.length(), vecB.length() ));
	for ( int irow=0; irow<vecA.length(); ++irow )
	for ( int icol=0; icol<vecB.length(); ++icol )
	{
		ret(irow,icol) = vecA[irow] * vecB[icol];
	}
	return ret;
}

[[nodiscard]] NumericVector pcoriaccel_sorted_unique( NumericVector vec )
{
	std::set<double> tmp;
	for ( double elem : vec ) tmp.emplace(elem);

	NumericVector ret( tmp.size() );
	int k = 0;
	for ( double val : tmp ) ret[k++]=val;

	std::ranges::sort(ret);

	return ret;
}
