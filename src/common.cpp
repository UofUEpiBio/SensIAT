// [[Rcpp::plugins(cpp20)]]

#include "common.h"
#include <algorithm> // std::sort



void negate_slot( List list, char const* slot_name ) noexcept
{
	auto slot = list.slot(slot_name);
	if      ( Rf_isMatrix(slot) ) list.slot(slot_name)=-as<NumericMatrix>(slot);
	else if ( Rf_isVector(slot) ) list.slot(slot_name)=-as<NumericVector>(slot);
	else                          list.slot(slot_name)=-as<double       >(slot);
}



[[nodiscard]] String pcoriaccel_hello()
{
	return String("Hello from the SensIAT Acceleration C++ sub-library!");
}

[[nodiscard]] NumericMatrix mmul( NumericMatrix matrA  , NumericMatrix matrB   )
{
	if ( matrA.ncol() != matrB.nrow() ) // [[unlikely]]
	{
		stop(
			"Matrix dimension mismatch " +
			    std::to_string(matrA.nrow()) + "x" + std::to_string(matrA.ncol()) + " by " +
			    std::to_string(matrB.nrow()) + "x" + std::to_string(matrB.ncol()) + "!\n"
		);
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
	if ( row_vec.length() != matr.nrow() ) // [[unlikely]]
	{
		stop(
			"Matrix dimension mismatch 1x" + std::to_string(row_vec.length()) + " by " +
			    std::to_string(matr.nrow()) + "x" + std::to_string(matr.ncol()) + "!\n"
		);
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
	if ( matr.ncol() != col_vec.length() ) // [[unlikely]]
	{
		stop(
			"Matrix dimension mismatch " +
			    std::to_string(matr.nrow()) + "x" + std::to_string(matr.ncol()) + " by " +
			    std::to_string(col_vec.length()) + "x1!\n"
		);
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
	if ( vecA.length() != vecB.length() ) // [[unlikely]]
	{
		stop(
			"Matrix dimension mismatch 1x" + std::to_string(vecA.length()) +
			    " by " + std::to_string(vecB.length()) + "x1!\n"
		);
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

#if __cplusplus >= 202002L
	std::ranges::sort(ret);
#else
	std::sort(ret.begin(), ret.end());
#endif
	return ret;
}
