#pragma once

#include <cmath>

#include <algorithm>
#include <limits>
#include <ranges>
#include <set>
#include <vector>

#include <Rcpp.h>



using namespace Rcpp;



/*
When the R-facing API changes, regenerate by running (from root dir):
	library("Rcpp");compileAttributes(".");devtools::document()

	R -e 'library("Rcpp");compileAttributes(".");devtools::document()' should work, but doesn't do anything???

Quick testing individual files (from "src/" dir):
	sourceCpp("common.cpp")
	sourceCpp("compute_influence_term_2_quadv_sim_via_matrix.cpp")
	sourceCpp("estimate_pmf.cpp")
	sourceCpp("integrate.cpp")
	sourceCpp("NW.cpp")
	sourceCpp("spline_basis.cpp")

Build package (from root dir):
	R CMD INSTALL .

Set directory:
	setwd("⟨path⟩")
*/



constexpr double RT_TWOPI_RECIP =  0.398942280401432677; // (2π)⁻¹⸍²
constexpr double NEG_HALF_LOG2E = -0.7213475204444817  ; // -½ log₂e

#if __cplusplus >= 202002L
// C++20 (and later) code
#include <format>

template<> struct std::formatter< std::vector<double>, char >
{
	template<class ParseCtx> [[nodiscard]] constexpr
	ParseCtx::iterator parse( ParseCtx& ctx ) { return ctx.begin(); }

	template<class FmtCtx>
	FmtCtx::iterator format( std::vector<double> const& vec, FmtCtx& ctx ) const
	{
		if ( vec.empty() ) return std::format_to( ctx.out(), "{{}}" );
		std::string tmp = "{ ";
		for ( double elem : vec ) tmp+=std::format("{}, ",elem);
		tmp[tmp.length()-2]=' '; tmp[tmp.length()-1]='}';
		return std::format_to( ctx.out(), "{}", tmp );
	}
};

template<class... Args> inline
void print( std::format_string<Args...> fmt, Args&&... args ) noexcept
{
	Rcout << std::format( fmt, std::forward<Args>(args)... );
}

#endif // __cplusplus >= 202002L

template<class T> [[nodiscard]] constexpr
T sq( T val ) noexcept { return val*val; }

[[nodiscard]] constexpr double clone( double val ) noexcept { return val; }

template<class T> [[nodiscard]]
T zero_with_shape_like( T const& shape_example=T() ) noexcept
{
	if      constexpr (std::is_same_v< T, double        >) return 0.0;
	else if constexpr (std::is_same_v< T, NumericVector >)
	{
		NumericVector ret = NumericVector( shape_example.length() );
		ret.fill(0.0);
		return ret;
	}
	else if constexpr (std::is_same_v< T, NumericMatrix >)
	{
		NumericMatrix ret = NumericMatrix(Dimension( shape_example.nrow(), shape_example.ncol() ));
		ret.fill(0.0);
		return ret;
	}
}

void negate_slot( List list, char const* slot_name ) noexcept;



inline NumericMatrix& operator*=( NumericMatrix& matr, double val ) noexcept
{
	for ( double& elem : matr ) elem*=val;
	return matr;
}
inline NumericMatrix& operator/=( NumericMatrix& matr, double val ) noexcept
{
	for ( double& elem : matr ) elem/=val;
	return matr;
}
inline NumericVector& operator*=( NumericVector& vec, double val ) noexcept
{
	for ( double& elem : vec ) elem*=val;
	return vec;
}
inline NumericVector& operator/=( NumericVector& vec, double val ) noexcept
{
	for ( double& elem : vec ) elem/=val;
	return vec;
}
inline NumericMatrix& operator-=( NumericMatrix& matrA, NumericMatrix const& matrB )
{
	if ( matrA.nrow()!=matrB.nrow() || matrA.ncol()!=matrB.ncol() ) // [[unlikely]]
	    stop("Matrix dimension mismatch in subtract (" +
	        std::to_string(matrA.nrow()) + "x" + std::to_string(matrA.ncol()) +
	        " vs. " +
	        std::to_string(matrB.nrow()) + "x" + std::to_string(matrB.ncol()) +
	        ")!");
	for ( int k=0; k<matrA.length(); ++k ) matrA[k]-=matrB[k];
	return matrA;
}
inline NumericMatrix& operator+=( NumericMatrix& matrA, NumericMatrix const& matrB )
{
    if (matrA.nrow() != matrB.nrow() || matrA.ncol() != matrB.ncol()) // [[unlikely]]
        stop("Matrix dimension mismatch in add (" +
            std::to_string(matrA.nrow()) + "x" + std::to_string(matrA.ncol()) +
            " vs. " +
            std::to_string(matrB.nrow()) + "x" + std::to_string(matrB.ncol()) +
            ")!");
    for ( int k=0; k<matrA.length(); ++k ) matrA[k]+=matrB[k];
	return matrA;
}
inline NumericVector& operator-=( NumericVector& vecA, NumericVector const& vecB )
{
	if ( vecA.length() != vecB.length() ) // [[unlikely]]
    	stop("Vector dimension mismatch in subtract (" +
    	    std::to_string(vecA.length()) + " vs. " +
    		std::to_string(vecB.length()) + ")!"
    	    );
	for ( int k=0; k<vecA.length(); ++k ) vecA[k]-=vecB[k];
	return vecA;
}
inline NumericVector& operator+=( NumericVector& vecA, NumericVector const& vecB )
{
	if ( vecA.length() != vecB.length() ) // [[unlikely]]
    	stop("Vector dimension mismatch in add (" +
            std::to_string(vecA.length()) + " vs. " +
            std::to_string(vecB.length()) + ")!"
        	);

	for ( int k=0; k<vecA.length(); ++k ) vecA[k]+=vecB[k];
	return vecA;
}

inline void exp_elems( NumericMatrix* matr ) noexcept
{
	for ( int k=0; k<matr->length(); ++k ) (*matr)[k]=std::exp((*matr)[k]);
}
inline void exp_elems( NumericVector* vec ) noexcept
{
	for ( int k=0; k<vec->length(); ++k ) (*vec)[k]=std::exp((*vec)[k]);
}

template<class T>
[[nodiscard]] inline bool all_finite( T const& val ) noexcept
{
	if constexpr (std::is_same_v< T, double >) return std::isfinite(val);
	else
	{
		for ( int k=0; k<val.length(); ++k ) if ( !std::isfinite(val[k]) ) return false;
		return true;
	}
}

template<class T>
[[nodiscard]] inline double max_diff( T const& val0, T const& val1 ) noexcept
{
	if constexpr (std::is_same_v< T, double >) return val0 - val1;
	else
	{
		double ret = 0.0;
		for ( int k=0; k<val0.length(); ++k ) ret=std::max( ret, val0[k]-val1[k] );
		return ret;
	}
}

template<class T>
[[nodiscard]] inline double max_abs_diff( T const& val0, T const& val1 ) noexcept
{
	if constexpr (std::is_same_v< T, double >) return std::abs( val0 - val1 );
	else
	{
		double ret = 0.0;
		for ( int k=0; k<val0.length(); ++k ) ret=std::max( ret, std::abs(val0[k]-val1[k]) );
		return ret;
	}
}

enum KernelType { Kt_normal, Kt_biweight2, Kt_biweight4 };


[[nodiscard]] constexpr double K_normal   ( double x, double h ) noexcept
{
	//https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/Normal
	x /= h;
	return RT_TWOPI_RECIP * std::exp( -0.5 * sq(x) );
}
[[nodiscard]] constexpr double K_biweight2( double x, double h ) noexcept
{
	if ( std::abs(x) > h ) return 0.0;
	x /= h;
	return (15.0/16.0) * sq(1.0-sq(x));
}
[[nodiscard]] constexpr double K_biweight4( double x, double h ) noexcept
{
	if ( std::abs(x) > h ) return 0.0;
	double y = sq( x / h );
	return (105.0/64.0) * (1.0-3.0*y) * sq(1.0-y);
}



[[nodiscard]] NumericMatrix mmul( NumericMatrix matrA  , NumericMatrix matrB   );
[[nodiscard]] NumericVector mmul( NumericVector row_vec, NumericMatrix matr    );
[[nodiscard]] NumericVector mmul( NumericMatrix matr   , NumericVector col_vec );

//' Multiplies two matrices.  If the first argument is a vector, it is interpreted as a row vector.
//' Otherwise, if the second argument is a vector, it is interpreted as a column vector.
//' @param matrA first matrix
//' @param matrB second matrix
//' @return The product of `matrA` and `matrB`.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
[[nodiscard]] SEXP pcoriaccel_mmul( SEXP matrA, SEXP matrB );

//' Inner product (dot product) of two vectors.
//' @param vecA first vector
//' @param vecB second vector
//' @return scalar product of `vecA` and `vecB`
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
[[nodiscard]] double pcoriaccel_inner_prod( NumericVector vecA, NumericVector vecB );

//' Outer sum of two vectors.
//' @param vecA first vector
//' @param vecB second vector
//' @return Matrix where each element of `vecA`(row) is added to the element of `vecB`(column).
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
[[nodiscard]] NumericMatrix pcoriaccel_outer_sum( NumericVector vecA, NumericVector vecB );

//' Outer product of two vectors.
//' @param vecA first vector
//' @param vecB second vector
//' @return matrix of the outer product of vectors `vecA` and `vecB`.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
[[nodiscard]] NumericMatrix pcoriaccel_outer_prod( NumericVector vecA, NumericVector vecB );

//' Returns the unique elements of a vector, sorted in ascending order.
//' @param vec the vector
//' @return `sort(unique(vec))`
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
[[nodiscard]] NumericVector pcoriaccel_sorted_unique( NumericVector vec );
