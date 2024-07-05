#pragma once

#include <cmath>

#include <algorithm>
#include <format>
#include <ranges>
#include <set>
#include <vector>

#include <Rcpp.h>



using namespace Rcpp;



/*
When the R-facing API changes, regenerate by running (from root dir):
	library("Rcpp")
	compileAttributes("../")
	devtools::document()

	R -e 'library("Rcpp");compileAttributes("../");devtools::document()' should work, but doesn't do anything???

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

template<class T> [[nodiscard]] constexpr
T sq( T val ) noexcept { return val*val; }



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
inline NumericMatrix& operator-=( NumericMatrix& matrA, NumericMatrix const& matrB )
{
	if ( matrA.nrow()!=matrB.nrow() || matrA.ncol()!=matrB.ncol() ) [[unlikely]] stop("Matrix dimension mismatch!");
	for ( int k=0; k<matrA.length(); ++k ) matrA[k]-=matrB[k];
	return matrA;
}
inline NumericMatrix& operator+=( NumericMatrix& matrA, NumericMatrix const& matrB )
{
	if ( matrA.nrow()!=matrB.nrow() || matrA.ncol()!=matrB.ncol() ) [[unlikely]] stop("Matrix dimension mismatch!");
	for ( int k=0; k<matrA.length(); ++k ) matrA[k]+=matrB[k];
	return matrA;
}

inline void exp_elems( NumericMatrix* matr ) noexcept
{
	for ( int k=0; k<matr->length(); ++k ) (*matr)[k]=std::exp((*matr)[k]);
}
inline void exp_elems( NumericVector* vec ) noexcept
{
	for ( int k=0; k<vec->length(); ++k ) (*vec)[k]=std::exp((*vec)[k]);
}



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



//' Returns a string, just as a basic check that the C++ plugin library is working.
//' @return hello string
//' @export
// [[Rcpp::export]]
[[nodiscard]] String pcoriaccel_hello();

[[nodiscard]] NumericMatrix mmul( NumericMatrix matrA  , NumericMatrix matrB   );
[[nodiscard]] NumericVector mmul( NumericVector row_vec, NumericMatrix matr    );
[[nodiscard]] NumericVector mmul( NumericMatrix matr   , NumericVector col_vec );
//' Multiplies two matrices.  If the first argument is a vector, it is interpreted as a row vector.
//' Otherwise, if the second argument is a vector, it is interpreted as a column vector.
//' @param matrA first matrix
//' @param matrB second matrix
//' @return matrA * matrB
//' @export
// [[Rcpp::export]]
[[nodiscard]] SEXP pcoriaccel_mmul( SEXP matrA, SEXP matrB );

//' Inner product (dot product) of two vectors.
//' @param vecA first vector
//' @param vecB second vector
//' @return vecAᵀ * vecB = vecA • vecB
//' @export
// [[Rcpp::export]]
[[nodiscard]] double pcoriaccel_inner_prod( NumericVector vecA, NumericVector vecB );
//' Outer sum of two vectors.
//' @param vecA first vector
//' @param vecB second vector
//' @return vecA ⊕ vecB
//' @examples
//' pcoriaccel_outer_sum( c(1,2,3,4,5), c(2,4,6) )
//' @export
// [[Rcpp::export]]
[[nodiscard]] NumericMatrix pcoriaccel_outer_sum( NumericVector vecA, NumericVector vecB );
//' Outer product of two vectors.
//' @param vecA first vector
//' @param vecB second vector
//' @return vecA * vecBᵀ = vecA ⊗ vecB
//' @examples
//' pcoriaccel_outer_prod( c(1,2,3,4,5), c(2,4,6) )
//' @export
// [[Rcpp::export]]
[[nodiscard]] NumericMatrix pcoriaccel_outer_prod( NumericVector vecA, NumericVector vecB );

//' Returns the unique elements of a vector, sorted in ascending order.
//' @param vec the vector
//' @return sort(unique(vec))
//' @export
// [[Rcpp::export]]
[[nodiscard]] NumericVector pcoriaccel_sorted_unique( NumericVector vec );
