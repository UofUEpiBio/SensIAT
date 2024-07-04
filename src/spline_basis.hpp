#pragma once

#include "common.hpp"



//Wraps S4 class `orthogonalsplinebasis::SplineBasis`
class SplineBasis final
{
	private:
		NumericVector _knots; // e.g. double[9]
		int           _order; // integer
		NumericVector _matrs; // e.g. double[4x5x2]
		IntegerVector _mdims; // e.g. {4,5,2}

	public:
		explicit SplineBasis( S4 backing );

		[[nodiscard]] inline double get_lo_knot() const noexcept
		{
			return _knots[ _order - 1 ];
		}
		[[nodiscard]] inline double get_hi_knot() const noexcept
		{
			return _knots[ _knots.length() - _order ];
		}

		[[nodiscard]] NumericVector evaluate( double x ) const noexcept;
};

//' Rcpp version of `evaluate_basis(⋯)` function
//'
//' @param spline_basis   Spline basis, S4 class `orthogonalsplinebasis::SplineBasis`
//' @param x              Evaluation point
//'
//' @return Vector fyxb
//'
// [[Rcpp::export]]
[[nodiscard]] NumericVector pcoriaccel_evaluate_basis(
	S4 spline_basis, double x
);
