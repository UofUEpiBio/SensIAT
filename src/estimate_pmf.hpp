#pragma once

#include "common.hpp"



//' Estimate the PMF directly with the K2_Biweight kernel.
//'
//' @param Xb    a vector (expected to be about 500 elements)
//' @param Y     a vector (same size as Xb)
//' @param xi    a scalar
//' @param y_seq a vector
//' @param h     a scalar, the bandwidth of kernel
//'
//' @return estimated PMF
//'
//' @export
// [[Rcpp::export]]
[[nodiscard]] NumericVector pcoriaccel_estimate_pmf(
	NumericVector X, NumericVector Y,
	double xi,
	NumericVector y_seq,
	double h
);
