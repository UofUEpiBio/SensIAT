#pragma once

#include "common.hpp"



//' Estimate the PMF directly with the K2_Biweight kernel.
//'
//' @param Xb vector (expected to be about 500 elements)
//' @param Y vector (same size as Xb)
//' @param xi vector
//' @param y_seq vector
//' @param h scalar bandwidth of kernel
//'
//' @return estimated PMF
//'
// [[Rcpp::export]]
[[nodiscard]] NumericVector pcoriaccel_estimate_pmf(
	NumericVector X, NumericVector Y,
	double xi,
	NumericVector y_seq,
	double h
);
