#pragma once

#include "common.h"




[[nodiscard]] NumericVector pcoriaccel_estimate_pmf(
	NumericVector X, NumericVector Y,
	double xi,
	NumericVector y_seq,
	double h,
	String kernel = "K2_Biweight"
);
