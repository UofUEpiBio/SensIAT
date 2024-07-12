// [[Rcpp::plugins(cpp20)]]

#include "common.h"



[[nodiscard]] NumericVector pcoriaccel_estimate_pmf(
	NumericVector X, NumericVector Y,
	double xi,
	NumericVector y_seq,
	double h
) {
	//auto K = K_normal;
	auto K = K_biweight2;
	//auto K = K_biweight4;

	/*
	Compute kernel applied to each pair of Xbⱼ and xᵢ

	Kxb[j,i] = K( Xb[j]-xb[i], h )

	Equivalent to `outer( Xb,xb, \(a,b){K(a-b,h)} )`
	*/

	// pmf_est[j] = sum_k K(Xb[k]-xi,h) * 1(Y[k] == y_seq[j]) / sum_k K(Xb[k]-xi,h)
	NumericVector pmf_est = NumericVector( y_seq.length() );
	double denom = 0.0;
	for ( int j=0; j<X.length(); ++j )
	{
		double Kxb_j = K( X[j]-xi, h );

		auto i = std::find( y_seq.cbegin(),y_seq.cend(), Y[j] );
		pmf_est[ std::distance(y_seq.cbegin(),i) ] += Kxb_j;

		denom += Kxb_j;
	}

	if ( denom == 0.0 ) [[unlikely]] pmf_est.fill(0.0);
	else
	{
		for ( double& elem : pmf_est ) elem/=denom;
	}

	return pmf_est;
}
