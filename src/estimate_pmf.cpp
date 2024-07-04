// [[Rcpp::plugins(cpp20)]]

#include "common.hpp"



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
	NumericVector Kxb = NumericVector(X.length());
	for ( int j=0; j<Kxb.length(); ++j )
	{
		Kxb(j) = K( X[j]-xi, h );
	}

	// est_pmf[j] = sum_k K(Xb[k]-xi,h) * 1(Y[k] == y_seq[j]) / sum_k K(Xb[k]-xi,h)
	NumericVector est_pmf = NumericVector( y_seq.length() );
	double denom = 0.0;
	for ( int j=0; j<Kxb.length(); ++j )
	{
		denom += Kxb(j);

		auto i = std::find( y_seq.begin(), y_seq.end(), Y[j] );

		est_pmf(std::distance(y_seq.begin(),i)) += Kxb(j);
	}

	if ( denom == 0.0 ) [[unlikely]]
	{
		est_pmf.fill(0.0);
	}
	else
	{
		for ( double& elem : est_pmf ) elem/=denom;
	}

	return est_pmf;
}
