// [[Rcpp::plugins(cpp20)]]

//#include "common.hpp"
#include "estimate_pmf.hpp"
#include "integrate.hpp"
#include "spline_basis.hpp"



//' Runs an optimized implementation of the `compute_influence_term_2_quadv_sim_via_matrix`
//' function.
//'
//' @param X              Matrix of all covariates, transformed as necessary by model
//' @param Y              Vector of all outcomes (same length as a column of `X`)
//' @param times          Vector of observation times for individual
//' @param individual_X   Matrix of covariates for individual rows correspond to times prepared for
//'                       inferences for integration.
//' @param x_slope        Vector of numeric(length(beta)) indicating how
//' @param alpha          Vector of sensitivity parameters
//' @param beta           Vector of coefficients of the outcome model
//' @param spline_basis   Spline basis object (`orthogonalsplinebasis::SplineBasis`)
//' @param bandwidth      Bandwidth for the kernel density estimate of the outcome model.
//' @param tol            Tolerance for integration
//' @param ...            Additional arguments passed to the pcoriaccel_NW function, not implemented.
//'
//' @return integration result
//'
//' @export
// [[Rcpp::export]]
[[nodiscard]] NumericMatrix pcoriaccel_compute_influence_term_2_quadv_sim_via_matrix(
	NumericMatrix X,            // e.g. num[1:453,1:5]
	NumericVector Y,            // e.g. num[1:453]

	NumericVector times,        // e.g. num[1:3]
	NumericMatrix individual_X, // e.g. num[1:3,1:5]
	NumericVector x_slope,      // e.g. num[1:5] numeric(length(beta))

	NumericVector alpha,        // e.g. num[1:3]
	NumericVector beta ,        // e.g. num[1:5]

	S4 spline_basis,

	double bandwidth,
	double tol = 0.0001220703 // .Machine$double.eps^(1/4)
) {
	#if 1
	if ( X.ncol()!=beta.length() || X.ncol()!=individual_X.ncol() ) [[unlikely]] stop(std::format(
		"Width of matrix `X` ({}), width of matrix `individual_X` ({}), "
		"and length of vector `beta` ({}) must match!",
		X.ncol(), individual_X.ncol(), beta.length()
	));
	if ( X.nrow() != Y.length() ) [[unlikely]] stop(std::format(
		"Height of matrix `X` ({}) and length of vector `Y` ({}) must match!",
		X.nrow(), Y.length()
	));
	if ( individual_X.ncol() != x_slope.length() ) [[unlikely]] stop(std::format(
		"Width of matrix `individual_X` ({}) and length of vector `x_slope` ({}) must match!",
		individual_X.ncol(), x_slope.length()
	));
	if ( individual_X.nrow() != times.length() ) [[unlikely]] stop(std::format(
		"Height of matrix `individual_X` ({}) and length of vector `times` ({}) must match!",
		individual_X.nrow(), times.length()
	));
	#endif

	//Spline basis wrapper
	SplineBasis basis(spline_basis);

	//Global integral bounds [a,b]
	double a = basis.get_lo_knot();
	double b = basis.get_hi_knot();

	//Time points determining sub-integrals: basically the list `times`, clipped to the endpoint
	//	knots, and ensuring we have the endpoint knots themselves.
	std::vector<double> period_times = { a };
	auto push_period_times = [ a,b, &period_times ]( double val )
	{
		if ( val<a || val>b || val==period_times.back() ) return;
		period_times.emplace_back(val);
	};
	for ( double elem : times ) push_period_times(elem);
	push_period_times(b);
	//print( "a={},b={}, period_times={}\n", a,b, period_times );

	//Compute globally used terms
	NumericVector distinct_Y = pcoriaccel_sorted_unique(Y); // e.g. num[1:35]
	NumericVector Xb = mmul( X, beta );                     // e.g. num[1:453,1]

	//Compute integrals for each period
	std::vector< IntegrateResult<NumericMatrix> > period_integrals( period_times.size() - 1 );
	for ( int period_ind=0; period_ind<(int)period_integrals.size(); ++period_ind )
	{
		double lower = period_times[ period_ind     ];
		double upper = period_times[ period_ind + 1 ];

		auto fn = [&]( double t )
		{
			/*
			The coefficient vector for the given individual is determined by extending the last
			observation to the current time point.  Time-dependent variables scale(t) and
			scale(delta_time) are determined by the centering parameters, which are global
			constants.  All other variables should be last observed values.  All this is encoded in
			the `x_slope` object.
			*/
			NumericVector x_at_time( x_slope.length() );
			double x_slope_sc = t - times[period_ind];
			for ( int k=0; k<x_slope.length(); ++k )
			{
				x_at_time[k] = individual_X(period_ind,k) + x_slope[k]*x_slope_sc;
			}

			//Compute the pmf for `Y` at the given time point (e.g. num[1:35]).
			NumericVector pmf = pcoriaccel_estimate_pmf(
				Xb,Y, pcoriaccel_inner(x_at_time,beta), distinct_Y, bandwidth
			);

			//Compute the expected value of the outcome at the given time point, given the
			//	sensitivity parameter `alpha` (e.g. num[1:35,1:3]).
			NumericMatrix eay = pcoriaccel_outer( distinct_Y, alpha );
			for ( double& elem : eay ) elem=std::exp(elem);
			//	numerator   = transpose( distinct_Y.elems * eay.elems ) * pmf
			//	denominator = transpose(                    eay       ) * pmf
			NumericVector ev( eay.ncol() );
			for ( int i=0; i<eay.ncol(); ++i )
			{
				double numer=0.0, denom=0.0;
				for ( int k=0; k<pmf.length(); ++k )
				{
					double denom_term = eay(k,i) * pmf[k];          //    exp(αyₖ) pmfₖ
					double numer_term = distinct_Y[k] * denom_term; // yₖ exp(αyₖ) pmfₖ
					numer += numer_term;
					denom += denom_term;
				}
				ev[i] = numer / denom;
			}

			//Evaluate the spline basis functions at the given time point.
			NumericVector B = basis.evaluate(t);

			//Returned value is a length(alpha) by ncol(B) matrix.
			NumericMatrix ret = pcoriaccel_outer( ev, B );
			//Rcout << ret << "\n";
			return ret;
		};
		//Rcout << "Integrand value   ind , time   =   " << period_ind << " , " << lower << ":\n";
		//fn(lower);

		auto integrated = integrate_trap( fn, lower,upper );
		//Rcout << "Integral value   ind , time   =   " << period_ind << " ,( " << lower << " , " << upper << "):\n" << integrated.Q;
		period_integrals[period_ind] = integrated;
	}

	//Sum integrals of periods to get entire integral (use first one as accumulator).
	//	Note by above construction we do have at least two points (one integral period).
	for ( size_t k=1; k<period_integrals.size(); ++k )
	{
		period_integrals[0] += period_integrals[k];
	}

	return period_integrals[0].Q;
}
