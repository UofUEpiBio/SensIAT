#pragma once

#include "common.h"



/*
Basic numerical integration with trapezoidal rule, for testing

∫ f(x) dx ≈ Δx ( ½f(x₁) + f(x₂) + f(x₃) + ⋯ + f(xₙ₋₁) + ½f(xₙ) )
*/
template<class Fn>
[[nodiscard]] inline List integrate_trap(
	Fn&& integrand,
	double lo, double hi,
	unsigned N = 1000
) {
	List ret(3);
	if      ( lo <  hi ) // [[likely]]
	    {} //continue below
	else if ( lo == hi ) return List::create(
		Named("Q"         ) = 0.0,
		Named("fcnt"      ) = 0,
		Named("estim.prec") = 0
	);
	else if ( hi <  lo )
	{
		List tmp = integrate_trap( integrand, lo,hi, N );
		negate_slot( tmp, "Q" );
		return tmp;
	}
	else stop("Invalid integral bounds!"); //NaNs

	double delx = ( hi - lo ) / (double)N;

	auto Q = integrand(lo);
	Q += integrand(hi);
	Q *= 0.5;
	unsigned fcnt = 2;

	for ( unsigned k=1; k<N-1; ++k )
	{
		double x = lo + k*delx; //Step from `lo` each time for accuracy + parallel; 1 FMADD anyway

		Q += integrand(x);
		++fcnt;
	}

	Q *= delx;

	return List::create(
		Named("Q"         ) = Q,
		Named("fcnt"      ) = fcnt,
		Named("estim.prec") = std::numeric_limits<double>::quiet_NaN()
	);
}

/*
Adaptive Simpson integration, ported almost directly from R `pracma::quadv(⋯)`:
	https://rdrr.io/cran/pracma/man/quadv.html
*/
template<class Fn>
[[nodiscard]] inline List integrate_simp(
	Fn&& integrand,
	double lo, double hi,
	double tol = 1.490116e-08 //.Machine$double.eps^(1/2)
) {
	using Tval = decltype(integrand(0.0));

	List ret(3);
	if      ( lo <  hi ) {} //continue below
	else if ( lo == hi ) return List::create(
		Named("Q"         ) = 0.0,
		Named("fcnt"      ) = 0,
		Named("estim.prec") = 0
	);
	else if ( hi <  lo )
	{
		List tmp = integrate_simp( integrand, lo,hi, tol );
		negate_slot( tmp, "Q" );
		return tmp;
	}
	else stop("Invalid integral bounds!"); //NaNs

	//Smallest ε such that 1+ε != 1   ( ε ≈ 2.220446⨯10⁻¹⁶ )
	constexpr double eps = std::numeric_limits<double>::epsilon();
	double const hmin = (hi-lo) * (eps/1024.0);

	//Initial values separated by unequal intervals
	//	Value apparently chosen so as not to alias with many integrands; see:
	//	https://www.mathworks.com/matlabcentral/answers/1894-why-0-13579
	double h = 0.13579 * (hi-lo);
	double xs[7] = { lo, lo+h, lo+2*h, 0.5*(lo+hi), hi-2*h, hi-h, hi };
	Tval   ys[7] = {
		integrand(xs[0]), integrand(xs[1]), integrand(xs[2]),
		integrand(xs[3]),
		integrand(xs[4]), integrand(xs[5]), integrand(xs[6])
	};
	Tval Qret = zero_with_shape_like(ys[0]);
	unsigned fcnt = 7;

	//Fudge endpoints to avoid infinities
	if ( !all_finite(ys[0]) )
	{
		ys[0] = integrand( lo + eps*(hi-lo) );
		++fcnt;
	}
	if ( !all_finite(ys[6]) )
	{
		ys[6] = integrand( hi - eps*(hi-lo) );
		++fcnt;
	}

	//Recursive calls to main integrator function
	auto helper = [ integrand, tol,hmin, &Qret,&fcnt ](
		auto& helper_ref, double a,double c,double e, Tval const& fa,Tval const& fc,Tval const& fe
	) {
		constexpr unsigned fcnt_max = 10000;
		if ( fcnt+2 > fcnt_max ) // [[unlikely]]
		    stop("Too many integrand evaluations; singularity likely.");
		//if ( depth > 5 ) stop("abort");

		//Calculate subintervals
		double h = e - a;
		if ( h<hmin || a==c || c==e ) // [[unlikely]]
		    stop("Minimum step size; singularity possible.");
		double b = 0.5*( a + c );
		double d = 0.5*( c + e );

		//Evaluate integrals at new midpoints
		Tval fb = integrand(b);
		Tval fd = integrand(d);
		fcnt += 2;

		//Simpson's 1/3 rule applied to points (a,c,e):   (h/6)( f(a) + 4f(c) + f(e) )
		Tval Q1=clone(fc); Q1*=4.0; Q1+=fa; Q1+=fe; Q1*=h*(1.0/6.0);
		//Same, for points (a,b,c) added to (c,d,e):   (h/12)( f(a) + 4f(b) + 2f(c) + 4f(d) + f(e) )
		Tval Q2=clone(fb); Q2+=fd; Q2*=2.0; Q2+=fc; Q2*=2.0; Q2+=fa; Q2+=fe; Q2*=h*(1.0/12.0);
		//One step of Romberg extrapolation:   Q = Q2 + (Q2-Q1)/15
		Tval Q=clone(Q2); Q-=Q1; Q*=1.0/15.0; Q+=Q2;
		if ( !all_finite(Q) ) // [[unlikely]]
		    stop("Improper integrand: ∞ or NaN");

		//If accurate enough, return.  Original was just max(Q2-Q)<tol; doesn't sound correct to me.
		//	TODO: can optimize with prev.
		//double diff = max_diff( Q2, Q );
		double diff = max_abs_diff( Q2, Q );
		//for (unsigned i=0;i<depth;++i) print("  ");
		//print("a={:.4f}, b={:.4f}, c={:.4f}, d={:.4f}, e={:.4f} (diff={:e}/{:e}, depth={})\n",a,b,c,d,e,diff,tol,depth);
		//Rcout << "Q1=\n"<<Q1 << "Q2=\n"<<Q2 << "Q=\n"<<Q << "\n\n";
		//Rcout << "f(a)=\n"<<fa << "f(b)=\n"<<fb << "f(c)=\n"<<fc << "f(d)=\n"<<fd << "f(e)=\n"<<fe << "\n";
		if ( diff < tol )
		{
			Qret += Q;
			return;
		}
		//Otherwise, recurse
		helper_ref( helper_ref, a,b,c, fa,fb,fc );
		helper_ref( helper_ref, c,d,e, fc,fd,fe );
	};
	helper( helper, xs[0],xs[1],xs[2], ys[0],ys[1],ys[2] );
	helper( helper, xs[2],xs[3],xs[4], ys[2],ys[3],ys[4] );
	helper( helper, xs[4],xs[5],xs[6], ys[4],ys[5],ys[6] );

	return List::create(
		Named("Q"         ) = Qret,
		Named("fcnt"      ) = fcnt,
		Named("estim.prec") = tol * 0.5 * (fcnt-7)
	);
}



//' Integrate function using adaptive Simpson quadrature.
//'
//' @param integrand   The integrand, must take scalar argument, may return scalar, vector, or matrix.
//' @param lo          Lower integration bound
//' @param hi          Upper integration bound
//' @param tol         Tolerance for integration, default `.Machine$double.eps^(1/2)`
//'
//' @return integration result, list with elements `$Q` (the integral estimate), `$fcnt` (the number
//' of function evaluations), and `$estim.prec` (a (pessimistic) estimate of the precision).
//'
//' @keywords internal
// [[Rcpp::export]]
[[nodiscard]] List pcoriaccel_integrate_simp(
	Function integrand,
	double lo, double hi,
	double tol = 1.490116e-08 //.Machine$double.eps^(1/2)
) {
	auto test = integrand(lo);
	List ret;
	if      ( Rf_isMatrix(test) ) ret=integrate_simp(
		[integrand]( double x ) { return as<NumericMatrix>(integrand(x)); }, lo,hi, tol
	);
	else if ( Rf_isVector(test) ) ret=integrate_simp(
		[integrand]( double x ) { return as<NumericVector>(integrand(x)); }, lo,hi, tol
	);
	else                          ret=integrate_simp(
		[integrand]( double x ) { return as<double>(integrand(x)); }, lo,hi, tol
	);
	ret["fcnt"] = as<int>(ret["fcnt"]) + 1;
	return ret;
}
