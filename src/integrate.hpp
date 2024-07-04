#pragma once

#include "common.hpp"



template<class Tval>
struct IntegrateResult final
{
	Tval     Q;          //The estimated value
	unsigned fcnt;       //The number of function calls
	//Tval     estim_prec; //Pessimistic precision

	IntegrateResult& operator+=( IntegrateResult const& other ) noexcept
	{
		Q    += other.Q;
		fcnt += other.fcnt;
		return *this;
	}
};

/*
Basic numerical integration with trapezoidal rule, for testing

∫ f(x) dx ≈ Δx ( ½f(x₁) + f(x₂) + f(x₃) + ⋯ + f(xₙ₋₁) + ½f(xₙ) )
*/
template<class Fn>
[[nodiscard]] inline auto integrate_trap(
	Fn&& integrand,
	double lo, double hi,
	unsigned N = 1000
) noexcept {
	using Tval = decltype(integrand(0.0));

	double delx = ( hi - lo ) / (double)N;

	IntegrateResult<Tval> ret;
	ret.Q =  integrand(lo);
	ret.Q += integrand(hi);
	ret.Q *= 0.5;
	ret.fcnt = 2;

	for ( unsigned k=1; k<N-1; ++k )
	{
		double x = lo + k*delx; //Step from `lo` each time for accuracy + parallel; 1 FMADD anyway

		ret.Q += integrand(x);
		++ret.fcnt;
	}

	ret.Q *= delx;

	return ret;
}
