// [[Rcpp::plugins(cpp20)]]

#include "spline_basis.h"



SplineBasis::SplineBasis( S4 backing )
{
	//Unpack the object slots once, here, so we don't have to do that for each access

	if ( !backing.hasSlot("knots"   ) ) // [[unlikely]]
    	stop(
    		"Spline basis expected to have slot `knots`!"
    	);
	_knots = (NumericVector)backing.slot("knots"   );

	if ( !backing.hasSlot("order"   ) ) // [[unlikely]]
    	stop(
    		"Spline basis expected to have slot `order`!"
    	);
	_order = (int          )backing.slot("order"   );

	if ( !backing.hasSlot("Matrices") ) // [[unlikely]]
    	stop(
    		"Spline basis expected to have slot `Matrices`!"
    	);
	_matrs = (NumericVector)backing.slot("Matrices");

	if ( !_matrs.hasAttribute("dim")  ) // [[unlikely]]
    	stop(
    		"Spline basis `Matrix` expected to have attribute `dim`!"
    	);
	_mdims = (IntegerVector)_matrs.attr("dim");

	//!backing.hasSlot("class") -> "SplineBasis" (apparently optional)
}

[[nodiscard]] NumericVector SplineBasis::evaluate( double x ) const noexcept
{
	double lo = get_lo_knot();
	double hi = get_hi_knot();

	//print( "Evaluate x={} in knots [lo={},hi={}], order {}\n", x, lo,hi, _order );

	NumericVector ret( _mdims[1] );

	//Out-of-range `x`; note test written so as to also catch NaN
	if (!( x>=lo && x<=hi )) // [[unlikely]]
	{
		//print("  -> out of range\n");
		ret.fill(NA_REAL);
		return ret;
	}

	//Find knots `x` is between; binary search not worth overhead at this size.  Loop test
	//	needed only to handle `x == hi` edge case.
	int iknot0;
	for ( iknot0=_order-1; iknot0+1<_knots.length()-_order; ++iknot0 )
	{
		if ( _knots[iknot0+1] > x ) break;
	}

	//Evaluate
	double u = ( x - _knots[iknot0] ) / ( _knots[iknot0+1] - _knots[iknot0] );
	NumericVector U(_order);
	for ( int k=0; k<_order; ++k ) U[k]=std::pow(u,k);
	//print( "  knots {} to {} with lerp {}\n", iknot0,iknot0+1, u );
	//Rcout << "  U=" << U << "\n";
	//Rcout << "  _matrs=" << _matrs << "\n";

	//TODO: use our `mmul(⋯)`?
	//	U * _matrs[i];
	for ( int imatrcol=0; imatrcol<_mdims[1]; ++imatrcol )
	{
		double dp = 0.0;
		for ( int k=0; k<_mdims[0]; ++k )
		{
			int ind = _mdims[0]*_mdims[1]*(iknot0-(_order-1)) + _mdims[0]*imatrcol + k;
			dp += U[k] * _matrs[ind];
		}
		ret[imatrcol] = dp;
	}

	return ret;
}



[[nodiscard]] NumericVector pcoriaccel_evaluate_basis(
	S4 spline_basis, double x
) {
	SplineBasis basis(spline_basis);
	return basis.evaluate(x);
}

[[nodiscard]] NumericMatrix pcoriaccel_evaluate_basis_mat(
        S4 spline_basis, NumericVector x
) {
    SplineBasis basis(spline_basis);

    NumericMatrix ret(x.length(), basis.ncol());

    for( int i=0; i<x.length(); ++i )
    {
        ret(i,_) = basis.evaluate(x[i]);
    }

    return ret;
}

