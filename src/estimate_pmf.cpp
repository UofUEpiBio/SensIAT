// [[Rcpp::plugins(cpp20)]]

#include "common.h"


template< class Tfloat, KernelType kernel > [[nodiscard]] inline static
    NumericVector _pcoriaccel_estimate_pmf(
            NumericVector const& Xb, NumericVector const& Y,
            Tfloat xi,
            NumericVector const& y_seq,
            Tfloat h
    ) noexcept {


    /*
     Compute kernel applied to each pair of Xbⱼ and xᵢ

     Kxb[j,i] = K( Xb[j]-xb[i], h )

     Equivalent to `outer( Xb,xb, \(a,b){K(a-b,h)} )`
     */

    // pmf_est[j] = sum_k K(Xb[k]-xi,h) * 1(Y[k] == y_seq[j]) / sum_k K(Xb[k]-xi,h)
    NumericVector pmf_est = NumericVector( y_seq.length() );
    Tfloat denom = 0.0;
    for ( int j=0; j<Xb.length(); ++j )
    {
        Tfloat Kxb_j;
        if constexpr( kernel == Kt_normal )
        {
            Kxb_j = K_normal( Xb[j]-xi, h );
        }
        else if constexpr( kernel == Kt_biweight2 )
        {
            Kxb_j = K_biweight2( Xb[j]-xi, h );
        }
        else if constexpr( kernel == Kt_biweight4 )
        {
            Kxb_j = K_biweight4( Xb[j]-xi, h );
        }

        auto i = std::find( y_seq.cbegin(),y_seq.cend(), Y[j] );
        pmf_est[ std::distance(y_seq.cbegin(),i) ] += Kxb_j;

        denom += Kxb_j;
    }

    if ( denom == 0.0 ) // [[unlikely]]
    {
        pmf_est.fill(0.0);
    }
    else
    {
        for ( double& elem : pmf_est ) elem/=denom;
    }

    return pmf_est;
}


//' Directly estimate the probability mass function of Y.
//'
//' @param Xb Numeric vector of individual linear predictors from the data
//' @param Y Numeric vector of individual responses from the data
//' @param xi value of the individuals linear predictor at the point of estimation
//' @param y_seq Numeric vector of unique values of `Y`.
//' @param h bandwidth of the kernel
//' @param kernel character string specifying the kernel to use, either `"dnorm"`, `"K2_Biweight"`, or `"K4_Biweight"`
//'
// [[Rcpp::export]]
[[nodiscard]] NumericVector pcoriaccel_estimate_pmf(
	NumericVector Xb, NumericVector Y,
	double xi,
	NumericVector y_seq,
	double h,
    String kernel = "K2_Biweight"
) {
    if ( Xb.length() != Y.length() ) // [[unlikely]]
    {
        stop("Lengths of arguments `Xb` and `Y` must match!");
    }

    using Tfloat = double;
    if      ( kernel == "dnorm" )
    {
        return _pcoriaccel_estimate_pmf< Tfloat, Kt_normal >( Xb,Y, xi, y_seq, (Tfloat)h );
    }
    else if ( kernel == "K2_Biweight" ) // [[likely]]
    {
        return _pcoriaccel_estimate_pmf< Tfloat, Kt_biweight2 >( Xb,Y, xi, y_seq, (Tfloat)h );
    }
    // else if ( kernel == "K4_Biweight" )
    // {
    //     return _pcoriaccel_estimate_pmf< Tfloat, Kt_biweight4 >( Xb,Y, xi, y_seq, (Tfloat)h );
    // }
    else // [[unlikely]]
    {
        stop("Invalid value for `kernel`: choices are { \"dnorm\", \"K2_Biweight\" }.");
    }
}
