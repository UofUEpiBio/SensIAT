#include <cmath>
#include <vector>

#include <Rcpp.h>



using namespace Rcpp;



/*
When the R-facing API changes, regenerate by running, from R in the project root directory:

	library("Rcpp")
	compileAttributes(".")

	R -e 'library("Rcpp"); compileAttributes(".")' #should work but doesn't?
*/



constexpr double RT_TWOPI_RECIP = 0.398942280401432677;

template<class T> [[nodiscard]] constexpr
T sq( T val ) noexcept { return val*val; }



//' Returns a string, just as a basic check that the C++ plugin library is working.
//'
//' @return hello string
// [[Rcpp::export]]
[[nodiscard]] String pcoriaccel_hello() noexcept
{
	return String("Hello from the PCORI Acceleration C++ sub-library!");
}



[[nodiscard]] constexpr double K_normal( double x, double h ) noexcept
{
	//https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/Normal
	x /= h;
	return RT_TWOPI_RECIP * std::exp( -0.5 * sq(x) );
}
[[nodiscard]] constexpr double K_biweight2( double x, double h ) noexcept
{
	if ( std::abs(x) > h ) return 0.0;
	x /= h;
	return (15.0/16.0) * sq(1.0-sq(x));
}
[[nodiscard]] constexpr double K_biweight4( double x, double h ) noexcept
{
	if ( std::abs(x) > h ) return 0.0;
	double y = sq( x / h );
	return (105.0/64.0) * (1.0-3.0*y) * sq(1.0-y);
}

//' Runs a *basic* implementation of the "NW" function with the "K2_Biweight" kernel, just as a
//' proof-of-concept.
//'
//' @param Xb vector
//' @param Y vector (same size as Xb)
//' @param xb vector
//' @param y_seq vector
//' @param h scalar bandwidth of kernel
//' @return Matrix fyxb
// [[Rcpp::export]]
[[nodiscard]] NumericMatrix pcoriaccel_NW_basic(
	NumericVector Xb, NumericVector Y, //Same length, about 500
	NumericVector xb,
	NumericVector y_seq, //unique values
	double h
) noexcept {
	//auto K = K_normal;
	auto K = K_biweight2;
	//auto K = K_biweight4;

	/*
	Compute kernel applied to each pair of Xbⱼ and xᵢ

	Kxb[j,i] = K( Xb[j]-xb[i], h )

	Equivalent to `outer( Xb,xb, \(a,b){K(a-b,h)} )`
	*/
	NumericMatrix Kxb = NumericMatrix(Dimension( Xb.length(), xb.length() ));
	for ( int j=0; j<Kxb.nrow(); ++j )
	for ( int i=0; i<Kxb.ncol(); ++i )
	{
		Kxb( j, i ) = K( Xb[j]-xb[i], h );
	}

	//Ylty[j,i] = Y[j] <= y_seq[i]
	NumericMatrix Ylty = NumericMatrix(Dimension( Y.length(), y_seq.length() ));
	for ( int j=0; j<Ylty.nrow(); ++j )
	for ( int i=0; i<Ylty.ncol(); ++i )
	{
		Ylty( j, i ) = Y[j]<=y_seq[i] ? 1.0 : 0.0;
	}

	//fyxb[j,i] = sum_k K(Xb[k]-xb[j],h) * 1(Y[k]<=y_seq[i]) / sum_k K(Xb[k]-xb[j],h)
	NumericMatrix fyxb = NumericMatrix(Dimension( xb.length(), y_seq.length() ));
	for     ( int j=0; j<fyxb.nrow(); ++j )
	{
		double denom = 0.0;
		for ( int k=0; k<Y.length(); ++k )
		{
			denom += Kxb(k,j);
		}

		if ( denom == 0.0 ) [[unlikely]]
		{
			for ( int i=0; i<fyxb.ncol(); ++i )
			{
				fyxb( j, i ) = 0.0;
			}
		}
		else
		{
			for ( int i=0; i<fyxb.ncol(); ++i )
			{
				double numer = 0.0;
				for ( int k=0; k<Y.length(); ++k )
				{
					numer += Kxb(k,j) * Ylty(k,i);
				}
				fyxb( j, i ) = numer / denom;
			}
		}
	}

	return fyxb;
}
