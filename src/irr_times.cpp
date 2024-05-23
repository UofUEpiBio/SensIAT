#include <cmath>
#include <vector>

#include <Rcpp.h>



using namespace Rcpp;



/*
When the R-facing API changes, regenerate by running, from R in the project root directory:

	library("Rcpp")
	compileAttributes(".")

	R -e 'library("Rcpp");compileAttributes(".")' #should work but doesn't?
*/



constexpr double RT_TWOPI_RECIP =  0.398942280401432677;
constexpr double NEG_HALF_LOG2E = -0.7213475204444817  ; // -½ log₂e

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
//' @param Xb vector (expected to be about 500 elements)
//' @param Y vector (same size as Xb)
//' @param xb vector
//' @param y_seq vector
//' @param h scalar bandwidth of kernel
//' @return Matrix fyxb
// [[Rcpp::export]]
[[nodiscard]] NumericMatrix pcoriaccel_NW_basic(
	NumericVector Xb, NumericVector Y,
	NumericVector xb,
	NumericVector y_seq,
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



template< class Tfloat, int kernel > [[nodiscard]] inline static
NumericMatrix _pcoriaccel_NW(
	NumericVector const& Xb, NumericVector const& Y,
	NumericVector const& xb,
	NumericVector const& y_seq,
	Tfloat h
) noexcept {
	/*
	The original code constructs a matrix `Kxb`, which is a kernel function applied to a simple
	function, and `Ylty`, a boolean matrix.  The code effectively computes a slightly fancy matrix
	multiplication.  Let `J` be the matrix of 1s of the same size as `Ylty`.  Then each element of
	the result is as the element in `Kxbᵀ Ylty` divided by the element in `Kxbᵀ J`.

	Besides the intrinsic benefits related to rewriting in C++ over R (compile-time optimizations,
	higher performance all around, etc.), we also optimize a lot explicitly here.

	The quotient divides out any constants in the kernel function, so we can omit them for more
	performance and accuracy.  Constructing the matrices, especially `Ylty`, explicitly in memory is
	wasteful.  We construct each row of `Kxbᵀ` in memory, and compute each element of `Ylty` and `J`
	on the fly (they are trivial).  The simpler loops open the door to vectorization (currently only
	autovectorization) and potential future multithreading.
	*/

	//Constant used for kernel function
	Tfloat C;
	if constexpr ( kernel == 1 ) C=NEG_HALF_LOG2E/sq(h);
	else                         C=sq(h);

	//Temporary row: a row of `Kxbᵀ`
	std::vector<Tfloat> KxbT_row( Xb.length() );

	//Result
	NumericMatrix fyxb = NumericMatrix(Dimension( xb.length(), y_seq.length() ));

	//Calculation
	for ( int j=0; j<fyxb.nrow(); ++j )
	{
		//Compute row of `Kxbᵀ`
		for ( int i=0; i<Xb.length(); ++i )
		{
			Tfloat x = (Tfloat)Xb[i] - (Tfloat)xb[j];

			//Omit constants in kernel functions; see above
			if      constexpr ( kernel == 1 ) //"dnorm"
			{
				//Note also `exp2(⋯)`; usually slightly faster than `exp(⋯)` (in theory and in practice)
				KxbT_row[i] = std::exp2( C * sq(x) );
			}
			else if constexpr ( kernel == 2 ) //"K2_Biweight"
			{
				if ( std::abs(x) >= h ) KxbT_row[i]=0;
				else                    KxbT_row[i]=sq( C - sq(x) );
			}
			else if constexpr ( kernel == 3 ) //"K4_Biweight"
			{
				if ( std::abs(x) >= h ) KxbT_row[i]=0;
				else
				{
					Tfloat y = sq(x);
					KxbT_row[i] = (C-3*y) * sq(C-y);
				}
			}
		}
		//Compute sum of that row
		//	TODO: maybe better to move into previous loop with manual vectorization
		Tfloat denom = 0;
		for ( int i=0; i<Xb.length(); ++i ) denom+=KxbT_row[i];

		//(Row of answer is zeros if denominator is zero)
		if ( denom == 0 ) [[unlikely]]
		{
			for ( int i=0; i<fyxb.ncol(); ++i )
			{
				fyxb( j, i ) = 0;
			}

			continue;
		}
		//(Otherwise...)

		Tfloat denom_recip = 1 / denom; //Repeated multiplication is faster than repeated division
		for ( int i=0; i<fyxb.ncol(); ++i )
		{
			//Compute element via definition of matrix multiplication
			//	TODO: branchy; avoid with mask and manual vectorization?
			Tfloat numer = 0;
			for ( int k=0; k<Y.length(); ++k )
			{
				if ( Y[k] > y_seq[i] ) continue;
				numer += KxbT_row[k];
			}

			fyxb( j, i ) = (double)( numer * denom_recip );
		}
	}

	//Done
	return fyxb;
}

//' Runs an optimized implementation of the "NW" function.
//'
//' @param Xb vector (expected to be about 500 elements)
//' @param Y vector (same size as Xb)
//' @param xb vector
//' @param y_seq vector
//' @param h scalar bandwidth of kernel
//' @return Matrix fyxb
// [[Rcpp::export]]
[[nodiscard]] NumericMatrix pcoriaccel_NW(
	NumericVector Xb, NumericVector Y, //Same length, about 500
	NumericVector xb,
	NumericVector y_seq, //unique values
	double h, //scalar bandwidth
	String kernel = "K2_Biweight"
) {
	if ( Xb.length() != Y.length() ) [[unlikely]]
	{
		stop("Lengths of arguments `Xb` and `Y` must match!");
	}

	using Tfloat = double;
	if      ( kernel == "dnorm" )
	{
		return _pcoriaccel_NW< Tfloat, 1 >( Xb,Y, xb, y_seq, (Tfloat)h );
	}
	else if ( kernel == "K2_Biweight" )
	{
		return _pcoriaccel_NW< Tfloat, 2 >( Xb,Y, xb, y_seq, (Tfloat)h );
	}
	else if ( kernel == "K4_Biweight" ) [[likely]]
	{
		return _pcoriaccel_NW< Tfloat, 3 >( Xb,Y, xb, y_seq, (Tfloat)h );
	}
	else [[unlikely]]
	{
		stop("Invalid value for `kernel`: choices are { \"dnorm\", \"K2_Biweight\", \"K4_Biweight\" }.");
	}
}
