// [[Rcpp::plugins(cpp20)]]

#include "common.h"



//' Runs a *basic* implementation of the "NW" function with the "K2_Biweight" kernel, just as a
//' proof-of-concept.
//'
//' @param Xb    a vector (expected to be about 500 elements)
//' @param Y     a vector (same size as `Xb`)
//' @param xb    a vector
//' @param y_seq a vector
//' @param h     a scalar, the bandwidth of kernel
//'
//' @return A matrix of the same size as `xb` by `y_seq`.
//'
//' @export
// [[Rcpp::export]]
[[nodiscard]] NumericMatrix pcoriaccel_NW_basic(
	NumericVector Xb, NumericVector Y,
	NumericVector xb,
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
//' @param Xb    a vector (expected to be about 500 elements)
//' @param Y     a vector (same size as `Xb`)
//' @param xb    a vector
//' @param y_seq a vector
//' @param h     a scalar, the bandwidth of kernel
//' @param kernel a string, denoting the kernel function to use, either `"dnorm"`, `"K2_Biweight"`, or `"K4_Biweight"`
//'
//' @return A matrix of the same size as `xb` by `y_seq`.
//'
//' @export
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
