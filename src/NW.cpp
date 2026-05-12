// [[Rcpp::plugins(cpp20)]]

#include "common.h"
#include <numeric> // std::iota



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
//' @keywords internal
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

		if ( denom == 0.0 ) // [[unlikely]]
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
	int const N = Xb.length();
	int const Q = xb.length();
	int const M = y_seq.length();

	// Constant used for kernel function (constants cancel in ratio, so omitted for biweight)
	Tfloat C;
	if constexpr ( kernel == 1 ) C = NEG_HALF_LOG2E / sq(h);
	else                          C = sq(h);

	// Pre-sort training indices by Y value once for all queries.
	// Enables O(N+M) CDF computation per query (merge scan) instead of O(N*M).
	// y_seq is guaranteed sorted (= sort(unique(Y)) from R), so we can advance a
	// single pointer through sort_order_Y rather than re-scanning N points per y.
	std::vector<int> sort_order_Y(N);
	std::iota(sort_order_Y.begin(), sort_order_Y.end(), 0);
	std::sort(sort_order_Y.begin(), sort_order_Y.end(),
	          [&](int a, int b){ return Y[a] < Y[b]; });

	// For compact-support kernels (biweight): pre-sort Xb indices so we can use
	// binary search to find the active window [xb[j]-h, xb[j]+h] and skip all
	// training points with zero weight.
	std::vector<int>    sort_order_Xb;
	std::vector<Tfloat> sorted_Xb_vals;
	if constexpr ( kernel == 2 || kernel == 3 ) {
		sort_order_Xb.resize(N);
		std::iota(sort_order_Xb.begin(), sort_order_Xb.end(), 0);
		std::sort(sort_order_Xb.begin(), sort_order_Xb.end(),
		          [&](int a, int b){ return Xb[a] < Xb[b]; });
		sorted_Xb_vals.resize(N);
		for ( int i = 0; i < N; ++i )
			sorted_Xb_vals[i] = (Tfloat)Xb[sort_order_Xb[i]];
	}

	// Kernel weights, reused across queries
	std::vector<Tfloat> KxbT_row(N, (Tfloat)0);

	// Result
	NumericMatrix fyxb = NumericMatrix(Dimension(Q, M));

	for ( int j = 0; j < Q; ++j )
	{
		Tfloat xb_j  = (Tfloat)xb[j];
		Tfloat denom = 0;

		if constexpr ( kernel == 2 || kernel == 3 ) {
			// Clear weights from previous query
			std::fill(KxbT_row.begin(), KxbT_row.end(), (Tfloat)0);

			// Binary search: only iterate points inside the support [xb_j-h, xb_j+h].
			// Points at exactly ±h yield weight 0 and are harmless to include.
			int lo = (int)( std::lower_bound(sorted_Xb_vals.begin(), sorted_Xb_vals.end(), xb_j - h) - sorted_Xb_vals.begin() );
			int hi = (int)( std::upper_bound(sorted_Xb_vals.begin(), sorted_Xb_vals.end(), xb_j + h) - sorted_Xb_vals.begin() );

			// Kernel + denom in one pass over the active window only
			for ( int idx = lo; idx < hi; ++idx ) {
				int    i = sort_order_Xb[idx];
				Tfloat x = (Tfloat)Xb[i] - xb_j;
				Tfloat w;
				if constexpr ( kernel == 2 ) {
					w = sq( C - sq(x) );
				} else {
					Tfloat y2 = sq(x);
					w = (C - 3*y2) * sq(C - y2);
				}
				KxbT_row[i]  = w;
				denom       += w;
			}
		} else {
			// Gaussian: full scan; merge denom accumulation into kernel loop
			for ( int i = 0; i < N; ++i ) {
				Tfloat x = (Tfloat)Xb[i] - xb_j;
				// exp2 is slightly faster than exp in practice
				KxbT_row[i]  = std::exp2( C * sq(x) );
				denom       += KxbT_row[i];
			}
		}

		if ( denom == 0 ) // [[unlikely]]
		{
			for ( int i = 0; i < M; ++i ) fyxb(j, i) = 0;
			continue;
		}

		Tfloat denom_recip = (Tfloat)1 / denom;

		// CDF merge scan: O(N+M) instead of O(N*M).
		// Walk sort_order_Y (Y ascending) and y_seq (also ascending) together.
		// Advance ptr while the next training Y is <= the current CDF threshold.
		Tfloat cumsum = 0;
		int    ptr    = 0;
		for ( int i = 0; i < M; ++i )
		{
			while ( ptr < N && Y[sort_order_Y[ptr]] <= y_seq[i] ) {
				cumsum += KxbT_row[sort_order_Y[ptr]];
				++ptr;
			}
			// Clamp to [0,1]: denom and cumsum sum the same values in different orders
			// (index order vs Y-sorted order), so floating-point rounding can push the
			// ratio infinitesimally outside [0,1] for kernels with full support (dnorm).
			double val = (double)( cumsum * denom_recip );
			fyxb(j, i) = val < 0.0 ? 0.0 : (val > 1.0 ? 1.0 : val);
		}
	}

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
//' @keywords internal
// [[Rcpp::export]]
[[nodiscard]] NumericMatrix pcoriaccel_NW(
	NumericVector Xb, NumericVector Y, //Same length, about 500
	NumericVector xb,
	NumericVector y_seq, //unique values
	double h, //scalar bandwidth
	String kernel = "K2_Biweight"
) {
	if ( Xb.length() != Y.length() ) // [[unlikely]]
	{
		stop("Lengths of arguments `Xb` and `Y` must match!");
	}

	using Tfloat = double;
	if      ( kernel == "dnorm" )
	{
		return _pcoriaccel_NW< Tfloat, 1 >( Xb,Y, xb, y_seq, (Tfloat)h );
	}
	else if ( kernel == "K2_Biweight" ) // [[likely]]
	{
		return _pcoriaccel_NW< Tfloat, 2 >( Xb,Y, xb, y_seq, (Tfloat)h );
	}
	// else if ( kernel == "K4_Biweight" )
	// {
	// 	return _pcoriaccel_NW< Tfloat, 3 >( Xb,Y, xb, y_seq, (Tfloat)h );
	// }
	else // [[unlikely]]
	{
		stop("Invalid value for `kernel`: choices are { \"dnorm\", \"K2_Biweight\" }.");
	}
}

// ---------------------------------------------------------------------------
// Vectorised conditional-expectation function
// ---------------------------------------------------------------------------
// Computes E[exp(alpha*Y)|xb[j]] and E[Y*exp(alpha*Y)|xb[j]] for all j
// in a single C++ call without materialising the Q×M CDF matrix.
//
// Key differences from pcoriaccel_NW:
//  - exp(alpha*Y[k]) is pre-computed once for all training points
//  - The weighted sums are accumulated directly (no merge-scan over y_seq)
//  - Returns a Q×2 matrix: col 0 = E_exp_alphaY, col 1 = E_Yexp_alphaY
//
// Complexity per query: O(K) for K2_Biweight (K = active-window size), O(N) for dnorm

template< class Tfloat, int kernel > [[nodiscard]] inline static
NumericMatrix _pcoriaccel_NW_expectations(
	NumericVector const& Xb,
	NumericVector const& Y,
	NumericVector const& xb,
	Tfloat alpha,
	Tfloat h
) noexcept {
	int const N = Xb.length();
	int const Q = xb.length();

	// Kernel constant (same convention as _pcoriaccel_NW)
	Tfloat C;
	if constexpr ( kernel == 1 ) C = NEG_HALF_LOG2E / sq(h);
	else                          C = sq(h);

	// Pre-compute exp(alpha*Y[k]) and Y[k]*exp(alpha*Y[k]) once for all queries.
	std::vector<Tfloat> eY(N), YeY(N);
	for ( int k = 0; k < N; ++k ) {
		Tfloat ey = std::exp( alpha * (Tfloat)Y[k] );
		eY[k]  = ey;
		YeY[k] = (Tfloat)Y[k] * ey;
	}

	// For compact-support kernels: pre-sort Xb once to enable binary-search active window.
	std::vector<int>    sort_order_Xb;
	std::vector<Tfloat> sorted_Xb_vals;
	if constexpr ( kernel == 2 || kernel == 3 ) {
		sort_order_Xb.resize(N);
		std::iota(sort_order_Xb.begin(), sort_order_Xb.end(), 0);
		std::sort(sort_order_Xb.begin(), sort_order_Xb.end(),
		          [&](int a, int b){ return Xb[a] < Xb[b]; });
		sorted_Xb_vals.resize(N);
		for ( int i = 0; i < N; ++i )
			sorted_Xb_vals[i] = (Tfloat)Xb[sort_order_Xb[i]];
	}

	// Result: Q rows × 2 cols  [E_exp_alphaY, E_Yexp_alphaY]
	NumericMatrix result( Q, 2 );

	for ( int j = 0; j < Q; ++j ) {
		Tfloat xb_j = (Tfloat)xb[j];
		Tfloat denom = 0, numer_exp = 0, numer_Yexp = 0;

		if constexpr ( kernel == 2 || kernel == 3 ) {
			int lo = (int)( std::lower_bound(sorted_Xb_vals.begin(), sorted_Xb_vals.end(), xb_j - h) - sorted_Xb_vals.begin() );
			int hi = (int)( std::upper_bound(sorted_Xb_vals.begin(), sorted_Xb_vals.end(), xb_j + h) - sorted_Xb_vals.begin() );
			for ( int idx = lo; idx < hi; ++idx ) {
				int    i = sort_order_Xb[idx];
				Tfloat x = (Tfloat)Xb[i] - xb_j;
				Tfloat w;
				if constexpr ( kernel == 2 ) {
					w = sq( C - sq(x) );
				} else {
					Tfloat y2 = sq(x);
					w = (C - 3*y2) * sq(C - y2);
				}
				denom      += w;
				numer_exp  += w * eY[i];
				numer_Yexp += w * YeY[i];
			}
		} else {
			// Gaussian: full scan
			for ( int i = 0; i < N; ++i ) {
				Tfloat x = (Tfloat)Xb[i] - xb_j;
				Tfloat w = std::exp2( C * sq(x) );
				denom      += w;
				numer_exp  += w * eY[i];
				numer_Yexp += w * YeY[i];
			}
		}

		if ( denom == 0 ) {
			result(j, 0) = 0.0;
			result(j, 1) = 0.0;
		} else {
			result(j, 0) = (double)( numer_exp  / denom );
			result(j, 1) = (double)( numer_Yexp / denom );
		}
	}

	return result;
}

//' Compute conditional expectations E[exp(alpha*Y)|xb] and E[Y*exp(alpha*Y)|xb]
//'
//' Faster alternative to calling [pcoriaccel_NW()] and post-processing: computes
//' both expectations in a single pass without materialising the Q x M CDF matrix.
//'
//' @param Xb     Numeric vector of training linear predictors (length N)
//' @param Y      Numeric vector of training outcomes (length N)
//' @param xb     Numeric vector of query linear predictors (length Q)
//' @param alpha  Sensitivity parameter (scalar)
//' @param h      Bandwidth (scalar)
//' @param kernel Kernel name: `"dnorm"` or `"K2_Biweight"`
//'
//' @return A Q x 2 numeric matrix. Column 1 is E[exp(alpha*Y)|xb[j]],
//'   column 2 is E[Y*exp(alpha*Y)|xb[j]], for j = 1..Q.
//'
//' @keywords internal
// [[Rcpp::export]]
[[nodiscard]] NumericMatrix pcoriaccel_NW_expectations(
	NumericVector Xb,
	NumericVector Y,
	NumericVector xb,
	double alpha,
	double h,
	String kernel = "K2_Biweight"
) {
	if ( Xb.length() != Y.length() )
		stop("Lengths of arguments `Xb` and `Y` must match!");

	using Tfloat = double;
	if      ( kernel == "dnorm" )
		return _pcoriaccel_NW_expectations< Tfloat, 1 >( Xb, Y, xb, (Tfloat)alpha, (Tfloat)h );
	else if ( kernel == "K2_Biweight" )
		return _pcoriaccel_NW_expectations< Tfloat, 2 >( Xb, Y, xb, (Tfloat)alpha, (Tfloat)h );
	else
		stop("Invalid value for `kernel`: choices are { \"dnorm\", \"K2_Biweight\" }.");
}
