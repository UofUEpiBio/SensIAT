library("assertthat")
library("inline")
library("pcoriRPackage")
library("Rcpp")
library("splines")
library("tidyverse") #`map`

load("testing.object.20240607.RData")

pcoriaccel_estimate_pmf_orig <- rcpp(
	sig = c(
		X = "numeric",
		Y = "numeric",
		xi = "numeric",
		y_seq = "numeric",
		h = "numeric"
	),
	body = "return pcoriaccel_estimate_pmf_orig( X,Y, Rcpp::as<double>(xi), y_seq, Rcpp::as<double>(h) );",
	includes = {c("
#include <cmath>
#include <vector>
#include <string>

constexpr double RT_TWOPI_RECIP =  0.398942280401432677;
constexpr double NEG_HALF_LOG2E = -0.7213475204444817  ; // -½ log₂e

template<class T> [[nodiscard]] constexpr
T sq( T val ) noexcept { return val*val; }
[[nodiscard]] constexpr double K_biweight2( double x, double h ) noexcept
{
	if ( std::abs(x) > h ) return 0.0;
	x /= h;
	return (15.0/16.0) * sq(1.0-sq(x));
}

//' Estimate the PMF directly with the K2_Biweight kernel.
//'
//' @param Xb vector (expected to be about 500 elements)
//' @param Y vector (same size as Xb)
//' @param xb vector
//' @param y_seq vector
//' @param h scalar bandwidth of kernel
//' @return estimated PMF
// [[Rcpp::export]]
[[nodiscard]] NumericVector pcoriaccel_estimate_pmf_orig(
	NumericVector X,
	NumericVector Y,
	double xi,
	NumericVector y_seq,
	double h
) noexcept {
	//auto K = K_normal;
	auto K = K_biweight2;
	//auto K = K_biweight4;

	/*
	Compute kernel applied to each pair of Xbⱼ and xᵢ

	Kxb[j,i] = K( Xb[j]-xb[i], h )

	Equivalent to `outer( Xb,xb, function(a,b){K(a-b,h)} )`
	*/
	NumericVector Kxb = NumericVector(X.length());
	for ( int j=0; j<Kxb.length(); ++j )
	{
		Kxb(j) = K( X[j]-xi, h );
	}

	// estimated_pmf[j] = sum_k K(Xb[k]-xi,h) * 1(Y[k] == y_seq[j]) / sum_k K(Xb[k]-xi,h)
	NumericVector estimated_pmf = NumericVector( y_seq.length() );
	double denom = 0.0;
	for ( int j=0; j<Kxb.length(); ++j )
	{
		denom += Kxb(j);

		auto i = std::find( y_seq.begin(), y_seq.end(), Y[j] );

		estimated_pmf(std::distance(y_seq.begin(),i)) += Kxb(j);
	}

	if ( denom == 0.0 ) [[unlikely]]
	{
		for ( int i=0; i<estimated_pmf.length(); ++i )
		{
			estimated_pmf(i) = 0.0;
		}
	}
	else
	{
		for ( int i=0; i<estimated_pmf.length(); ++i )
		{
			estimated_pmf(i) /= denom;
		}
	}
	return estimated_pmf;
}
")},

)

compute_influence_term_2_quadv_sim_via_matrix <-
function(
	X,                               #< Matrix of all covariates, transformed as necessary by model
	Y,                               #< Vector of all outcomes

	times,                           #< Vector of observation times for individual
	individual_X,                    #< Matrix of covariates for individual rows correspond to times
	                                 #  prepared for inferences for integration.
	x_slope,                         #< Vector of numeric(length(beta)) indicating how

	alpha,                           #< Vector of sensitivity parameters
	beta,                            #< Vector of coefficients of the outcome model

	spline_basis,                    #< Spline basis object (`orthogonalsplinebasis::SplineBasis`)

	bandwidth,                       #< Bandwidth for the kernel density estimate of the outcome model.
	tol = .Machine$double.eps^(1/4), #< Tolerance for integration

	...                              #< Additional arguments passed to the pcoriaccel_NW function, not implemented
) {
	# input validation checks
	assertthat::assert_that(
		is.matrix(X), is.vector(Y),
		ncol(X) == length(beta),
		ncol(X) == ncol(individual_X),
		nrow(X) == length(Y),
		length(x_slope) == ncol(individual_X),
		length(times) == nrow(individual_X),
		is( spline_basis, "SplineBasis" ),
		rlang::is_scalar_double(bandwidth),
		rlang::is_scalar_double(tol),
		is.numeric(alpha)
	)

	# obtain global integral bounds [a,b]
	a <- spline_basis@knots[[spline_basis@order]]
	b <- tail( spline_basis@knots, spline_basis@order )[[1]]

	# find time points determining sub integrals.
	period.times <- unique(c( a, times[a<times&times<b], b ))

	# Compute globally used terms.
	distinct.y <- sort( unique(Y) )
	Xb <- X %*% beta

	# compute integrals for each period.
	integrand <- \(period_ind,time)
	{
		# The coefficient vector for the given individual is determined by extending the last
		# observation to the current time point.  Time-dependent variables scale(time) and
		# scale(delta_time) are determined by the centering parameters which are global constants.
		# All other variables should be last observed values.  All this is encoded in the `x_slope`
		# object.
		x_at_time <- individual_X[period_ind,] + x_slope * (time-times[period_ind])

		# compute the pmf for Y at the given time point.
		pmf <- pcoriaccel_estimate_pmf_orig( Xb,Y, x_at_time%*%beta, distinct.y, h=bandwidth )

		# compute the expected value of the outcome at the given time point
		# given the sensitivity parameter alpha.
		eay <- exp(outer(distinct.y, alpha))
		numerator   <- crossprod(distinct.y*eay, pmf ) # sumₖ yₖ exp(αyₖ) pmfₖ
		denominator <- crossprod(           eay, pmf ) # sumₖ    exp(αyₖ) pmfₖ
		ev <- numerator/denominator

		# evaluate the spline basis functions at the given time point.
		B <- evaluate_basis(spline_basis, time)
		ret <- ev %*% B # returned value is a length(alpha) by ncol(B) matrix.
		return(ret)
	}
	#cat("Value at   ind , time   =  ",1,",",lower,"  :\n")
	#print(integrand(1,lower))

	period.integrals <-
		map(seq.int(length(period.times)-1), \(period.index){
			# Operating on period index.

			# lower and upper bound of the integral.
			lower <- period.times[period.index]
			upper <- period.times[period.index+1]

			# perform the integration through vectorized quadrature.
			int <- pracma::quadv( \(time){ integrand(period.index,time) }, lower,upper, tol=tol )
			#cat("Integral value   ind , time   =  ",period.index,",(",lower,",",upper,"):\n")
			#print(int$Q)
			return(int)
		})

	# Sum integrals of periods to get entire integral.
	map(period.integrals, getElement, "Q") %>%
		reduce(`+`) |>
		structure(    # this is collecting and carrying forward the convergence information of the period integrals.
			fcnt       = purrr::map_int( period.integrals, getElement, "fcnt"       ),
			estim.prec = purrr::map_dbl( period.integrals, getElement, "estim.prec" )
		)
}

evaluate_basis <- function(
	object, # Spline object which contains Matrices, an array of dimension
	x
) {
	# Extract components of the spline basis.
	#
	# /* M is a three dimensional array with dimension(
	#        * order,
	#        * length(knots)-order(i.e. number of basis functions),
	#        * length(knots)-2*order+1( Number of distinct internal intervals created by the knots)
	# )*/
	# const NumericArray M = object.slot('Matrices');
	M     <- object@Matrices
	# // Vector of spline knots
	# const NumericVector knots = object.slot('knots');
	knots <- object@knots
	# // Order of the spline, one more that the degree of the polynomial
	# // Could also be derived as dimension(M)[0]
	# const uint order =    object.slot('order');
	order <- object@order

	# Early kickout for out of bounds x
	if ( x<knots[order] || x>knots[length(knots)-order+1] ) return(rep(NA,dim(object)[2]))

	# Find the interval in which x lies
	if ( x == knots[length(knots)-order+1] )
		ind <- x <= knots # Special Case when x == upper bound of valid knots.
	else
		ind <- x <  knots
	if( all(ind) || all(!ind) )
	{
		if (x == knots[length(knots)-order+1] )
			return( rep(1,order) %*% matrix(M[,,dim(M)[3]],nrow=order) )
		else
			return( rep(0,ncol(object)) )
	}
	i <- which(ind)[1] - 1
	u <- (x-knots[i]) / (knots[i+1]-knots[i])
	U <- u^(0:(order-1))

	# auto i = order-1;
	# while(x>knots[i] && i < length(knots)-order-1) ++i;
	# u = (x-knots[i])/(knots[i+1]-knots[i])
	# U = NumericVector(u^0, u^1, ..., u^(order-1)
	# return( U %*% matrix(M[,,i-order+1], nrow=order))

	return( U %*% matrix(M[,,i-order+1],nrow=order) )
}



# test the function
alpha <- c(-0.5, 0, 0.5)

# Extract data from the fitted outcome model.
X <- model.matrix( terms(object$outcome.model), data=object$data )
Y <- model.response(model.frame( terms(object$outcome.model), data=object$data ))
beta <- coef(object$outcome.model)
bandwidth <- object$outcome.model$bandwidth

# The data for the integration must be manipulated for proper integration, and make sure that time
# and id are extracted as well.
integration_data <- rlang::inject(
	model.frame(
		terms(object$outcome.model),
		data = object$data |>
			dplyr::select(-!!object$variables$prev_outcome) |>
			dplyr::mutate( "{object$variables$prev_outcome}" := !!object$variables$outcome )
		,
		id   = !!object$variables$id,
		time = !!object$variables$time
	)
)



# From the adjusted data compute the covariate matrix for the integration.
integration_X <- model.matrix( object$outcome.model, data=integration_data )
assertthat::assert_that( nrow(integration_X) == nrow(integration_data) )

time <- integration_data[["(time)"]]
ids  <- integration_data[["(id)"  ]]

# Just picking an ID to use for testing and development any id works.
# All ids will need to be done.
set.seed(42)
id    <- sample( unique(ids), 1 )
times <- time[ ids == id ]

individual_X <- integration_X[ ids==id, ]

# Create the slope parameter from the centering statistics.
# This is model-specific, and will need to be generalized.
x_slope <- rep( 0, ncol(X) )
names(x_slope) <- colnames(X)
x_slope[4] <- 1/object$outcome.model.centering[[2]]
x_slope[5] <- 1/object$outcome.model.centering[[4]]


t0 <- Sys.time()
ref  = compute_influence_term_2_quadv_sim_via_matrix(
	X,Y, times,individual_X,x_slope, alpha,beta, object$base, bandwidth, 1e-8
)
t1 <- Sys.time()
test = pcoriaccel_compute_influence_term_2_quadv_sim_via_matrix(
	X,Y, times,individual_X,x_slope, alpha,beta, object$base, bandwidth
)
t2 <- Sys.time()

cat( "Reference result, in time ", t1-t0, ":\n", sep="" )
print(ref )
cat( "\nTest result, in time "   , t2-t1, ":\n", sep="" )
print(test)

cat( "\nSpeedup factor:", as.numeric(t1-t0)/as.numeric(t2-t1) )
