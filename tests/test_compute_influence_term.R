library("assertthat")
requireNamespace("inline")
library("pcoriRPackage")
library("Rcpp")
library("splines")

object <-
    fit_PCORI_within_group_model(
        group.data = PCORI_example_data,
        outcome_modeler = PCORI_sim_outcome_modeler,
        alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
        id.var = Subject_ID,
        outcome.var = Outcome,
        time.var = Time,
        End = 830,
        knots = c(60,60,60,60,260,460,460,460,460),
        control = pcori_control('quadv')
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
		pmf <- pcoriaccel_estimate_pmf( Xb,Y, x_at_time%*%beta, distinct.y, h=bandwidth )

		# compute the expected value of the outcome at the given time point
		# given the sensitivity parameter alpha.
		eay <- exp(outer(distinct.y, alpha))
		numerator   <- crossprod(distinct.y*eay, pmf ) # sumₖ yₖ exp(αyₖ) pmfₖ
		denominator <- crossprod(           eay, pmf ) # sumₖ    exp(αyₖ) pmfₖ
		ev <- numerator/denominator

		# evaluate the spline basis functions at the given time point.
		B <- orthogonalsplinebasis::evaluate(spline_basis, time)
		ret <- ev %*% B # returned value is a length(alpha) by ncol(B) matrix.
		return(ret)
	}
	#cat("Value at   ind , time   =  ",1,",",lower,"  :\n")
	#print(integrand(1,lower))

	period.integrals <-
		purrr::map(seq.int(length(period.times)-1), \(period.index){
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
	purrr::map(period.integrals, getElement, "Q") %>%
		reduce(`+`) |>
		structure(    # this is collecting and carrying forward the convergence information of the period integrals.
			fcnt       = purrr::map_int( period.integrals, getElement, "fcnt"       ),
			estim.prec = purrr::map_dbl( period.integrals, getElement, "estim.prec" )
		)
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
