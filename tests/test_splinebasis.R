#source("test_splinebasis.R")

library("pcoriRPackage")

load("testing.object.20240607.RData")

evaluate_basis <- function(
		object, # Spline object which contains Matrices, an array of dimension
		x
	){

	# Extract components of the spline basis.
	#
	# /* M is a three dimensional array with dimension(
	#    * order,
	#    * length(knots)-order(i.e. number of basis functions),
	#    * length(knots)-2*order+1( Number of distinct internal intervals created by the knots)
	# )*/
	# const NumericArray M = object.slot('Matrices');
	M    <-object@Matrices
	# // Vector of spline knots
	# const NumericVector knots = object.slot('knots');
	knots<-object@knots
	# // Order of the spline, one more that the degree of the polynomial
	# // Could also be derived as dimension(M)[0]
	# const uint order =  object.slot('order');
	order<-object@order

	# Early kickout for out of bounds x
	if(x < knots[order] | x>knots[length(knots)-order+1])return(rep(NA,dim(object)[2]))

	# Find the interval in which x lies
	if(x==knots[length(knots)-order+1]){
		# Special Case when x == upper bound of valid knots.
		ind <- x<=knots
	} else {
		ind <- x<knots
	}
	if(all(ind)|all(!ind))  {
		if(x==knots[length(knots)-order+1])
			return(rep(1,order)%*%matrix(M[,,dim(M)[3]],nrow=order))
		else
			return(rep(0,ncol(object)))
	}
	i<-which(ind)[1]-1
	u<-(x-knots[i])/(knots[i+1]-knots[i])
	U <-u^(0:(order-1))

	# auto i = order-1;
	# while(x>knots[i] && i < length(knots)-order-1) ++i;
	# u = (x-knots[i])/(knots[i+1]-knots[i])
	# U = NumericVector(u^0, u^1, ..., u^(order-1)
	# return( U %*% matrix(M[,,i-order+1], nrow=order))

	return(U%*%matrix(M[,,i-order+1],nrow=order))
}

test_basis <- function(basis){
	vals <- c( 59.0,60.0,61.0, 259.0,260.0,261.0, 459.0,460.0,461.0 )

	test_fn <- function(fn){
		print(cbind(
			fn(vals[[1]]), fn(vals[[2]]), fn(vals[[3]]),
			fn(vals[[4]]), fn(vals[[5]]), fn(vals[[6]]),
			fn(vals[[7]]), fn(vals[[8]]), fn(vals[[9]])
		))
	}

	t0 <- Sys.time()
	test_fn(function(val){
		return(matrix( evaluate_basis(basis,val), ncol=1 ))
	})
	t1 <- Sys.time()
	test_fn(function(val){
		return(matrix( pcoriaccel_evaluate_basis(basis,val), ncol=1 ))
	})
	t2 <- Sys.time()
	print( t1 - t0 )
	print( t2 - t1 )
}

test_basis(object$base)
