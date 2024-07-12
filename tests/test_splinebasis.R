#source("test_splinebasis.R")

library("pcoriRPackage")
library("orthogonalsplinebasis")

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
		return(matrix( evaluate(basis,val), ncol=1 ))
	})
	t1 <- Sys.time()
	test_fn(function(val){
		return(matrix( pcoriaccel_evaluate_basis(basis,val), ncol=1 ))
	})
	t2 <- Sys.time()

	print( t1 - t0 )
	print( t2 - t1 )
}

basis <- SplineBasis(c(60,60,60,60,260,460,460,460,460))
test_basis(basis)
