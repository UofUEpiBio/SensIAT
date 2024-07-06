library("pcoriRPackage")
library("pracma")



profile <- function( name, integrand, lo,hi, tol=.Machine$double.eps^(1/2) )
{
	t0 <- Sys.time()
	ref  <- quadv( integrand, lo,hi, tol )
	t1 <- Sys.time()
	test <- pcoriaccel_integrate_simp( integrand, lo,hi, tol )
	t2 <- Sys.time()

	cat("============================================\n\n")
	cat(name)
	cat("\n")
	cat( "Reference result, in time ", t1-t0, ":\n", sep="" )
	print(ref )
	cat( "Test result, in time "     , t2-t1, ":\n", sep="" )
	print(test)
	cat( "Speedup factor: ", as.numeric(t1-t0)/as.numeric(t2-t1), "\n\n", sep="" )
}



#Compare with:
#	https://rdrr.io/cran/pracma/man/quadv.html
#	Reference seems a little more accurate, but test is 20⨯ faster.
profile( "⟨sin(x),cos(x)⟩, 0≤x≤π", \(x){ c( sin(x), cos(x) ) }, 0,pi )

#Compare with:
#	https://rdrr.io/cran/pracma/man/quad.html
#	Reference `quadv` takes 21 samples and gets a very wrong answer; for some reason it doesn't
#	detect it needs more.  Test takes 358 samples and gets the right answer (1.28213), and does it
#	in ≈1/8th the time.
profile( "x*cos(exp(x)/10)*sin(πexp(x)/10), 0≤x≤4", \(x){ x * cos(0.1*exp(x)) * sin(0.1*pi*exp(x)) }, 0,4 )
