NW_test <- function( Xb,Y, xb, y, h, kernel="K2_Biweight" )
{
    # select kernel
    if      ( kernel == "dnorm"       )
    {
        K <- function(x,h) { dnorm(x/h,0,1) } # Gaussian
    }
    else if ( kernel == "K2_Biweight" )
    {
        K <- function(x,h) { 15/16 * (1-(x/h)^2)^2 * (abs(x)<=h) } # K2_biweight
    }
    else if ( kernel=="K4_Biweight"   )
    {
        K <- function(x,h) { 105/64 * (1-3*((x/h)^2)) * (1-(x/h)^2)^2 * (abs(x)<=h) } # K4_biweight
    }

    # Compute the kernel applied to each pair of Xb_i and x_j
    # equivalent to outer(Xb, xb, \(a,b){K(a-b, h)})
    # Kxb[i,j] = K(Xb[i] - xb[j], h)
    Kxb <- sapply( xb, function(x,Xb) K(Xb-x,h), Xb=Xb )

    # Ylty[i,k] = Y[i] <= y[k]
    Ylty <- sapply( y, function(x,Y) 1*(Y<=x), Y=Y )

    # denom[j] = sum_i K(Xb[i] - xb[j], h)
    denom <- colSums(Kxb)

    # fyxb[j,k] = sum_i K(Xb[i] - xb[j], h) * 1(Y[i] <= y[k]) / sum_i K(Xb[i] - xb[j], h)
    fyxb <- (denom!=0) * crossprod(Kxb,Ylty) / (denom+(denom==0))

    return(fyxb)
}
test_that("NW/K2_Biweight", {


    # Define the test data
    set.seed(20240522)
    n <- 100
    Xb <- rnorm(n)
    Y <- round( rnorm(n), 1 )
    xb <- seq( -2.5, 2.5, 0.5 )
    y <- sort( unique(Y) )

    # parameters
    n_y <- length(y)
    n_xb <- length(xb)

    Fhat  <- NW_test(Xb,Y, xb, y, h=0.6, "K2_Biweight" )
    Fhat2 <- pcoriaccel_NW( Xb,Y, xb, y, h=0.6, "K2_Biweight" )

    expect_equal(dim(Fhat2), c(n_xb,n_y))
    expect_true( all( apply( Fhat2, 1, \(f){all(head(f, -1) <= tail(f, -1))} ) ))
    expect_true(all( Fhat2 >= 0 ))
    expect_true(all( Fhat2 <= 1 ))
    expect_equal(Fhat2[,ncol(Fhat2)], rep(1, nrow(Fhat2)) )
    expect_equal(Fhat, Fhat2)
})
# test_that("NW/K4_Biweight", {
#
#
#     # Define the test data
#     set.seed(20240522)
#     n <- 100
#     Xb <- rnorm(n)
#     Y <- round( rnorm(n), 1 )
#     xb <- seq( -2.5, 2.5, 0.5 )
#     y <- sort( unique(Y) )
#
#     # parameters
#     n_y <- length(y)
#     n_xb <- length(xb)
#
#     Fhat  <- NW_test(Xb,Y, xb, y, h=0.6, "K4_Biweight" )
#     Fhat2 <- pcoriaccel_NW( Xb,Y, xb, y, h=0.6, "K4_Biweight" )
#
#     expect_equal(dim(Fhat2), c(n_xb,n_y))
#     expect_true( all( apply( Fhat2, 1, \(f){all(head(f, -1) <= tail(f, -1))} ) ))
#     expect_true(all( Fhat2 >= 0 ))
#     expect_true(all( Fhat2 <= 1 ))
#     expect_equal(Fhat2[,ncol(Fhat2)], rep(1, nrow(Fhat2)) )
#     expect_equal(Fhat, Fhat2)
# })
test_that("NW/dnorm", {


    # Define the test data
    set.seed(20240522)
    n <- 100
    Xb <- rnorm(n)
    Y <- round( rnorm(n), 1 )
    xb <- seq( -2.5, 2.5, 0.5 )
    y <- sort( unique(Y) )

    # parameters
    n_y <- length(y)
    n_xb <- length(xb)

    Fhat  <- NW_test(Xb,Y, xb, y, h=0.6, "dnorm" )
    Fhat2 <- pcoriaccel_NW( Xb,Y, xb, y, h=0.6, "dnorm" )

    expect_equal(dim(Fhat2), c(n_xb,n_y))
    expect_true( all( apply( Fhat2, 1, \(f){all(head(f, -1) <= tail(f, -1))} ) ))
    expect_true(all( Fhat2 >= 0 ))
    expect_true(all( Fhat2 <= 1 ))
    expect_equal(Fhat2[,ncol(Fhat2)], rep(1, nrow(Fhat2)) )
    expect_equal(Fhat, Fhat2)
})
