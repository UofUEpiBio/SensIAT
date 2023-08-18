#####  Function for drawing a value of Y (support = 0, 1/6, 2/6, ..., 6)

#####  Inputs:  a prediction from a negative binomial regression model on the
#####     'response' scale mu, and the theta (size) parameter from this model

#####  Output:  1/6 times a draw from a negative binomial distribution with mean mu
#####     and size theta, truncated to have support 0,1,2,..., 36
Y_draw_fn_single <- function(X, Y, x, beta, bandwidth){

    y <- seq(0, 6, by = 1/6)
    # conditional distribution
    Fhat <- NW(X = X %*% beta, Y = Y, x = x %*% beta,
               regression = "distribution",
               kernel = "dnorm",
               y = y,
               bandwidth = bandwidth)

    # density function
    Fhat1 <- c(0, Fhat[1:(length(y) - 1)])
    pmf <- Fhat - Fhat1

    Y_new <- sample(y, size=1, prob=pmf)
    return(Y_new)

}
