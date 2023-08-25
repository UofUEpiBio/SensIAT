#################################### - change this part for single index model
#####   Inputs:  a mean mu and size theta, and sens parameter value alpha
#####   Output:  a list containing E[ Y exp(alpha Y ]/E[ exp(alpha Y] and E[ exp(alpha Y) ]
#####      for Y a truncated version of a negative binomial having mean mu and size theta

#' Compute conditional mean for single index model
#'
#' @param alpha
#' @param X
#' @param Y
#' @param x
#' @param beta
#' @param bandwidth
#'
#' @return List with elements
#'
#'   * `E(y)` and
#'   * `E(exp(alpha*y))`
#'
#' @export
#'
#' @examples
compute_conditional_mean_SIM <-
function(alpha, X, Y, x, beta, bandwidth){

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

    E_exp_alphaY <- sum( exp(alpha*y)*pmf )

    E_Yexp_alphaY <- sum( y*exp(alpha*y)*pmf )

    E_Y_past <- E_Yexp_alphaY/E_exp_alphaY

    return(list(
        'E(y)' = E_Y_past,
        'E(exp(alpha*y))' = E_exp_alphaY
    ))

}
