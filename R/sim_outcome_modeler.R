#' Outcome Modeler for PCORI Single Index Model.
#'
#' @param formula The outcome model formula
#' @param data The data to fit the outcome model to.
#' @param ... Currently ignored, included for future compatibility.
#'
#' @return TODO
#' @export
#'
#' @examples
PCORI_sim_outcome_modeler <- function(formula, data, kernel = "dnorm", method = "nmk", ...){
  mf <- model.frame(formula, data = data)
  Xi <- model.matrix(formula, data = data)

  Yi <-model.response(mf)

  SDR1 <- cumuSIR_new(X = Xi, Y = Yi)
  structure( SIDRnew(X = Xi, Y = Yi, initial = SDR1$basis[, 1], kernel = kernel, method = method),
             class = c('PCORI::Single-index-model'))
}

cumuSIR_new <- function(X, Y, eps = 1e-7)
{
    X <- as.matrix(X)
    Y <- as.matrix(Y)

    number_n <- dim(X)[1]
    number_p <- dim(X)[2]

    Y.CP <- matrix(
        Y[rep(1:number_n, times = number_n), ]<=
            Y[rep(1:number_n, each = number_n), ],
        nrow = number_n, ncol = number_n
    )

    # centralizing covariates
    X.cs <- t(t(X)-colMeans(X))

    # calculating m(y)=\E[X_i 1(Y_i\leq y)]
    m.y <- t(X.cs) %*% Y.CP/number_n
    # calculating K=\E[m(Y_i)m(Y_i)^T]
    Km <- m.y %*% t(m.y)/number_n

    Bhat <- eigen(solve(var(X) + eps*diag(number_p), Km))$vectors

    return(list(basis = Bhat))
}


##### This function defined by MingYueh
#' Multi-index distribution regression
#'
#' @param X
#' @param Y
#' @param Y.CP
#' @param initial
#' @param kernel It could be any of
#' `K2_Biweight`, `dnorm`, `K4_Biweight`, `K2_Biweight`, or `K4_Biweight`
#' @param bandwidth
#' @param wi.boot
#' @param fun.
#' @param method
#' @param optim_method It could be any of
#' `Nelder-Mead`, `BFGS`, `CG`, `L-BFGS-B`, or `Brent`
#' @param abs.tol
#'
#' @export
#'
SIDRnew <- function(X, Y,
                    Y.CP = NULL,
                    initial = NULL,
                    kernel = "dnorm",
                    method = "optim",
                    optim_method = "BFGS",
                    abs.tol = 1e-4,
                    bandwidth = NULL,
                    wi.boot = NULL)
{
    X <- as.matrix(X)
    Y <- as.matrix(Y)

    number_n <- dim(X)[1]
    number_p <- dim(X)[2]

    if (is.null(initial))
    {
        initial <- c(1, rep(0, number_p-1))
    }else
    {
        initial <- as.vector(initial)
        initial <- initial/initial[1]
    }

    if (is.null(bandwidth))
    {
        if (kernel=="K2_Biweight")
        {
            if (is.null(wi.boot))
            {
                Eij3 <- function(parameter){
                    K <- function(x, h) 15/16*(1-(x/h)^2)^2 * (abs(x) <= h) * (x!=0)

                    b <- c(1, parameter[1:(number_p-1)])
                    h <- exp(parameter[number_p])

                    x <- c(X%*%b)
                    y <- Y

                    n <- length(y)
                    yo <- order(y)
                    ys <- y[yo]
                    uy <- rle(ys)[[1]]
                    cols <- cumsum(uy)
                    ei <- rep(0, n)
                    for (i in 1:n){
                        Kih <- K(x-x[i],h=h)
                        denom <- sum(Kih)
                        ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
                    }
                    return(sum(ei)/n^2)
                }
                # cv.bh <- function(parameter)
                # {
                #   b <- c(1, parameter[1:(number_p-1)])
                #   h <- exp(parameter[number_p])
                #   cv <- mean((Y.CP-NWcv_K2B_rcpp(X = X %*% b, Y = Y.CP,
                #                                  h = h))^2)
                #   return(cv)
                # }
            }else
            {
                stop("There's no weighted version of the K2_Biweight kernel.")
            }
        }else if (kernel == "dnorm")
        {
            if (is.null(wi.boot))
            {
                Eij3 <- function(parameter){
                    K <- function(x, h) dnorm(x/h, 0, 1)* (x!=0)

                    b <- c(1, parameter[1:(number_p-1)])
                    h <- exp(parameter[number_p])

                    x <- c(X%*%b)
                    y <- Y

                    n <- length(y)
                    yo <- order(y)
                    ys <- y[yo]
                    uy <- rle(ys)[[1]]
                    cols <- cumsum(uy)
                    ei <- rep(0, n)
                    for (i in 1:n){
                        Kih <- K(x-x[i],h=h)
                        denom <- sum(Kih)
                        ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
                    }
                    return(sum(ei)/n^2)
                }
                # cv.bh <- function(parameter)
                # {
                #   b <- c(1, parameter[1:(number_p-1)])
                #   h <- exp(parameter[number_p])
                #   cv <- mean((Y.CP-NWcv_dnorm_rcpp(X = X %*% b, Y = Y.CP,
                #                                    h = h))^2)
                #   return(cv)
                # }
            }else
            {
                # wi.boot <- as.vector(wi.boot)
                # cv.bh <- function(parameter)
                # {
                #   b <- c(1, parameter[1:(number_p-1)])
                #   h <- exp(parameter[number_p])
                #   cv <- mean((Y.CP-NWcv_K2B_w_rcpp(X = X %*% b, Y = Y.CP,
                #                                    h = h, w = wi.boot))^2)
                #   return(cv)
                # }
                stop("There's no weighted version of the dnorm kernel.")
            }
        }else if (kernel=="K4_Biweight")
        {
            if (is.null(wi.boot))
            {
                Eij3 <- function(parameter){
                    K <- function(x, h) 105/64*(1-3*((x/h)^2))*(1-(x/h)^2)^2 * (abs(x) <= h) * (x!=0)

                    b <- c(1, parameter[1:(number_p-1)])
                    h <- exp(parameter[number_p])

                    x <- c(X%*%b)
                    y <- Y

                    n <- length(y)
                    yo <- order(y)
                    ys <- y[yo]
                    uy <- rle(ys)[[1]]
                    cols <- cumsum(uy)
                    ei <- rep(0, n)
                    for (i in 1:n){
                        Kih <- K(x-x[i],h=h)
                        denom <- sum(Kih)
                        ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
                    }
                    return(sum(ei)/n^2)
                }
                # cv.bh <- function(parameter)
                # {
                #   b <- c(1, parameter[1:(number_p-1)])
                #   h <- exp(parameter[number_p])
                #   cv <- mean((Y.CP-pmin(pmax(NWcv_K4B_rcpp(X = X %*% b, Y = Y.CP,
                #                                            h = h), 0), 1))^2)
                #   return(cv)
                # }
            }else
            {
                stop("There's no weighted version of the K4_Biweight kernel.")
            }
        }

        if(method == "nlminb")
        {
            esti <- nlminb(start = c(initial[-1], 0),
                           objective = Eij3,
                           control = list(abs.tol = abs.tol))
        }else if (method == "optim")
        {
            # the new optimize function using optim, you can change the lower and upper
            esti <- optim(par = c(initial[-1], 0),
                          fn = Eij3,
                          method = optim_method,
                          control = list(abstol = abs.tol))
        }else if (method == "nmk")
        {
            assertthat::assert_that(requireNamespace("dfoptim", quietly = TRUE))
            esti <- dfoptim::nmk(par = c(initial[-1], 0),
                        fn = Eij3,
                        control = list(tol = abs.tol))
        }

        results <- list(coef = c(1, esti$par[1:(number_p-1)]),
                        bandwidth = exp(esti$par[number_p]),
                        details = esti)
    }else
    {
        if (kernel=="K2_Biweight")
        {
            if (is.null(wi.boot))
            {
                Eij3 <- function(parameter){
                    K <- function(x, h) 15/16*(1-(x/h)^2)^2 * (abs(x) <= h) * (x!=0)

                    b <- c(1, parameter[1:(number_p-1)])
                    h <- bandwidth

                    x <- c(X%*%b)
                    y <- Y

                    n <- length(y)
                    yo <- order(y)
                    ys <- y[yo]
                    uy <- rle(ys)[[1]]
                    cols <- cumsum(uy)
                    ei <- rep(0, n)
                    for (i in 1:n){
                        Kih <- K(x-x[i],h=h)
                        denom <- sum(Kih)
                        ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
                    }
                    return(sum(ei)/n^2)
                }

                # cv.b <- function(parameter)
                # {
                #   b <- c(1, parameter[1:(number_p-1)])
                #   cv <- mean((Y.CP-NWcv_K2B_rcpp(X = X %*% b, Y = Y.CP,
                #                                  h = bandwidth))^2)
                #   return(cv)
                # }
            }else
            {
                stop("There's no weighted version of the K2_Biweight kernel.")
            }
        }else if (kernel=="dnorm")
        {
            if (is.null(wi.boot))
            {
                Eij3 <- function(parameter){
                    K <- function(x, h) dnorm(x/h,0,1) * (x!=0)

                    b <- c(1, parameter[1:(number_p-1)])
                    h <- bandwidth

                    x <- c(X%*%b)
                    y <- Y

                    n <- length(y)
                    yo <- order(y)
                    ys <- y[yo]
                    uy <- rle(ys)[[1]]
                    cols <- cumsum(uy)
                    ei <- rep(0, n)
                    for (i in 1:n){
                        Kih <- K(x-x[i],h=h)
                        denom <- sum(Kih)
                        ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
                    }
                    return(sum(ei)/n^2)
                }

                # cv.b <- function(parameter)
                # {
                #   b <- c(1, parameter[1:(number_p-1)])
                #   cv <- mean((Y.CP-NWcv_dnorm_rcpp(X = X %*% b, Y = Y.CP,
                #                                    h = bandwidth))^2)
                #   return(cv)
                # }
            }else
            {
                # wi.boot <- as.vector(wi.boot)
                # cv.b <- function(parameter)
                # {
                #   b <- c(1, parameter[1:(number_p-1)])
                #   cv <- mean((Y.CP-NWcv_K2B_w_rcpp(X = X %*% b, Y = Y.CP,
                #                                    h = bandwidth, w = wi.boot))^2)
                #   return(cv)
                # }
                stop("There's no weighted version of the dnorm kernel.")
            }
        }else if (kernel=="K4_Biweight")
        {
            if (is.null(wi.boot))
            {
                Eij3 <- function(parameter){
                    K <- function(x, h) 105/64*(1-3*((x/h)^2))*(1-(x/h)^2)^2 * (abs(x) <= h) * (x!=0)

                    b <- c(1, parameter[1:(number_p-1)])
                    h <- bandwidth

                    x <- c(X%*%b)
                    y <- Y

                    n <- length(y)
                    yo <- order(y)
                    ys <- y[yo]
                    uy <- rle(ys)[[1]]
                    cols <- cumsum(uy)
                    ei <- rep(0, n)
                    for (i in 1:n){
                        Kih <- K(x-x[i],h=h)
                        denom <- sum(Kih)
                        ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
                    }
                    return(sum(ei)/n^2)
                }
                # cv.b <- function(parameter)
                # {
                #   b <- c(1, parameter[1:(number_p-1)])
                #   cv <- mean((Y.CP-pmin(pmax(NWcv_K4B_rcpp(X = X %*% b, Y = Y.CP,
                #                                            h = bandwidth), 0), 1))^2)
                #   return(cv)
                # }
            }else
            {
                stop("There's no weighted version of the K4_Biweight kernel.")
            }
        }

        if(method == "nlminb")
        {
            esti <- nlminb(start = initial[-1],
                           objective = Eij3,
                           control = list(abs.tol = abs.tol))
        }else if (method == "optim")
        {
            # the new optimize function using optim, you can change the lower and upper
            esti <- optim(par = initial[-1],
                          fn = Eij3,
                          method = optim_method,
                          control = list(abstol = abs.tol))
        }else if (method == "nmk")
        {
            assertthat::assert_that(requireNamespace("dfoptim", quietly = TRUE))
            esti <- dfoptim::nmk(par = initial[-1],
                        fn = Eij3,
                        control = list(tol = abs.tol))
        }else rlang::abort("Invalid optimization method")
        results <- list(coef = c(1, esti$par[1:(number_p-1)]),
                        bandwidth = bandwidth,
                        details = esti)
    }

    return(results)
}


NW_new <- function(Xb, Y, xb, y, h, kernel = "dnorm"){

    if(kernel == "dnorm"){
        K <- function(x, h){dnorm(x/h, 0, 1)} # Gaussian
    } else if(kernel == "K2_Biweight"){
        K <- function(x, h){15/16*(1-(x/h)^2)^2 * (abs(x) <= h)} # K2_biweight
    } else if(kernel=="K4_Biweight"){
        K <- function(x, h){105/64*(1-3*((x/h)^2))*(1-(x/h)^2)^2 * (abs(x) <= h) }# K4_biweight
    }

    Kxb <- sapply(xb, function(x, Xb) K(Xb-x, h), Xb=Xb)

    Ylty <- sapply(y, function(x, Y) 1*(Y <= x), Y=Y)

    denom <- colSums(Kxb)

    fyxb <- (denom!=0)*crossprod(Kxb, Ylty)/(denom + (denom==0))

    return(fyxb)

}

