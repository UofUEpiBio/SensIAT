##### This function defined by MingYueh
SIDRnew <- function(X, Y,
                    Y.CP = NULL,
                    initial = NULL,
                    kernel = "K2_Biweight",
                    method = c("optim", "nlminb", "nmk"),
                    optim_method = "BFGS",
                    abs.tol = 1e-4,
                    bandwidth = NULL,
                    wi.boot = NULL)
{
    X <- as.matrix(X)
    Y <- as.matrix(Y)

    method <- match.arg(method)

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

        results <- list(coefficients = c(1, esti$par[1:(number_p-1)]),
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
        results <- list(coefficients = c(1, esti$par[1:(number_p-1)]),
                        bandwidth = bandwidth,
                        details = esti)
    }

    return(results)
}

SIDRnew_fixed_bandwidth <-
function(X, Y, ids,
         Y.CP = NULL,
         initial = NULL,
         kernel = "K2_Biweight",
         method = c("optim", "nlminb", "nmk"),
         optim_method = "BFGS",
         abs.tol = 1e-4
    ){
    X <- as.matrix(X)
    Y <- as.matrix(Y)

    method <- match.arg(method)

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

    if (kernel=="K2_Biweight")
        K <- function(x) 15/16*(1-(x)^2)^2 * (abs(x) <= 1)
    else if (kernel=="dnorm")
        K <- function(x) dnorm(x,0,1)
    else if (kernel=="K4_Biweight")
        K <- function(x) 105/64*(1-3*((x)^2))*(1-(x)^2)^2 * (abs(x) <= 1)
    else rlang::abort("Bad Kernel")

    Eij3 <- function(parameter){
        x <- c(X %*% parameter)
        y <- Y

        n <- length(y)
        yo <- order(y)
        ys <- y[yo]
        uy <- rle(ys)[[1]]
        cols <- cumsum(uy)
        ei <- rep(0, n)
        for (i in 1:n){
            Kih <- ifelse(ids==ids[i], 0, K(x-x[i]))
            denom <- sum(Kih)
            ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
        }
        return(sum(ei)/n^2)
    }


    if(method == "nlminb")
    {
        esti <- nlminb(start = initial,
                       objective = Eij3,
                       control = list(abs.tol = abs.tol))
    }else if (method == "optim")
    {
        # the new optimize function using optim, you can change the lower and upper
        esti <- optim(par = initial,
                      fn = Eij3,
                      method = optim_method,
                      control = list(abstol = abs.tol))
    }else if (method == "nmk")
    {
        assertthat::assert_that(requireNamespace("dfoptim", quietly = TRUE))
        esti <- dfoptim::nmk(par = initial,
                             fn = Eij3,
                             control = list(tol = abs.tol))
    }else rlang::abort("Invalid optimization method")
    results <- list(coefficients = esti$par,
                    bandwidth = 1,
                    details = esti)
}
