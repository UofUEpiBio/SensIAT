SIDR_Ravinew <- function(X, Y,
                         Y.CP = NULL,
                         initial = NULL,
                         kernel = "dnorm",
                         method = "optim",
                         optim_method = "BFGS",
                         abs.tol = 1e-4,
                         bandwidth = NULL,
                         wi.boot = NULL,
                         index_ID)
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
        if(initial[1] != 0)
            initial <- initial/initial[1]
    }

    if (is.null(bandwidth))
    {
        if (kernel=="K2_Biweight")
        {
            if (is.null(wi.boot))
            {
                Eij3 <- function(parameter){
                    K <- function(x, h) 15/16*(1-(x/h)^2)^2 * (abs(x) <= h)

                    b <- c(1, parameter[1:(number_p-1)])
                    h <- exp(parameter[number_p])

                    x <- c(X%*%b)
                    y <- Y

                    n <- length(y)
                    yo <- order(y)
                    ys <- y[yo]
                    uy <- rle(ys)[[1]]
                    cols <- cumsum(uy)  # Ravi
                    ei <- rep(0, n)
                    for (i in 1:n){
                        Kih <- K(x-x[i],h=h)

                        # remove all obs of ith obs's patient
                        index_remove <- which(index_ID == index_ID[i])
                        Kih[index_remove] <- 0

                        # Kih[i] <- 0                       # the fix
                        denom <- sum(Kih)
                        ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
                    }
                    return(sum(ei)/n^2)
                }
            }else
            {
                stop("There's no weighted version of the K2_Biweight kernel.")
            }
        }else if (kernel == "dnorm")
        {
            if (is.null(wi.boot))
            {
                Eij3 <- function(parameter){
                    K <- function(x, h) dnorm(x/h, 0, 1)

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

                        # remove all obs of ith obs's patient
                        index_remove <- which(index_ID == index_ID[i])
                        Kih[index_remove] <- 0

                        # Kih[i] <- 0       # the fix
                        denom <- sum(Kih)
                        ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
                    }
                    return(sum(ei)/n^2)
                }
            }else
            {
                stop("There's no weighted version of the dnorm kernel.")
            }
        # }else if (kernel=="K4_Biweight")
        # {
        #     if (is.null(wi.boot))
        #     {
        #         Eij3 <- function(parameter){
        #             K <- function(x, h) 105/64*(1-3*((x/h)^2))*(1-(x/h)^2)^2 * (abs(x) <= h)
        #
        #             b <- c(1, parameter[1:(number_p-1)])
        #             h <- exp(parameter[number_p])
        #
        #             x <- c(X%*%b)
        #             y <- Y
        #
        #             n <- length(y)
        #             yo <- order(y)
        #             ys <- y[yo]
        #             uy <- rle(ys)[[1]]
        #             cols <- cumsum(uy)
        #             ei <- rep(0, n)
        #             for (i in 1:n){
        #                 Kih <- K(x-x[i],h=h)
        #
        #                 # remove all obs of ith obs's patient
        #                 index_remove <- which(index_ID == index_ID[i])
        #                 Kih[index_remove] <- 0
        #
        #                 # Kih[i] <- 0                       # the fix
        #                 denom <- sum(Kih)
        #                 ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
        #             }
        #             return(sum(ei)/n^2)
        #         }
        #     }else
        #     {
        #         stop("There's no weighted version of the K4_Biweight kernel.")
        #     }
        } else
            rlang::abort("The kernel is not supported.")

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
                    K <- function(x, h) 15/16*(1-(x/h)^2)^2 * (abs(x) <= h)

                    b <- c(1, parameter[1:(number_p-1)])
                    h <- bandwidth

                    x <- c(X%*%b)
                    y <- Y

                    n <- length(y)
                    yo <- order(y)
                    ys <- y[yo]
                    uy <- rle(ys)[[1]]
                    cols <- cumsum(uy)  # Ravi
                    ei <- rep(0, n)
                    for (i in 1:n){
                        Kih <- K(x-x[i],h=h)

                        # remove all obs of ith obs's patient
                        index_remove <- which(index_ID == index_ID[i])
                        Kih[index_remove] <- 0

                        # Kih[i] <- 0                       # the fix
                        denom <- sum(Kih)
                        ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
                    }
                    return(sum(ei)/n^2)
                }

            }else
            {
                stop("There's no weighted version of the K2_Biweight kernel.")
            }
        }else if (kernel=="dnorm")
        {
            if (is.null(wi.boot))
            {
                Eij3 <- function(parameter){
                    K <- function(x, h) dnorm(x/h,0,1)

                    b <- c(1, parameter[1:(number_p-1)])
                    h <- bandwidth

                    x <- c(X%*%b)
                    y <- Y

                    n <- length(y)
                    yo <- order(y)
                    ys <- y[yo]
                    uy <- rle(ys)[[1]]
                    cols <- cumsum(uy)  # Ravi
                    ei <- rep(0, n)
                    for (i in 1:n){
                        Kih <- K(x-x[i],h=h)

                        # remove all obs of ith obs's patient
                        index_remove <- which(index_ID == index_ID[i])
                        Kih[index_remove] <- 0

                        # Kih[i] <- 0                       # the fix
                        denom <- sum(Kih)
                        ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
                    }
                    return(sum(ei)/n^2)
                }
            }else
            {
                stop("There's no weighted version of the dnorm kernel.")
            }
        # }else if (kernel=="K4_Biweight")
        # {
        #     if (is.null(wi.boot))
        #     {
        #         Eij3 <- function(parameter){
        #             K <- function(x, h) 105/64*(1-3*((x/h)^2))*(1-(x/h)^2)^2 * (abs(x) <= h)
        #
        #             b <- c(1, parameter[1:(number_p-1)])
        #             h <- bandwidth
        #
        #             x <- c(X%*%b)
        #             y <- Y
        #
        #             n <- length(y)
        #             yo <- order(y)
        #             ys <- y[yo]
        #             uy <- rle(ys)[[1]]
        #             cols <- cumsum(uy)  # Ravi
        #             ei <- rep(0, n)
        #             for (i in 1:n){
        #                 Kih <- K(x-x[i],h=h)
        #
        #                 # remove all obs of ith obs's patient
        #                 index_remove <- which(index_ID == index_ID[i])
        #                 Kih[index_remove] <- 0
        #
        #                 # Kih[i] <- 0                       # the fix
        #                 denom <- sum(Kih)
        #                 ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
        #             }
        #             return(sum(ei)/n^2)
        #         }
        #     }else
        #     {
        #         stop("There's no weighted version of the K4_Biweight kernel.")
        #     }
        } else
            rlang::abort("The kernel is not supported.")

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
        }
        results <- list(coef = c(1, esti$par[1:(number_p-1)]),
                        bandwidth = bandwidth,
                        details = esti)
    }

    return(results)
}
