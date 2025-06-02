#' Single Index Model using Mave and Optimizing Bandwidth.
#'
#' @param formula The outcome model formula
#' @param data The data to fit the outcome model to.
#'             Should only include follow-up data, i.e. time > 0.
#' @param kernel The kernel to use for the outcome model.
#' @param mave.method The method to use for the MAVE estimation.
#' @param id The patient identifier variable for the data.
#' @param bw.selection The method for bandwidth selection, either 'ise' for Integrated Squared Error or 'mse' for Mean Squared Error.
#' @param ... Additional arguments to be passed to [optim].
#'
#' @return Object of class `SensIAT::Single-index-outcome-model` which contains the outcome model portion.
#' @export
SensIAT_sim_outcome_modeler_mave <-
function(formula, data,
         kernel = "K2_Biweight",
         mave.method = "meanMAVE",
         id = ..id..,
         bw.selection = c('ise', 'mse'),
         ...
){
    id <- ensym(id)
    mf <- rlang::inject(model.frame(formula, data = data, id = !!id))
    X <- model.matrix(formula, data = mf)
    Y <- model.response(mf)
    ids <- mf[['(id)']]
    unique_y <- sort(unique(Y))

    if (kernel=="K2_Biweight")
        K <- function(x) 15/16*(1-(x)^2)^2 * (abs(x) <= 1)
    else if (kernel=="dnorm")
        K <- function(x) dnorm(x,0,1)
    else if (kernel=="K4_Biweight")
        K <- function(x) 105/64*(1-3*((x)^2))*(1-(x)^2)^2 * (abs(x) <= 1)
    else
        stop("Unknown kernel type. Please use either 'K2_Biweight', 'dnorm', or 'K4_Biweight'.")


    assertthat::assert_that(requireNamespace("MAVE", quietly = TRUE))
    mave_fit <- MAVE::mave.compute(X, Y, max.dim = 1, method = mave.method)

    beta_hat <- mave_fit$dir[[1]][,1,drop=TRUE]
    Xbeta <- (X %*% beta_hat)[,]
    D <- outer(Xbeta, Xbeta, FUN = `-`)
    Id_neq <- outer(ids, ids, \(x,y)as.numeric(x!=y))
    Imat <- outer(Y, unique_y, `<=`)

    bw.selection <- match.arg(bw.selection)
    if(bw.selection == 'ise'){
        # Use Integrated Squared Error (ISE) for bandwidth selection
        err <- function(log_bandwidth){
            W <- K(D/exp(log_bandwidth)) * Id_neq
            denom <- rowSums(W)
            Fhat <- sweep(W %*% Imat, 1, denom, "/")
            Fhat[is.nan(Fhat)] <- 0
            return (sum((Imat - Fhat)^2)/length(Fhat)^2)
        }
    } else if(bw.selection == 'mse'){
        # Use Mean Squared Error (MSE) for bandwidth selection
        err <- function(log_bandwidth){
            W <- K(D/exp(log_bandwidth)) * Id_neq
            denom <- rowSums(W)
            mean_Y <- c(W %*% Y) / denom
            return (mean((Y - mean_Y)^2, na.rm = TRUE))
        }
    } else {
        stop("Unknown bw.selection type. Please use either 'ise' or 'mse'.")
    }

    # bw_opt <- optimize(err,
    #                    interval = c(log(.Machine$double.xmin), log(.Machine$double.xmax)),
    #                    ...,
    #                    lower = log(.Machine$double.xmin),
    #                    upper = log(.Machine$double.xmax)
    #                    )
    bw_opt <- optim(0, err, method = "BFGS",  ...)
    structure(
        list(
            coefficients = beta_hat,
            bandwidth = exp(bw_opt$par),
            details = bw_opt,
            frame = mf,
            data = data
        )
        , class = c('SensIAT::outcome-model', 'SensIAT::Single-index-outcome-model')
        , kernel = kernel
        , terms = terms(mf)
    )
}
