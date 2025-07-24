#' Single Index Model using MAVE and Optimizing Bandwidth.
#'
#' Single index model estimation using minimum average variance estimation (MAVE).
#' A direction is estimated using MAVE, and then the bandwidth is selected by
#' minimization of the cross-validated pseudo-integrated squared error.
#' Optionally, the initial coefficients of the outcome model can be re-estimated
#' by optimization on a spherical manifold.  This option requires the
#' [ManifoldOptim][ManifoldOptim::manifold.optim] package.
#'
#' @param formula The outcome model formula
#' @param data The data to fit the outcome model to.
#'             Should only include follow-up data, i.e. time > 0.
#' @param kernel The kernel to use for the outcome model.
#' @param mave.method The method to use for the MAVE estimation.
#' @param id The patient identifier variable for the data.
#' @param bw.selection The criteria for bandwidth selection, either `'ise'` for Integrated Squared Error or `'mse'` for Mean Squared Error.
#' @param bw.method The method for bandwidth selection, either `'optim'` for using optimization or `'grid'` for grid search.
#' @param reestimate.coef Logical indicating whether to re-estimate the coefficients of the outcome model after bandwidth selection.
#' @param bw.range A numeric vector of length 2 indicating the range of bandwidths to consider for selection as a multiple of the standard deviation of the single index predictor.
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
         bw.method = c('optim', 'grid', 'optimize'),
         bw.range = c(0.01, 1.5),
         reestimate.coef = FALSE,
         ...
){
    id <- ensym(id)
    mf <- rlang::inject(model.frame(formula, data = data, id = !!id))
    X <- model.matrix(formula, data = mf)
    Y <- model.response(mf)
    ids <- mf[['(id)']]
    unique_y <- sort(unique(Y))

    if (kernel=="K2_Biweight"){
        K <- function(x) 15/16*(1-(x)^2)^2 * (abs(x) <= 1)
        K1<- function(x) 3/4*(1-(x)^2) * (abs(x) <= 1)
    } else if (kernel=="dnorm") {
        K <- function(x) dnorm(x,0,1)
        K1 <- function(x) -x*dnorm(x, 0, 1)
    } else if (kernel=="K4_Biweight"){
        K <- function(x) 105/64*(1-3*((x)^2)) * (1-(x)^2)^2 * (abs(x) <= 1)
    } else{
        stop("Unknown kernel type. Please use either 'K2_Biweight', 'dnorm', or 'K4_Biweight'.")
    }

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

    bw.method <- match.arg(bw.method)
    sigma <- sd(Xbeta)
    if(bw.method == 'grid'){
        # Use grid search for bandwidth selection
        log_bw_seq <- log(seq(min(bw.range), max(bw.range), length.out = 100) * sigma)
        err_values <- purrr::map_dbl(log_bw_seq, err)
        bw_opt <- log_bw_seq[which.min(err_values)]
        bw.details <- list(par = bw_opt,
                       value = min(err_values),
                       convergence = 0,
                       message = "Grid search completed successfully.",
                       bw.method = bw.method,
                       log_bw_seq = log_bw_seq,
                       err_values = err_values
                       )
    } else if(bw.method == 'optim'){
        # Use optimization for bandwidth selection
        initial <- log(sigma * 0.30)
        bw.details <- optim(initial, err, method = "L-BFGS-B", lower = log(sigma * min(bw.range)), upper = log(sigma * max(bw.range)), ...)
        bw.details$initial = initial
        bw.details$bw.method = bw.method
        bw_opt <- bw.details$par
    } else if(bw.method == 'optimize'){
        # Use optimize for bandwidth selection
        result <- stats::optimize(err,
                           interval = c(log(sigma * min(bw.range)), log(sigma * max(bw.range))),
                           ...)
        bw_opt <- result$minimum
        bw.details <- list(minimum = result$minimum,
                           value = result$objective,
                           bw.method  = bw.method,
                           interval = c(log(sigma * 0.05), log(sigma * 1.5)))
    } else{
        stop("Unknown bw.method type. Please use either 'optim' or 'grid'.")
    }

    if(!reestimate.coef) {
        return(
            structure(
                list(
                    coefficients = beta_hat,
                    bandwidth = exp(bw_opt),
                    details = bw.details,
                    frame = mf,
                    data = data
                )
                , class = c('SensIAT::outcome-model', 'SensIAT::Single-index-outcome-model')
                , id = id
                , kernel = kernel
                , terms = terms(mf)
            )
        )
    }

    if(bw.selection == "ise") {
        objFun <- function(beta){
            Xbeta <- (X %*% beta)[,]
            D <- outer(Xbeta, Xbeta, FUN = `-`)
            W <- K(D/exp(bw_opt)) * Id_neq
            denom <- rowSums(W)
            Fhat <- sweep(W %*% Imat, 1, denom, "/")
            Fhat[is.nan(Fhat)] <- 0
            return (sum((Imat - Fhat)^2)/length(Fhat)^2)
        }
    } else if(bw.selection == "mse") {
        objFun <- function(beta){
            Xbeta <- (X %*% beta)[,]
            D <- outer(Xbeta, Xbeta, FUN = `-`)
            W <- K1(D/exp(bw_opt)) * Id_neq
            denom <- rowSums(W)
            mean_Y <- c(W %*% Y) / denom
            return (mean((Y - mean_Y)^2, na.rm = TRUE))
        }
    } else {
        stop("Unknown bw.selection type. Please use either 'ise' or 'mse'.")
    }
    rlang::check_installed("ManifoldOptim")
    requireNamespace("ManifoldOptim", quietly = TRUE)
    mod <- Rcpp::Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
    prob <- new(mod$RProblem, objFun)
    mani.params <- ManifoldOptim::get.manifold.params(IsCheckParams = FALSE)
    solver.params <- ManifoldOptim::get.solver.params(IsCheckParams = FALSE)
    maniDef <- ManifoldOptim::get.sphere.defn(length(beta_hat))

    res <- ManifoldOptim::manifold.optim(prob, maniDef,
                          mani.params = mani.params,
                          solver.params = solver.params, x0 = beta_hat)
    beta_hat <- res$xopt

    return(
        structure(
            list(
                coefficients = as.vector(res$xopt),
                bandwidth = exp(bw_opt),
                details = bw.details,
                frame = mf,
                data = data,
                details.refit = res
            )
            , class = c('SensIAT::outcome-model', 'SensIAT::Single-index-outcome-model')
            , id = id
            , kernel = kernel
            , terms = terms(mf)
        )
    )

}
