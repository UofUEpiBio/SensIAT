#' Outcome Modeler for `SensIAT` Single Index Model.
#'
#' @param formula The outcome model formula
#' @param data The data to fit the outcome model to.
#'             Should only include follow-up data, i.e. time > 0.
#' @param kernel The kernel to use for the outcome model.
#' @param method The optimization method to use for the outcome model, either `"optim"`, `"nlminb"`, or `"nmk"`.
#' @param id The patient identifier variable for the data.
#' @param initial Either a vector of initial values or a function to estimate initial values.
#'      If NULL (default), the initial values are estimated using the `MAVE::mave.compute` function.
#' @param ... Currently ignored, included for future compatibility.
#'
#' @return Object of class `SensIAT::Single-index-outcome-model` which contains the outcome model portion.
#' @export
#' @examples
#' \donttest{
#' # A basic example using fixed intensity bandwidth.
#' object <-
#'     fit_SensIAT_within_group_model(
#'         group.data = SensIAT_example_data,
#'         outcome_modeler = fit_SensIAT_single_index_fixed_bandwidth_model,
#'         id = Subject_ID,
#'         outcome = Outcome,
#'         time = Time,
#'         knots = c(60, 260, 460),
#'         End = 830,
#'         intensity.args = list(bandwidth = 30)
#'     )
#'
#' # A basic example using variable bandwidth but with fixed first coefficient.
#' object.bw <-
#'     fit_SensIAT_within_group_model(
#'         group.data = SensIAT_example_data,
#'         outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
#'         id = Subject_ID,
#'         outcome = Outcome,
#'         time = Time,
#'         knots = c(60, 260, 460),
#'         End = 830,
#'         intensity.args = list(bandwidth = 30)
#'     )
#' }
fit_SensIAT_single_index_fixed_coef_model <-
    function(formula, data, kernel = "K2_Biweight", method = "nmk", id = ..id.., initial = NULL, ...) {
        id <- ensym(id)
        mf <- rlang::inject(model.frame(formula, data = data, id = !!id))
        Xi <- model.matrix(formula, data = mf)

        Yi <- model.response(mf)

        if (is.null(initial)) {
            requireNamespace("MAVE", quietly = TRUE)
            initial <- coef(MAVE::mave.compute(Xi, Yi, max.dim = 1), 1)
        } else if (is.function(initial)) {
            initial <- initial(Xi, Yi)
        } else if (is.numeric(initial)) {
            initial <- initial
        } else {
            stop("initial must be a function, a numeric vector, or NULL")
        }

        if (initial[1] < 0) initial <- -initial

        val <- SIDR_Ravinew(
            X = Xi, Y = Yi, index_ID = mf[["(id)"]],
            initial = initial,
            kernel = kernel,
            method = method,
            ...
        )
        structure(
            append(
                val,
                list(
                    frame = mf,
                    data = data
                )
            ),
            class = c("SensIAT::outcome-model", "SensIAT::Single-index-outcome-model"),
            kernel = kernel,
            terms = terms(mf),
            id = id,
            initial = initial,
            restriction = "fixed_coef"
        )
    }

#' @export
`model.frame.SensIAT::Single-index-outcome-model` <-
    function(formula, data = NULL, ...) {
        if (is.null(data)) {
            data <- formula$data
        }
        NextMethod("model.frame", data = data, ...)
    }
#' @export
`model.matrix.SensIAT::Single-index-outcome-model` <-
    function(object, data = model.frame(object), ...) {
        model.matrix(terms(object), data = data, ...)
    }
#' @export
`formula.SensIAT::Single-index-outcome-model` <-
    function(x, ...) {
        as.formula(terms(x))
    }
#' @export
`coef.SensIAT::Single-index-outcome-model` <-
    function(object, ...) object$coef

# Helper to detect if running in knitr context
in_knitr <- function() {
    isTRUE(getOption("knitr.in.progress"))
}

#' @export
`print.SensIAT::Single-index-outcome-model` <-
    function(x, digits = max(3L, getOption("digits") - 3L), markdown = in_knitr(), ...) {
        cf <- coef(x)
        mm <- model.matrix(x)
        if (length(cf) == ncol(mm)) {
            names(cf) <- colnames(mm)
        }

        restriction <- attr(x, "restriction")
        restriction_label <- switch(restriction,
            fixed_coef = "Fixed first coefficient (= 1)",
            fixed_bandwidth = "Fixed bandwidth",
            if (!is.null(restriction)) restriction else "Unknown"
        )

        if (markdown) {
            cat("\n### Single-Index Outcome Model\n\n")
            cat("**Formula:** `", deparse(formula(x)), "`\n\n", sep = "")
            cat("| Property | Value |\n")
            cat("|:---------|:------|\n")
            cat("| Restriction | ", restriction_label, " |\n", sep = "")
            cat("| Kernel | ", attr(x, "kernel"), " |\n", sep = "")
            cat("| Bandwidth | ", format(x$bandwidth, digits = digits), " |\n\n", sep = "")
            cat("**Coefficients:**\n\n")
            cat("| Term | Estimate |\n")
            cat("|:-----|--------:|\n")
            for (i in seq_along(cf)) {
                nm <- if (!is.null(names(cf))) names(cf)[i] else paste0("[", i, "]")
                cat("| ", nm, " | ", format(cf[i], digits = digits), " |\n", sep = "")
            }
            cat("\n")
        } else {
            cat("\nSingle-Index Outcome Model\n\n")
            cat("Formula:", deparse(formula(x)), "\n")
            cat("Restriction:", restriction_label, "\n")
            cat("Kernel: ", attr(x, "kernel"), "\n")
            cat("Bandwidth:", format(x$bandwidth, digits = digits), "\n\n")
            cat("Coefficients:\n")
            print.default(format(cf, digits = digits), print.gap = 2L, quote = FALSE)
            cat("\n")
        }
        invisible(x)
    }

#' @export
`summary.SensIAT::Single-index-outcome-model` <-
    function(object, ...) {
        cf <- coef(object)
        # Get coefficient names from model matrix
        mm <- model.matrix(object)
        if (length(cf) == ncol(mm)) {
            names(cf) <- colnames(mm)
        }

        ans <- list(
            call = match.call(),
            formula = formula(object),
            coefficients = cf,
            bandwidth = object$bandwidth,
            kernel = attr(object, "kernel"),
            restriction = attr(object, "restriction"),
            nobs = nrow(object$frame),
            convergence = object$details$convergence,
            objective = if (!is.null(object$details$value)) object$details$value else object$details$objective
        )
        class(ans) <- "summary.SensIAT::Single-index-outcome-model"
        ans
    }

#' @export
`print.summary.SensIAT::Single-index-outcome-model` <-
    function(x, digits = max(3L, getOption("digits") - 3L), markdown = in_knitr(), ...) {
        restriction <- x$restriction
        restriction_label <- switch(restriction,
            fixed_coef = "Fixed first coefficient (= 1)",
            fixed_bandwidth = "Fixed bandwidth",
            if (!is.null(restriction)) restriction else "Unknown"
        )

        if (markdown) {
            cat("\n### Single-Index Outcome Model Summary\n\n")
            cat("**Formula:** `", deparse(x$formula), "`\n\n", sep = "")
            cat("| Property | Value |\n")
            cat("|:---------|:------|\n")
            cat("| Restriction | ", restriction_label, " |\n", sep = "")
            cat("| Kernel | ", x$kernel, " |\n", sep = "")
            cat("| Bandwidth | ", format(x$bandwidth, digits = digits), " |\n", sep = "")
            cat("| Observations | ", x$nobs, " |\n\n", sep = "")
            cat("**Coefficients:**\n\n")
            cat("| Term | Estimate |\n")
            cat("|:-----|--------:|\n")
            cf <- x$coefficients
            for (i in seq_along(cf)) {
                nm <- if (!is.null(names(cf))) names(cf)[i] else paste0("[", i, "]")
                cat("| ", nm, " | ", format(cf[i], digits = digits), " |\n", sep = "")
            }
            if (!is.null(x$convergence)) {
                cat("\n**Optimization:**\n\n")
                cat("| Property | Value |\n")
                cat("|:---------|:------|\n")
                cat("| Convergence | ", if (x$convergence == 0) "Yes" else paste("No (code", x$convergence, ")"), " |\n", sep = "")
                if (!is.null(x$objective)) {
                    cat("| Objective | ", format(x$objective, digits = digits), " |\n", sep = "")
                }
            }
            cat("\n")
        } else {
            cat("\nSingle-Index Outcome Model Summary\n")
            cat(rep("=", 40), "\n", sep = "")
            cat("\nFormula:", deparse(x$formula), "\n")
            cat("Restriction:", restriction_label, "\n")
            cat("Kernel: ", x$kernel, "\n")
            cat("Bandwidth:", format(x$bandwidth, digits = digits), "\n")
            cat("Observations:", x$nobs, "\n\n")

            cat("Coefficients:\n")
            print.default(format(x$coefficients, digits = digits), print.gap = 2L, quote = FALSE)

            if (!is.null(x$convergence)) {
                cat("\nOptimization:\n")
                cat("  Convergence:", if (x$convergence == 0) "Yes" else paste("No (code", x$convergence, ")"), "\n")
                if (!is.null(x$objective)) {
                    cat("  Objective value:", format(x$objective, digits = digits), "\n")
                }
            }
            cat("\n")
        }
        invisible(x)
    }

#' @export
`vcov.SensIAT::Single-index-outcome-model` <-
    function(object, ...) {
        # Single-index models do not have closed-form variance estimates.
        # Use jackknife() on the parent within_group_model for variance estimation.
        warning("vcov() is not available for Single-index-outcome-model. ",
                "Use jackknife() on the fitted SensIAT_within_group_model for variance estimation.",
                call. = FALSE)
        NULL
    }

#' @export
`predict.SensIAT::Single-index-outcome-model` <-
    function(object,
             newdata = NULL,
             type = c("lp", "response", "terms"),
             ...) {
        if (is.null(newdata)) newdata <- model.frame(object)
        type <- match.arg(type)

        frame <- model.frame(object, data = newdata)

        X_new <- model.matrix(terms(object), data = frame)

        if (type == "terms") {
            return(X_new)
        }

        lp <- X_new %*% object$coef
        if (type == "lp") {
            return(lp)
        }

        response <- vector("numeric", nrow(X_new))


        lp0 <- model.matrix(terms(object), object$data) %*% object$coef
        Y <- model.response(model.frame(object))
        y <- sort(unique(Y))

        for (i in 1:nrow(X_new)) {
            Fhat <- pcoriaccel_NW(
                Xb = lp0, Y = Y,
                xb = lp[i], y_seq = y,
                h = object$bandwidth,
                kernel = attr(object, "kernel")
            )
            pmf <- diff(c(0, Fhat))
            response[i] <- sum(y * pmf)
        }
        return(response)
    }

estimate_starting_coefficients <- function(X, Y, eps = 1e-7) {
    X <- as.matrix(X)
    Y <- as.vector(Y)

    number_n <- dim(X)[1]
    number_p <- dim(X)[2]

    Y.CP <- outer(Y, Y, "<=")

    # centralizing covariates
    X.cs <- t(t(X) - colMeans(X))

    # calculating m(y)=\E[X_i 1(Y_i\leq y)]
    # m.y <- (t(X)-colMeans(X)) %*% Y.CP/number_n
    # calculating K=\E[m(Y_i)m(Y_i)^T]
    m.y <- crossprod(X.cs, Y.CP) / number_n
    Km <- tcrossprod(m.y) / number_n

    eigen(solve(var(X) + eps * diag(number_p), Km))$vectors[, 1]
}


K2_Biweight_kernel <- function(x, h) {
    15 / 16 * (1 - (x / h)^2)^2 * (abs(x) <= h)
}
K4_Biweight_kernel <- function(x, h) {
    105 / 64 * (1 - 3 * ((x / h)^2)) * (1 - (x / h)^2)^2 * (abs(x) <= h)
}

NW_new <-
    function(Xb, Y, xb, y, h, kernel = "K2_Biweight") {
        if (kernel == "dnorm") {
            K <- function(x, h) {
                dnorm(x / h, 0, 1)
            } # Gaussian
        } else if (kernel == "K2_Biweight") {
            K <- function(x, h) {
                15 / 16 * (1 - (x / h)^2)^2 * (abs(x) <= h)
            } # K2_biweight
            # } else if(kernel=="K4_Biweight"){
            #     K <- function(x, h){105/64*(1-3*((x/h)^2))*(1-(x/h)^2)^2 * (abs(x) <= h) }# K4_biweight
        } else {
            stop("Kernel not recognized")
        }

        Kxb <- sapply(xb, function(x, Xb) K(Xb - x, h), Xb = Xb)

        Ylty <- sapply(y, function(x, Y) 1 * (Y <= x), Y = Y)

        denom <- colSums(Kxb)

        fyxb <- (denom != 0) * crossprod(Kxb, Ylty) / (denom + (denom == 0))

        return(fyxb)
    }

Cond_mean_fn_single2 <-
    function(alpha #< sensitivity parameter
             , X #< Matrix of covariates for all observations, including the spline basis as well as other covariates such as lag(time) and lag(outcome)
             , Y #< Outcome vector for all observations
             , x #< vector of covariates for the observation of interest
             , beta,
             bandwidth,
             ... #< for passing kernel forward
    ) {
        y <- sort(unique(Y))

        # conditional distribution
        # start <- Sys.time()
        Fhat_matrix <- pcoriaccel_NW(
            Xb = X %*% beta, Y = Y,
            xb = x %*% beta, y_seq = y,
            h = bandwidth,
            ...
        )
        # Convert matrix to vector (should be 1 row, length(y) columns)
        Fhat <- as.vector(Fhat_matrix)
        # end <- Sys.time()
        # end - start

        # density function
        Fhat1 <- c(0, Fhat[1:(length(y) - 1)])
        pmf <- Fhat - Fhat1

        # Question: Are we assuming Y is finite with support range_y or are we approximating an integral here?
        E_exp_alphaY <- sum(exp(alpha * y) * pmf)

        E_Yexp_alphaY <- sum(y * exp(alpha * y) * pmf)

        E_Y_past <- E_Yexp_alphaY / E_exp_alphaY

        return(list(
            E_Y_past = E_Y_past,
            E_exp_alphaY = E_exp_alphaY,
            E_Yexp_alphaY = E_Yexp_alphaY
        ))
    }


#' @export
`compute_SensIAT_expected_values.SensIAT::Single-index-outcome-model` <-
    function(model,
             alpha,
             # gamma,
             new.data = model.frame(model),
             ...) {
        assert_that(
            is(model, "SensIAT::Single-index-outcome-model"),
            is.numeric(alpha)
        )
        if (length(alpha) > 1) {
            return(
                purrr::map_dfr(
                    alpha,
                    `compute_SensIAT_expected_values.SensIAT::Single-index-outcome-model`,
                    model = model, new.data = new.data,
                    ...
                )
            )
        }
        if (nrow(new.data) == 0) {
            return(mutate(
                new.data,
                alpha = alpha,
                E_Y_past = numeric(0),
                E_exp_alphaY = numeric(0),
                E_Yexp_alphaY = numeric(0)
            ))
        }

        Xi <- model.matrix(terms(model), model$data)
        Yi <- model.response(model.frame(model))
        for (var in setdiff(all.vars(terms(model)), tbl_vars(new.data))) {
            new.data[[var]] <- NA
        }
        Xi_new <- model.matrix(terms(model), data = new.data)

        if (nrow(Xi_new) == 0) {
            return(mutate(
                new.data,
                alpha = alpha,
                E_Y_past = NA_real_,
                E_exp_alphaY = NA_real_,
                E_Yexp_alphaY = NA_real_
            ))
        }

        E_Y_past <- numeric(nrow(Xi_new))
        E_exp_alphaY <- numeric(nrow(Xi_new))
        E_Yexp_alphaY <- numeric(nrow(Xi_new))

        for (k in 1:nrow(Xi_new)) {
            # df_k <- new.data[k, ]
            # x = model.matrix(terms(model), data = df_k)
            temp <- Cond_mean_fn_single2(alpha,
                X = Xi,
                Y = Yi,
                x = Xi_new[k, , drop = FALSE],
                beta = model$coef,
                bandwidth = model$bandwidth,
                kernel = attr(model, "kernel")
            )

            E_Y_past[k] <- temp$E_Y_past
            E_exp_alphaY[k] <- temp$E_exp_alphaY
            E_Yexp_alphaY[k] <- temp$E_Yexp_alphaY
        }

        tibble(new.data, alpha, E_Y_past, E_Yexp_alphaY, E_exp_alphaY)
    }


#' @describeIn fit_SensIAT_single_index_fixed_coef_model for fitting with a fixed bandwidth
#' @export
fit_SensIAT_single_index_fixed_bandwidth_model <-
    function(formula, data, kernel = "K2_Biweight", method = "nmk", id = ..id..,
             initial = NULL, ...) {
        id <- ensym(id)
        mf <- rlang::inject(model.frame(formula, data = data, id = !!id))
        Xi <- model.matrix(formula, data = mf)

        Yi <- model.response(mf)

        force(initial)
        if (is.null(initial)) {
            requireNamespace("MAVE", quietly = TRUE)
            initial <- coef(MAVE::mave.compute(Xi, Yi, max.dim = 1), 1)
        } else if (is.function(initial)) {
            initial <- initial(Xi, Yi)
        } else if (is.numeric(initial)) {
            initial <- initial
        } else {
            stop("initial must be a function, a numeric vector, or NULL")
        }

        if (initial[1] < 0) initial <- -initial

        val <- SIDRnew_fixed_bandwidth(
            X = Xi, Y = Yi, ids = mf[["(id)"]],
            initial = initial,
            kernel = kernel,
            method = method,
            ...
        )
        structure(
            append(
                val,
                list(
                    frame = mf,
                    data = data
                )
            ),
            class = c("SensIAT::outcome-model", "SensIAT::Single-index-outcome-model"),
            kernel = kernel,
            id = id,
            terms = terms(mf),
            initial = initial,
            restriction = "fixed_bandwidth"
        )
    }
