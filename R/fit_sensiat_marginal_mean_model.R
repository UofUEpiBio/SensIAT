#' Fit the Marginal Means Model
#'
#' @param data Data for evaluation of the model. Should match the data used to fit the intensity and outcome models.
#' @param id The subject identifier variable in the data. Lazy evaluation is used, so it can be a symbol or a string.
#' @param alpha Sensitivity parameter, a vector of values.
#' @param knots Location of spline knots. If a `SplineBasis` object is provided, it is used directly.
#' @param intensity.model The assessment time intensity model.
#' @param outcome.model The observed effects model.
#' @param spline.degree The degree of the spline basis, default is 3 (cubic splines).
#' @param ... Additional arguments passed to `compute_influence_terms`.
#'
#' @returns a list with the fitted model, including the coefficients and their variances for each alpha value.
#' @export
#'
#' @examples
#' # Note: example takes approximately 30 seconds to run.
#' \donttest{
#' library(survival)
#' library(dplyr)
#' library(splines)
#' # Create followup data with lags
#' # added variables `..prev_time..`, `..delta_time..` and `..prev_outcome..`
#' # have special interpretations when computing the influence.
#' data_with_lags <- SensIAT_example_data |>
#'     dplyr::group_by(Subject_ID) |>
#'     dplyr::mutate(
#'         ..prev_outcome.. = dplyr::lag(Outcome, default = NA_real_, order_by = Time),
#'         ..prev_time.. = dplyr::lag(Time, default = 0, order_by = Time),
#'         ..delta_time.. = Time - dplyr::lag(.data$Time, default = NA_real_, order_by = Time)
#'     )
#'
#' # Create the observation time intensity model
#' intensity.model <-
#'     coxph(Surv(..prev_time.., Time, !is.na(Outcome)) ~ ..prev_outcome.. + strata(Visit),
#'         data = data_with_lags |> dplyr::filter(.data$Time > 0)
#'     )
#'
#' # Create the observed outcome model
#' outcome.model <-
#'     fit_SensIAT_single_index_fixed_coef_model(
#'         Outcome ~ ns(..prev_outcome.., df = 3) + ..delta_time.. - 1,
#'         id = Subject_ID,
#'         data = data_with_lags |> filter(Time > 0)
#'     )
#'
#' # Fit the marginal outcome model
#' mm <- fit_SensIAT_marginal_mean_model(
#'     data = data_with_lags,
#'     id = Subject_ID,
#'     alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'     knots = c(60, 260, 460),
#'     intensity.model = intensity.model,
#'     time.vars = c("..delta_time.."),
#'     outcome.model = outcome.model
#' )
#' }
fit_SensIAT_marginal_mean_model <-
    function(data,
             id,
             alpha,
             knots,
             outcome.model,
             intensity.model,
             spline.degree = 3L,
             ...) {
        if (is(knots, "SplineBasis")) {
            base <- knots
        } else {
            knots <- c(
                rep(head(knots, 1), spline.degree),
                knots,
                rep(tail(knots, 1), spline.degree)
            )
            base <- SplineBasis(knots, order = spline.degree + 1L)
        }
        V_inverse <- solve(GramMatrix(base))


        needed.vars <- c(
            all.vars(terms(intensity.model)),
            all.vars(terms(outcome.model))
        )
        stopifnot(all(needed.vars %in% names(data)))

        for_each_alpha <- \(a){
            compute_influence_terms(
                data = data, id = {{ id }},
                base = base,
                alpha = a,
                outcome.model = outcome.model,
                intensity.model = intensity.model,
                ...
            )
        }
        influence.terms <- purrr::map(alpha, for_each_alpha)

        # Results
        Beta <- map(influence.terms, \(IT){
            uncorrected.beta_hat <- (colSums(IT$term1) + colSums(IT$term2)) / length(IT$id)
            estimate <- as.vector(V_inverse %*% uncorrected.beta_hat)
            variance <- tcrossprod(V_inverse %*% (t(IT$term1 + IT$term2) - uncorrected.beta_hat)) / c(length(IT$id)^2)
            list(estimate = estimate, variance = variance)
        })


        structure(
            list(
                models = list(
                    intensity = intensity.model,
                    outcome = outcome.model
                ),
                data = data,
                influence = influence.terms,
                alpha = alpha,
                coefficients = map(Beta, getElement, "estimate"),
                coefficient.variance = map(Beta, getElement, "variance"),
                influence.args = list(...),
                base = base,
                V_inverse = V_inverse
            ),
            class = "SensIAT_marginal_mean_model",
            call = match.call(expand.dots = TRUE),
            link = "identity",
            loss = "lp_mse"
        )
    }

# Utility methods for SensIAT_marginal_mean_model ----------------------------

#' @export
coef.SensIAT_marginal_mean_model <- function(object, alpha = NULL, ...) {
    if (is.null(alpha)) {
        # Return named list for all alpha values
        setNames(object$coefficients, paste0("alpha=", object$alpha))
    } else {
        # Find matching alpha and return coefficients
        idx <- match(alpha, object$alpha)
        if (any(is.na(idx))) {
            stop("alpha value(s) not found: ", paste(alpha[is.na(idx)], collapse = ", "))
        }
        if (length(idx) == 1) {
            object$coefficients[[idx]]
        } else {
            setNames(object$coefficients[idx], paste0("alpha=", alpha))
        }
    }
}

#' @export
vcov.SensIAT_marginal_mean_model <- function(object, alpha = NULL, ...) {
    if (is.null(alpha)) {
        setNames(object$coefficient.variance, paste0("alpha=", object$alpha))
    } else {
        idx <- match(alpha, object$alpha)
        if (any(is.na(idx))) {
            stop("alpha value(s) not found: ", paste(alpha[is.na(idx)], collapse = ", "))
        }
        if (length(idx) == 1) {
            object$coefficient.variance[[idx]]
        } else {
            setNames(object$coefficient.variance[idx], paste0("alpha=", alpha))
        }
    }
}

#' @export
print.SensIAT_marginal_mean_model <- function(x, digits = max(3L, getOption("digits") - 3L),
                                               markdown = isTRUE(getOption("knitr.in.progress")), ...) {
    link <- attr(x, "link") %||% "identity"
    loss <- attr(x, "loss") %||% "lp_mse"
    n_alpha <- length(x$alpha)
    n_coef <- length(x$coefficients[[1]])

    if (markdown) {
        cat("\n### SensIAT Marginal Mean Model\n\n")
        cat("| Property | Value |\n")
        cat("|:---------|:------|\n")
        cat("| Link | ", link, " |\n", sep = "")
        cat("| Loss | ", loss, " |\n", sep = "")
        cat("| Alpha values | ", n_alpha, " (", paste(x$alpha, collapse = ", "), ") |\n", sep = "")
        cat("| Spline coefficients | ", n_coef, " |\n\n", sep = "")
    } else {
        cat("\nSensIAT Marginal Mean Model\n\n")
        cat("Link:", link, "\n")
        cat("Loss:", loss, "\n")
        cat("Alpha values:", n_alpha, "(", paste(x$alpha, collapse = ", "), ")\n")
        cat("Spline coefficients:", n_coef, "\n\n")
    }
    invisible(x)
}

#' @export
summary.SensIAT_marginal_mean_model <- function(object, ...) {
    # Compute summary statistics for each alpha
    coef_summary <- purrr::map2(object$coefficients, object$coefficient.variance, function(cf, vr) {
        se <- sqrt(diag(vr))
        data.frame(
            Estimate = cf,
            Std.Error = se,
            row.names = paste0("B", seq_along(cf))
        )
    })
    names(coef_summary) <- paste0("alpha=", object$alpha)

    ans <- list(
        call = attr(object, "call"),
        alpha = object$alpha,
        link = attr(object, "link") %||% "identity",
        loss = attr(object, "loss") %||% "lp_mse",
        n_obs = nrow(object$data),
        coefficients = coef_summary
    )
    class(ans) <- "summary.SensIAT_marginal_mean_model"
    ans
}

#' @export
print.summary.SensIAT_marginal_mean_model <- function(x, digits = max(3L, getOption("digits") - 3L),
                                                       markdown = isTRUE(getOption("knitr.in.progress")), ...) {
    if (markdown) {
        cat("\n### SensIAT Marginal Mean Model Summary\n\n")
        cat("| Property | Value |\n")
        cat("|:---------|:------|\n")
        cat("| Link | ", x$link, " |\n", sep = "")
        cat("| Loss | ", x$loss, " |\n", sep = "")
        cat("| Observations | ", x$n_obs, " |\n\n", sep = "")

        for (nm in names(x$coefficients)) {
            cat("**Coefficients (", nm, "):**\n\n", sep = "")
            cf <- x$coefficients[[nm]]
            cat("| Term | Estimate | Std.Error |\n")
            cat("|:-----|--------:|----------:|\n")
            for (i in seq_len(nrow(cf))) {
                cat("| ", rownames(cf)[i], " | ",
                    format(cf$Estimate[i], digits = digits), " | ",
                    format(cf$Std.Error[i], digits = digits), " |\n", sep = "")
            }
            cat("\n")
        }
    } else {
        cat("\nSensIAT Marginal Mean Model Summary\n")
        cat(rep("=", 40), "\n", sep = "")
        cat("\nLink:", x$link, "\n")
        cat("Loss:", x$loss, "\n")
        cat("Observations:", x$n_obs, "\n\n")

        for (nm in names(x$coefficients)) {
            cat("Coefficients (", nm, "):\n", sep = "")
            print(x$coefficients[[nm]], digits = digits)
            cat("\n")
        }
    }
    invisible(x)
}

# Null-coalescing operator if not available
`%||%` <- function(x, y) if (is.null(x)) y else x
