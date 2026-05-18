# Helper: format a list of coefficient vectors as a matrix, indexed by alpha
.coef_matrix_by_alpha <- function(coef_list, alpha) {
    mat <- do.call(cbind, coef_list)
    n_coefs <- nrow(mat)
    rownames(mat) <- paste0("b", seq_len(n_coefs))
    colnames(mat) <- format(alpha, digits = 4, trim = TRUE)
    mat
}

# Helper: try to extract a formula from a model object, return NULL silently if not possible
.try_formula <- function(model) {
    tryCatch(formula(model), error = function(e) NULL)
}

#' Print method for `SensIAT_within_group_model`
#'
#' @param x A `SensIAT_within_group_model` object.
#' @param ... Currently ignored.
#'
#' @return `x`, invisibly.
#' @export
print.SensIAT_within_group_model <- function(x, ...) {
    cat("SensIAT Within-Group Model\n")

    cl <- attr(x, "call")
    if (!is.null(cl)) {
        cat("\nCall:\n  ", deparse(cl, width.cutoff = 70L), "\n", sep = "")
    }

    if (!is.null(x$link) && x$link != "identity") {
        cat("\nLink:", x$link, "\n")
    }

    intensity_formula <- .try_formula(x$models$intensity)
    if (!is.null(intensity_formula)) {
        cat("\nIntensity model:\n  ", deparse(intensity_formula), "\n", sep = "")
    }

    outcome_formula <- .try_formula(x$models$outcome)
    if (!is.null(outcome_formula)) {
        cat("\nOutcome model:\n  ", deparse(outcome_formula), "\n", sep = "")
    }

    alpha <- x$alpha
    coefs <- x$coefficients
    cat("\nCoefficients (", length(alpha), " alpha value(s), alpha by column):\n", sep = "")
    mat <- .coef_matrix_by_alpha(coefs, alpha)
    print(mat, digits = 4)

    invisible(x)
}

#' Print method for `SensIAT_fulldata_model`
#'
#' @param x A `SensIAT_fulldata_model` object.
#' @param ... Currently ignored.
#'
#' @return `x`, invisibly.
#' @export
print.SensIAT_fulldata_model <- function(x, ...) {
    cat("SensIAT Full-Data Model\n")
    cat("\n--- Control group ---\n")
    print(x$control, ...)
    cat("\n--- Treatment group ---\n")
    print(x$treatment, ...)
    invisible(x)
}

#' Print method for `SensIAT_marginal_outcome_model`
#'
#' @param x A `SensIAT_marginal_outcome_model` object.
#' @param ... Currently ignored.
#'
#' @return `x`, invisibly.
#' @export
print.SensIAT_marginal_outcome_model <- function(x, ...) {
    cat("SensIAT Marginal Outcome Model\n")

    cl <- attr(x, "call")
    if (!is.null(cl)) {
        cat("\nCall:\n  ", deparse(cl, width.cutoff = 70L), "\n", sep = "")
    }

    intensity_formula <- .try_formula(x$models$intensity)
    if (!is.null(intensity_formula)) {
        cat("\nIntensity model:\n  ", deparse(intensity_formula), "\n", sep = "")
    }

    outcome_formula <- .try_formula(x$models$outcome)
    if (!is.null(outcome_formula)) {
        cat("\nOutcome model:\n  ", deparse(outcome_formula), "\n", sep = "")
    }

    alpha <- x$alpha
    coefs <- x$coefficients
    cat("\nCoefficients (", length(alpha), " alpha value(s), alpha by column):\n", sep = "")
    mat <- .coef_matrix_by_alpha(coefs, alpha)
    print(mat, digits = 4)

    invisible(x)
}

#' Print method for `SensIAT_marginal_mean_model_generalized`
#'
#' @param x A `SensIAT_marginal_mean_model_generalized` object.
#' @param ... Currently ignored.
#'
#' @return `x`, invisibly.
#' @export
print.SensIAT_marginal_mean_model_generalized <- function(x, ...) {
    cat("SensIAT Marginal Mean Model (Generalized)\n")

    cl <- attr(x, "call")
    if (!is.null(cl)) {
        cat("\nCall:\n  ", deparse(cl, width.cutoff = 70L), "\n", sep = "")
    }

    intensity_formula <- .try_formula(x$models$intensity)
    if (!is.null(intensity_formula)) {
        cat("\nIntensity model:\n  ", deparse(intensity_formula), "\n", sep = "")
    }

    outcome_formula <- .try_formula(x$models$outcome)
    if (!is.null(outcome_formula)) {
        cat("\nOutcome model:\n  ", deparse(outcome_formula), "\n", sep = "")
    }

    alpha <- x$alpha
    coefs <- x$coefficients
    cat("\nCoefficients (", length(alpha), " alpha value(s), alpha by column):\n", sep = "")
    mat <- .coef_matrix_by_alpha(coefs, alpha)
    print(mat, digits = 4)

    invisible(x)
}

#' Print method for `SensIAT::Single-index-outcome-model`
#'
#' @param x A `SensIAT::Single-index-outcome-model` object.
#' @param ... Currently ignored.
#'
#' @return `x`, invisibly.
#' @export
`print.SensIAT::Single-index-outcome-model` <- function(x, ...) {
    cat("SensIAT Single-Index Outcome Model\n")

    trms <- attr(x, "terms")
    if (!is.null(trms)) {
        cat("\nFormula:\n  ", deparse(as.formula(trms)), "\n", sep = "")
    }

    kernel <- attr(x, "kernel")
    if (!is.null(kernel)) {
        cat("Kernel:", kernel, "\n")
    }

    if (!is.null(x$bandwidth)) {
        cat("Bandwidth:", x$bandwidth, "\n")
    }

    cat("\nSingle-index coefficients:\n")
    coef_vec <- x$coef
    names(coef_vec) <- paste0("b", seq_along(coef_vec))
    print(coef_vec, digits = 4)

    invisible(x)
}
