#' @describeIn fit_SensIAT_within_group_model
#' Fit the Marginal Mean Sensitivity Analysis with Tilting Assumption for Both Treatment and Control Groups.
#'
#' @param data the full data set.
#' @param trt  an expression that determine what is treated as the treatment.
#'              Everything not treatment is considered control.
#' @param ... common arguments passed to `fit_SensIAT_within_group_model`.
#'
#' @return a list with class `SensIAT-fulldata-fitted-model` with two components,
#'      `control` and `treatment`, each of which is an independently fitted
#'      `SensIAT-within-group-fitted-model` fit with the fit_within_group_model
#'      function.
#' @export
fit_SensIAT_fulldata_model <- function(data, trt, ...) {
    trt <- rlang::enquo(trt) #< diffuses evaluation for tidy expressions

    structure(list(
        control   = fit_SensIAT_within_group_model(filter(data, !as.logical(!!trt)), ...),
        treatment = fit_SensIAT_within_group_model(filter(data, as.logical(!!trt)), ...)
    ), class = "SensIAT_fulldata_model")
}

#' @describeIn predict.SensIAT_within_group_model
#' For each combination of `time` and `alpha` estimate the mean response and
#' variance for each group as well as estimate the mean treatment effect and
#' variance.
#' @export
`predict.SensIAT_fulldata_model` <-
    function(object, time, ...) {
        control.predicted <- predict(object$control, time, ...)
        treatment.predicted <- predict(object$treatment, time, ...)

        full_join(
            control.predicted,
            treatment.predicted,
            by = c("time"),
            suffix = c("_control", "_treatment")
        ) |>
            mutate(
                mean_effect = mean_treatment - mean_control,
                var_effect = var_treatment + var_control
            )
    }
globalVariables(c("mean_effect", "mean_treatment", "mean_control", "var_effect", "var_treatment", "var_control"))

# Utility methods for SensIAT_fulldata_model ----------------------------------

#' @export
coef.SensIAT_fulldata_model <- function(object, alpha = NULL, ...) {
    list(
        control = coef(object$control, alpha = alpha, ...),
        treatment = coef(object$treatment, alpha = alpha, ...)
    )
}

#' @export
vcov.SensIAT_fulldata_model <- function(object, alpha = NULL, ...) {
    list(
        control = vcov(object$control, alpha = alpha, ...),
        treatment = vcov(object$treatment, alpha = alpha, ...)
    )
}

#' @export
print.SensIAT_fulldata_model <- function(x, digits = max(3L, getOption("digits") - 3L),
                                          markdown = isTRUE(getOption("knitr.in.progress")), ...) {
    n_alpha <- length(x$control$alpha)
    n_coef <- length(x$control$coefficients[[1]])
    n_control <- length(unique(x$control$data[[as.character(x$control$variables$id)]]))
    n_treatment <- length(unique(x$treatment$data[[as.character(x$treatment$variables$id)]]))

    if (markdown) {
        cat("\n### SensIAT Full Data Model (Treatment vs Control)\n\n")
        cat("| Property | Control | Treatment |\n")
        cat("|:---------|--------:|----------:|\n")
        cat("| Subjects | ", n_control, " | ", n_treatment, " |\n", sep = "")
        cat("| End time | ", x$control$End, " | ", x$treatment$End, " |\n\n", sep = "")
        cat("| Property | Value |\n")
        cat("|:---------|:------|\n")
        cat("| Alpha values | ", n_alpha, " (", paste(x$control$alpha, collapse = ", "), ") |\n", sep = "")
        cat("| Spline coefficients | ", n_coef, " |\n\n", sep = "")
    } else {
        cat("\nSensIAT Full Data Model (Treatment vs Control)\n\n")
        cat("Control subjects:", n_control, "\n")
        cat("Treatment subjects:", n_treatment, "\n")
        cat("Alpha values:", n_alpha, "(", paste(x$control$alpha, collapse = ", "), ")\n")
        cat("Spline coefficients:", n_coef, "\n")
        cat("End time:", x$control$End, "\n\n")
    }
    invisible(x)
}

#' @export
summary.SensIAT_fulldata_model <- function(object, ...) {
    ans <- list(
        control = summary(object$control, ...),
        treatment = summary(object$treatment, ...)
    )
    class(ans) <- "summary.SensIAT_fulldata_model"
    ans
}

#' @export
print.summary.SensIAT_fulldata_model <- function(x, digits = max(3L, getOption("digits") - 3L),
                                                  markdown = isTRUE(getOption("knitr.in.progress")), ...) {
    if (markdown) {
        cat("\n### SensIAT Full Data Model Summary\n\n")
        cat("#### Control Group\n\n")
        cat("| Property | Value |\n")
        cat("|:---------|:------|\n")
        cat("| Subjects | ", x$control$n_subjects, " |\n", sep = "")
        cat("| Observations | ", x$control$n_obs, " |\n\n", sep = "")

        for (nm in names(x$control$coefficients)) {
            cat("**Control Coefficients (", nm, "):**\n\n", sep = "")
            cf <- x$control$coefficients[[nm]]
            cat("| Term | Estimate | Std.Error |\n")
            cat("|:-----|--------:|----------:|\n")
            for (i in seq_len(nrow(cf))) {
                cat("| ", rownames(cf)[i], " | ",
                    format(cf$Estimate[i], digits = digits), " | ",
                    format(cf$Std.Error[i], digits = digits), " |\n", sep = "")
            }
            cat("\n")
        }

        cat("#### Treatment Group\n\n")
        cat("| Property | Value |\n")
        cat("|:---------|:------|\n")
        cat("| Subjects | ", x$treatment$n_subjects, " |\n", sep = "")
        cat("| Observations | ", x$treatment$n_obs, " |\n\n", sep = "")

        for (nm in names(x$treatment$coefficients)) {
            cat("**Treatment Coefficients (", nm, "):**\n\n", sep = "")
            cf <- x$treatment$coefficients[[nm]]
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
        cat("\nSensIAT Full Data Model Summary\n")
        cat(rep("=", 40), "\n", sep = "")

        cat("\n--- Control Group ---\n")
        cat("Subjects:", x$control$n_subjects, "\n")
        cat("Observations:", x$control$n_obs, "\n\n")
        for (nm in names(x$control$coefficients)) {
            cat("Control Coefficients (", nm, "):\n", sep = "")
            print(x$control$coefficients[[nm]], digits = digits)
            cat("\n")
        }

        cat("\n--- Treatment Group ---\n")
        cat("Subjects:", x$treatment$n_subjects, "\n")
        cat("Observations:", x$treatment$n_obs, "\n\n")
        for (nm in names(x$treatment$coefficients)) {
            cat("Treatment Coefficients (", nm, "):\n", sep = "")
            print(x$treatment$coefficients[[nm]], digits = digits)
            cat("\n")
        }
    }
    invisible(x)
}
