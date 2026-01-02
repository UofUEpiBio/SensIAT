#' Fit SensIAT Marginal Mean Model (Unified)
#'
#' This function fits a marginal mean model for sensitivity analysis, supporting both
#' single-index (linear) and generalized linear outcome models.
#'
#' @param data Data frame containing all required variables.
#' @param outcome.model Outcome model object or formula (GLM, single-index, or LM).
#' @param intensity.model Intensity model object (e.g., coxph).
#' @param id Subject identifier variable.
#' @param time Time variable.
#' @param alpha Sensitivity parameter(s).
#' @param knots Spline knot locations.
#' @param ... Additional arguments passed to underlying methods.
#' @return A fitted SensIAT marginal mean model object.
#' @export
fit_marginal_model <- function(
  data,
  outcome.model,
  intensity.model,
  id,
  time,
  alpha = 0,
  knots,
  ...
) {
  UseMethod("fit_marginal_model", outcome.model)
}
