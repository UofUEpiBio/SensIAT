# Compute term2 influence for a single patient using different methods
#
# These functions compute the term 2 influence contribution for a single patient
# in the marginal mean model estimation. Both methods produce identical results
# but differ in performance characteristics.

#' Compute term2 influence for a patient (original method)
#'
#' @param patient_data data.frame with patient's observations (Time, Outcome, and lag variables)
#' @param outcome_model The fitted Single-index outcome model
#' @param base SplineBasis for marginal mean model
#' @param alpha Sensitivity parameter
#' @param marginal_beta Coefficients (beta) for the marginal mean spline basis
#' @param V_inv Inverse Gram matrix for base
#' @param tmin Lower integration bound
#' @param tmax Upper integration bound
#' @param impute_fn Function to impute data at time t: impute_fn(t, patient_data) -> data.frame
#' @param inv_link Inverse link function (e.g., exp for log link)
#'
#' @return Numeric vector of length ncol(base) with term2 influence values
#'
#' @keywords internal
compute_term2_influence_original <- function(
    patient_data,
    outcome_model,
    base,
    alpha,
    marginal_beta,
    V_inv,
    tmin,
    tmax,
    impute_fn,
    inv_link
) {
  # Weight function W(t, beta)
  W <- function(t, beta) {
    B <- as.vector(pcoriaccel_evaluate_basis(base, t))
    mu <- sum(B * beta)
    (V_inv %*% B) / mu
  }

  # Integrand
  term2_integrand <- function(t) {
    weight <- W(t, marginal_beta)
    expected <- compute_SensIAT_expected_values(
      model = outcome_model,
      alpha = alpha,
      new.data = impute_fn(t, patient_data)
    )
    B <- pcoriaccel_evaluate_basis(base, t)
    weight * as.numeric(
      expected$E_Yexp_alphaY / expected$E_exp_alphaY -
        inv_link(crossprod(B, marginal_beta))
    )
  }

  rslt <- pcoriaccel_integrate_simp(
    term2_integrand,
    tmin,
    tmax
  )
  rslt$Q
}


#' Compute term2 influence for a patient (fast method)
#'
#' @inheritParams compute_term2_influence_original
#'
#' @return Numeric vector of length ncol(base) with term2 influence values
#'
#' @keywords internal
compute_term2_influence_fast <- function(
    patient_data,
    outcome_model,
    base,
    alpha,
    marginal_beta,
    V_inv,
    tmin,
    tmax,
    impute_fn,
    inv_link
) {
  # Build fast integrand using closure-based optimization
  integrand_fast <- make_term2_integrand_fast(
    outcome.model = outcome_model,
    base = base,
    alpha = alpha,
    patient_times = patient_data$Time,
    patient_outcomes = patient_data$Outcome,
    marginal_beta = marginal_beta,
    V_inv = V_inv
  )

  rslt <- pcoriaccel_integrate_simp(
    integrand_fast,
    tmin,
    tmax
  )
  rslt$Q
}
