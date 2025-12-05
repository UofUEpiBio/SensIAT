# Compute term2 influence for a single patient using different methods
#
# These functions compute the term 2 influence contribution for a single patient
# in the marginal mean model estimation. Both methods produce identical results
# but differ in performance characteristics and generality.
#
# ## Comparison with Simulation Code Approaches
# 
# Simulation code (e.g., inst/3. Simu_point_estimate_fast.R) uses a pre-computation
# strategy that is highly efficient for specific single-index model formulas but
# lacks the generality of this implementation:
# 
# **Simulation Approach Advantages:**
# - Pre-computes expected values once per alpha
# - 2-5x faster for repeated optimization iterations  
# - Uses trapezoidal rule integration over fixed grids
# 
# **Package Approach Advantages (this file):**
# - Supports any outcome model with compute_SensIAT_expected_values method
# - Uses adaptive Simpson's rule for better numerical accuracy
# - Handles varying observation times and data structures robustly
# - Modular design allows method dispatch based on model type
#
# ## Integration Methods
# 
# Both methods use adaptive Simpson's quadrature (`pcoriaccel_integrate_simp`)
# which is more accurate than fixed-step trapezoidal rule used in simulation code.
# 
# ## Numerical Improvements
# 
# Weight functions now use exp(-μ) multiplication instead of division for
# better numerical stability with large linear predictors.

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
    as.vector((V_inv %*% B) * exp(-mu))
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
  # Extract outcome variable name from the model formula
  outcome_var <- as.character(rlang::f_lhs(formula(outcome_model)))
  
  # Get times and outcomes from patient_data
  # The patient_data should have the outcome variable
  patient_outcomes <- patient_data[[outcome_var]]
  
  # For time, we need to find it. Common patterns: Time, time, t, T
  # Or we can look for the column that's not the outcome
  time_candidates <- c("Time", "time", "t", "T", "obstime", "obs_time")
  time_var <- NULL
  for (candidate in time_candidates) {
    if (candidate %in% names(patient_data)) {
      time_var <- candidate
      break
    }
  }
  
  # If still not found, take the first numeric column that's not the outcome
  if (is.null(time_var)) {
    numeric_cols <- names(patient_data)[sapply(patient_data, is.numeric)]
    time_var <- setdiff(numeric_cols, outcome_var)[1]
  }
  
  if (is.null(time_var) || is.na(time_var)) {
    stop("Could not identify time variable in patient_data")
  }
  
  patient_times <- patient_data[[time_var]]
  
  # Remove NA observations
  valid_idx <- !is.na(patient_outcomes)
  patient_times <- patient_times[valid_idx]
  patient_outcomes <- patient_outcomes[valid_idx]
  
  # Build fast integrand using closure-based optimization
  integrand_fast <- make_term2_integrand_fast(
    outcome.model = outcome_model,
    base = base,
    alpha = alpha,
    patient_times = patient_times,
    patient_outcomes = patient_outcomes,
    marginal_beta = marginal_beta,
    V_inv = V_inv
  )

  # Split integration into segments at observation times to handle discontinuities
  # Integration segments: [tmin, t1], [t1, t2], ..., [t_{n-1}, tmax]
  obs_times_in_range <- patient_times[patient_times >= tmin & patient_times <= tmax]
  
  # Create breakpoints for integration
  breakpoints <- sort(unique(c(tmin, obs_times_in_range, tmax)))
  
  # Integrate over each segment and sum
  total_integral <- rep(0, ncol(base))
  
  for (i in seq_len(length(breakpoints) - 1)) {
    seg_min <- breakpoints[i]
    seg_max <- breakpoints[i + 1]
    
    # Skip zero-length segments
    if (abs(seg_max - seg_min) < 1e-10) next
    
    rslt <- pcoriaccel_integrate_simp(
      integrand_fast,
      seg_min,
      seg_max
    )
    total_integral <- total_integral + rslt$Q
  }
  
  total_integral
}
