# Utilities for fast Term 2 integrand computation (reduced allocations)
#
# The main entry point is make_term2_integrand_fast(), which returns a closure
# integrand(t) with minimal overhead:
# - Avoids building data.frames per t
# - Avoids model.matrix per t by caching the spline basis for prev_outcome
# - Uses pcoriaccel_estimate_pmf directly with precomputed Xb for the outcome model
# - Uses findInterval() to locate the most recent observation efficiently
#
# Assumptions:
# - outcome.model is a Single-index-outcome-model with formula like
#   Outcome ~ ns(..prev_outcome.., df=3) + ..delta_time.. - 1
# - Only ..prev_outcome.. and ..delta_time.. vary per t for a fixed patient
#
# If the formula differs substantially, this utility falls back to a safe path
# via model.matrix() with a 1-row data.frame, but still avoids copying the full
# patient data each time.


#' Build a fast term2 integrand for a single patient (closure)
#'
#' @param outcome.model The fitted Single-index outcome model
#' @param base SplineBasis for marginal mean model
#' @param alpha Sensitivity parameter
#' @param patient_times Numeric vector of observation times for the patient (sorted asc)
#' @param patient_outcomes Numeric vector of outcomes aligned with patient_times
#' @param marginal_beta Coefficients (beta) for the marginal mean spline basis
#' @param V_inv Precomputed inverse Gram matrix for base (optional; computed if NULL)
#' @return A function f(t) computing the term2 integrand at scalar t
#' @export
make_term2_integrand_fast <- function(
  outcome.model,
  base,
  alpha,
  patient_times,
  patient_outcomes,
  marginal_beta,
  V_inv = NULL
){
  stopifnot(is.numeric(patient_times), is.numeric(patient_outcomes))
  stopifnot(length(patient_times) == length(patient_outcomes))

  # Precompute constants
  if (is.null(V_inv)) {
    V_inv <- GramMatrix(base)
    V_inv <- solve(V_inv)
  }
  kernel <- attr(outcome.model, "kernel")
  bandwidth <- outcome.model$bandwidth

  # Outcome model design (training data)
  Xi <- model.matrix(terms(outcome.model), outcome.model$data)
  Yi <- model.response(model.frame(outcome.model))
  beta_outcome <- outcome.model$coef

  # Precompute items for pmf estimation
  Xb_all <- as.vector(Xi %*% beta_outcome)
  y_seq <- sort(unique(Yi))

  # Identify which column corresponds to ..delta_time.. in model.matrix
  # We'll construct a template row with delta_time=0 for each unique prev_outcome
  # IMPORTANT: drop the response to avoid needing 'Outcome' in single-row temp frames
  term_spec <- delete.response(terms(outcome.model))
  # Build a 1-row data.frame template with the necessary variables
  template_df <- data.frame(
    ..prev_outcome.. = 0,
    ..delta_time.. = 0,
    check.names = FALSE
  )
  template_row <- model.matrix(term_spec, data = template_df)
  mm_colnames <- colnames(template_row)
  # Best-effort lookup of delta_time column
  idx_delta <- which(mm_colnames == "..delta_time..")
  if (length(idx_delta) == 0L) {
    # Fallback: search by pattern
    idx_delta <- grep("delta_time", mm_colnames, fixed = TRUE)
  }
  if (length(idx_delta) != 1L) {
    # If still ambiguous, we will rebuild the entire row via model.matrix later
    idx_delta <- NA_integer_
  }

  # Cache ns-basis row for each unique observed prev_outcome value
  unique_prev_outcomes <- unique(patient_outcomes)
  ns_cache <- new.env(parent = emptyenv())
  for (val in unique_prev_outcomes) {
    df1 <- data.frame(..prev_outcome.. = val, ..delta_time.. = 0, check.names = FALSE)
    row0 <- model.matrix(term_spec, data = df1)
    ns_cache[[as.character(val)]] <- row0
  }

  # Fast helpers
  # Weight function W(t, beta) used by lp_mse/log link in current implementation
  W_fun <- function(t) {
    B <- as.vector(pcoriaccel_evaluate_basis(base, t))
    mu <- sum(B * marginal_beta) # linear predictor on link scale per current code
    (V_inv %*% B) / mu
  }

  # Build integrand closure (no data.frame allocations per call)
  function(t) {
    # Locate last observed index efficiently
    idx <- findInterval(t, patient_times, left.open = FALSE)
    if (idx < 1L) idx <- 1L
    if (idx > length(patient_times)) idx <- length(patient_times)

    prev_outcome <- patient_outcomes[idx]
    delta_time <- t - patient_times[idx]

    # Build outcome model row efficiently
    key <- as.character(prev_outcome)
    if (!is.null(ns_cache[[key]]) && !is.na(idx_delta)) {
      x_row <- ns_cache[[key]][1, , drop = TRUE]
      x_row[idx_delta] <- delta_time
    } else {
      # Fallback: construct a single-row model.matrix (still cheap)
      df1 <- data.frame(..prev_outcome.. = prev_outcome, ..delta_time.. = delta_time, check.names = FALSE)
      x_row <- model.matrix(term_spec, data = df1)[1, , drop = TRUE]
    }

    # Single-index xb for pmf estimation
    xb <- sum(x_row * beta_outcome)

    # Estimate pmf and expectations directly
    pmf <- pcoriaccel_estimate_pmf(
      Xb = Xb_all,
      Y = Yi,
      xi = xb,
      y_seq = y_seq,
      h = bandwidth,
      kernel = kernel
    )
    if (all(pmf == 0)) return(rep(0, ncol(base)))

    E_exp_alphaY <- sum(exp(alpha * y_seq) * pmf)
    E_Yexp_alphaY <- sum(y_seq * exp(alpha * y_seq) * pmf)

    # Marginal mean pieces
    B <- pcoriaccel_evaluate_basis(base, t)
    weight <- W_fun(t)
    eta <- sum(B * marginal_beta)
    mu_hat <- exp(eta) # inv.link for log link

    as.numeric(weight) * as.numeric(E_Yexp_alphaY / E_exp_alphaY - mu_hat)
  }
}
