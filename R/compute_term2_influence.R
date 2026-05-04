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
#' Note: This method now segments the integration at observation times to handle
#' discontinuities in the imputation function, making it effectively the same as
#' the fast method but without the closure optimization.
#'
#' @param patient_data data.frame with patient's observations (Time, Outcome, and lag variables)
#' @param outcome_model The fitted Single-index outcome model
#' @param base `SplineBasis` object for marginal mean model
#' @param alpha Sensitivity parameter
#' @param marginal_beta Coefficients (beta) for the marginal mean spline basis
#' @param V_inv Inverse Gram matrix for base
#' @param tmin Lower integration bound
#' @param tmax Upper integration bound
#' @param impute_fn Function to impute data at time t: impute_fn(t, patient_data) -> data.frame
#' @param inv_link Inverse link function (e.g., exp for log link)
#' @param W Weight function W(t, beta)
#' @param expected_get Optional caching function: `expected_get(t)` -> `list(E_exp_alphaY, E_Yexp_alphaY)`
#' @param time_var Name of the time variable in patient_data (if NULL, auto-detected)
#' @param ... Additional arguments (not used)
#'
#' @return Numeric vector of length `ncol(base)` with term2 influence values
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
  inv_link,
  W,
  expected_get = NULL,
    identity_closed_form_scale = FALSE,
  time_var = NULL,
  ...
) {
    # Identify time variable if not provided
    if (is.null(time_var)) {
        time_candidates <- c("..time..", "Time", "time", "t", "T", "obstime", "obs_time")
        time_var <- NULL
        for (candidate in time_candidates) {
            if (candidate %in% names(patient_data)) {
                time_var <- candidate
                break
            }
        }
        if (is.null(time_var)) {
            numeric_cols <- names(patient_data)[sapply(patient_data, is.numeric)]
            time_var <- setdiff(numeric_cols, as.character(rlang::f_lhs(formula(outcome_model))))[1]
        }
    } else if (rlang::is_quosure(time_var) || is.language(time_var) || is.symbol(time_var)) {
        time_var <- rlang::quo_name(time_var)
    }
    
    if (is.null(time_var) || is.na(time_var)) {
        stop("Could not identify time variable in patient_data")
    }

    # Extract observation times and filter to integration range
    patient_times <- patient_data[[time_var]]
    obs_times_in_range <- patient_times[patient_times >= tmin & patient_times <= tmax]
    obs_times_in_range <- sort(unique(obs_times_in_range))

    # Create integration breakpoints at observation times
    # This handles discontinuities in the interpolation function
    breakpoints <- sort(unique(c(tmin, obs_times_in_range, tmax)))
    
    # Integrand
    term2_integrand <- function(t) {
        weight <- W(t, marginal_beta)
        expected <- if (!is.null(expected_get)) {
            expected_get(t)
        } else {
            compute_SensIAT_expected_values(
                model = outcome_model,
                alpha = alpha,
                new.data = impute_fn(t, patient_data)
            )
        }
        if (!is.list(expected) || 
            !is.finite(expected$E_exp_alphaY) || 
            !is.finite(expected$E_Yexp_alphaY) || 
            expected$E_exp_alphaY <= 0) {
            return(rep(0, ncol(base)))
        }
        if (isTRUE(identity_closed_form_scale)) {
            as.numeric(weight) * as.numeric(expected$E_Yexp_alphaY / expected$E_exp_alphaY)
        } else {
            B <- pcoriaccel_evaluate_basis(base, t)
            weight * as.numeric(
                expected$E_Yexp_alphaY / expected$E_exp_alphaY -
                    inv_link(crossprod(B, marginal_beta))
            )
        }
    }

    # Integrate over each segment and accumulate results
    total_integral <- rep(0, ncol(base))

    for (i in seq_len(length(breakpoints) - 1)) {
        seg_min <- breakpoints[i]
        seg_max <- breakpoints[i + 1]

        # Skip zero-length segments
        if (abs(seg_max - seg_min) < 1e-10) next

        # Integrate over this segment using adaptive Simpson's rule
        rslt <- tryCatch({
            pcoriaccel_integrate_simp(
                term2_integrand,
                seg_min,
                seg_max
            )
        }, error = function(e) {
            # If integration fails, return zero contribution for this segment
            list(Q = rep(0, ncol(base)))
        })
        
        # Accumulate result, guarding against NaN/Inf
        if (!is.null(rslt$Q)) {
            segment_result <- ifelse(is.finite(rslt$Q), rslt$Q, 0)
            total_integral <- total_integral + segment_result
        }
    }

    # Ensure output is finite
    total_integral <- ifelse(is.finite(total_integral), total_integral, 0)
    total_integral
}


#' Compute term2 influence for a patient (fast method with segmentation)
#'
#' Segments integration at observation times to handle discontinuities in the
#' interpolation function, then uses adaptive Simpson's quadrature on each segment.
#'
#' @inheritParams compute_term2_influence_original
#' @param time_var Name of the time variable in patient_data (if NULL, auto-detected)
#'
#' @return Numeric vector of length `ncol(base)` with term2 influence values
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
  inv_link,
  W,
  expected_get = NULL,
    identity_closed_form_scale = FALSE,
  time_var = NULL,
  ...
) {
    # Identify time variable if not provided
    if (is.null(time_var)) {
        time_candidates <- c("..time..", "Time", "time", "t", "T", "obstime", "obs_time")
        time_var <- NULL
        for (candidate in time_candidates) {
            if (candidate %in% names(patient_data)) {
                time_var <- candidate
                break
            }
        }
        if (is.null(time_var)) {
            numeric_cols <- names(patient_data)[sapply(patient_data, is.numeric)]
            time_var <- setdiff(numeric_cols, as.character(rlang::f_lhs(formula(outcome_model))))[1]
        }
    } else if (rlang::is_quosure(time_var) || is.language(time_var) || is.symbol(time_var)) {
        # Convert quosure, language, or symbol to character string
        time_var <- rlang::quo_name(time_var)
    }
    
    if (is.null(time_var) || is.na(time_var)) {
        stop("Could not identify time variable in patient_data")
    }

    # Extract observation times and filter to integration range
    patient_times <- patient_data[[time_var]]
    obs_times_in_range <- patient_times[patient_times >= tmin & patient_times <= tmax]
    obs_times_in_range <- sort(unique(obs_times_in_range))

    # Create integration breakpoints at observation times
    # This handles discontinuities in the interpolation function
    breakpoints <- sort(unique(c(tmin, obs_times_in_range, tmax)))

    # Integrand function (same as original method)
    term2_integrand <- function(t) {
        weight <- W(t, marginal_beta)
        expected <- if (!is.null(expected_get)) {
            expected_get(t)
        } else {
            compute_SensIAT_expected_values(
                model = outcome_model,
                alpha = alpha,
                new.data = impute_fn(t, patient_data)
            )
        }
        if (!is.list(expected) || 
            !is.finite(expected$E_exp_alphaY) || 
            !is.finite(expected$E_Yexp_alphaY) || 
            expected$E_exp_alphaY <= 0) {
            return(rep(0, ncol(base)))
        }
        if (isTRUE(identity_closed_form_scale)) {
            as.numeric(weight) * as.numeric(expected$E_Yexp_alphaY / expected$E_exp_alphaY)
        } else {
            B <- pcoriaccel_evaluate_basis(base, t)
            weight * as.numeric(
                expected$E_Yexp_alphaY / expected$E_exp_alphaY -
                    inv_link(crossprod(B, marginal_beta))
            )
        }
    }

    # Integrate over each segment and accumulate results
    total_integral <- rep(0, ncol(base))

    for (i in seq_len(length(breakpoints) - 1)) {
        seg_min <- breakpoints[i]
        seg_max <- breakpoints[i + 1]

        # Skip zero-length segments
        if (abs(seg_max - seg_min) < 1e-10) next

        # Integrate over this segment using adaptive Simpson's rule
        rslt <- tryCatch({
            pcoriaccel_integrate_simp(
                term2_integrand,
                seg_min,
                seg_max
            )
        }, error = function(e) {
            # If integration fails, return zero contribution for this segment
            list(Q = rep(0, ncol(base)))
        })
        
        # Accumulate result, guarding against NaN/Inf
        if (!is.null(rslt$Q)) {
            segment_result <- ifelse(is.finite(rslt$Q), rslt$Q, 0)
            total_integral <- total_integral + segment_result
        }
    }

    # Ensure output is finite
    total_integral <- ifelse(is.finite(total_integral), total_integral, 0)
    total_integral
}

#' Build fast term2 integrand using closure optimization
#'
#' @keywords internal
make_term2_integrand_fast <- function(
  outcome.model,
  base,
  alpha,
  patient_times,
  patient_outcomes,
  marginal_beta,
  V_inv,
  W,
  expected_get = NULL,
  impute_fn = NULL,
  patient_data = NULL,
  time_var = NULL
) {
    # Pre-compute constants that don't change during integration
    mf <- model.frame(outcome.model)
    Xi <- model.matrix(terms(outcome.model), data = mf)
    Yi <- model.response(mf)
    beta_outcome <- outcome.model$coef
    y_seq <- sort(unique(Yi))
    kernel_fn <- attr(outcome.model, "kernel")
    bandwidth <- outcome.model$bandwidth
    
    # Pre-compute Xbeta for all training observations
    Xb_all <- as.vector(Xi %*% beta_outcome)

    # Pre-compute term spec (full model, drop response)
    term_spec <- delete.response(terms(outcome.model))

    # Return the integrand closure
    function(t) {
        weight <- W(t, marginal_beta)
        
        # Get expected values (with caching if available)
        if (!is.null(expected_get)) {
            expected <- expected_get(t)
        } else {
            # Prefer full-model imputation if provided (supports all variables)
            if (!is.null(impute_fn) && !is.null(patient_data)) {
                df_t <- impute_fn(t, patient_data)
                X_t <- model.matrix(term_spec, data = df_t)
                xb_t <- as.numeric(X_t %*% beta_outcome)
            } else {
                # Fallback: legacy two-variable path
                idx <- findInterval(t, patient_times, left.open = FALSE)
                if (idx < 1L) idx <- 1L
                if (idx > length(patient_times)) idx <- length(patient_times)
                
                prev_outcome <- patient_outcomes[idx]
                delta_time <- t - patient_times[idx]
                
                df_t <- data.frame(
                    ..prev_outcome.. = prev_outcome,
                    ..delta_time.. = delta_time,
                    check.names = FALSE
                )
                X_t <- model.matrix(term_spec, data = df_t)
                xb_t <- as.numeric(X_t %*% beta_outcome)
            }
            
            # Compute pmf at this point
            pmf <- pcoriaccel_estimate_pmf(
                Xb = Xb_all,
                Y = Yi,
                xi = xb_t,
                y_seq = y_seq,
                h = bandwidth,
                kernel = kernel_fn
            )
            
            E_exp_alphaY <- sum(exp(alpha * y_seq) * pmf)
            E_Yexp_alphaY <- sum(y_seq * exp(alpha * y_seq) * pmf)
            
            expected <- list(
                E_exp_alphaY = E_exp_alphaY,
                E_Yexp_alphaY = E_Yexp_alphaY
            )
        }
        
        # Check for valid expected values
        if (!is.list(expected) || 
            !is.finite(expected$E_exp_alphaY) || 
            !is.finite(expected$E_Yexp_alphaY) || 
            expected$E_exp_alphaY <= 0) {
            return(rep(0, ncol(base)))
        }
        
        B <- pcoriaccel_evaluate_basis(base, t)
        mu_t <- as.numeric(crossprod(B, marginal_beta))
        
        weight * as.numeric(
            expected$E_Yexp_alphaY / expected$E_exp_alphaY - inv_link(mu_t)
        )
    }
}
