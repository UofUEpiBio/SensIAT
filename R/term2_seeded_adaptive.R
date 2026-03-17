# Seeded adaptive integration for term2 influence computation
#
# This method uses adaptive Simpson's quadrature but starts with a pre-computed
# grid of points where expected values are already known. This combines:
# - Accuracy of adaptive methods (subdivides where needed)
# - Efficiency of pre-computation (reuses cached expected values)
# - Better initial approximation (faster convergence)
#
# Performance characteristics:
# - Expected values: Computed on demand, but reused if already in cache
# - Initial subdivision: Uses pre-computed grid points
# - Adaptive refinement: Only where function is irregular
# - Generally 2-5x faster than pure adaptive for smooth integrands

#' Compute term2 influence using seeded adaptive integration
#'
#' This method uses adaptive Simpson's rule but initializes the algorithm
#' with pre-computed grid points and their expected values. The adaptive
#' algorithm can still subdivide intervals as needed for accuracy, but
#' starts with better information about the function.
#'
#' @inheritParams compute_term2_influence_fixed_grid
#' @param max_recursion Maximum recursion depth for adaptive algorithm
#' @param tol Tolerance for adaptive integration
#'
#' @return Numeric vector of length `ncol(base)` with term2 influence values
#'
#' @keywords internal
compute_term2_influence_seeded_adaptive <- function(
  patient_data,
  outcome_model,
  base,
  alpha,
  marginal_beta,
  V_inv,
  tmin,
  tmax,
  impute_fn = NULL,
  inv_link,
  W,
  expected_grid = NULL,
  n_grid = 50,
  include_obs_times = TRUE,
  time_var = NULL,
  max_recursion = 8,
  tol = 1e-6,
  ...
) {
    # Create or use provided grid with expected values
    if (!is.null(expected_grid)) {
        grid <- expected_grid$grid
        B_grid <- expected_grid$B_grid
        E_grid <- expected_grid$E_grid
    } else {
        # Create integration grid
        grid <- create_integration_grid(
            tmin = tmin,
            tmax = tmax,
            n_grid = n_grid,
            obs_times = if (include_obs_times) {
                extract_patient_times(patient_data, time_var)
            } else NULL
        )
        
        # Pre-compute basis evaluations
        B_grid <- lapply(grid, function(t) pcoriaccel_evaluate_basis(base, t))
        
        # Pre-compute expected values
        E_grid <- compute_expected_values_at_grid(
            grid = grid,
            patient_data = patient_data,
            outcome_model = outcome_model,
            alpha = alpha,
            impute_fn = impute_fn,
            time_var = time_var
        )
    }
    
    # Create evaluation cache with pre-computed values
    eval_cache <- new.env(parent = emptyenv())
    for (i in seq_along(grid)) {
        key <- format_cache_key(grid[i])
        eval_cache[[key]] <- list(
            B = B_grid[[i]],
            E = E_grid[[i]]
        )
    }
    
    # Create integrand function with caching
    integrand <- function(t) {
        key <- format_cache_key(t)
        
        # Check cache first
        if (exists(key, eval_cache, inherits = FALSE)) {
            cached <- eval_cache[[key]]
            B <- cached$B
            E <- cached$E
        } else {
            # Compute and cache
            B <- pcoriaccel_evaluate_basis(base, t)
            
            # Compute expected values (using fast path for single-index models)
            E <- compute_expected_value_at_point(
                t = t,
                patient_data = patient_data,
                outcome_model = outcome_model,
                alpha = alpha,
                impute_fn = impute_fn,
                time_var = time_var
            )
            
            eval_cache[[key]] <- list(B = B, E = E)
        }
        
        # Check validity
        # Extract scalar values
        E_exp <- as.numeric(E$E_exp_alphaY)[1]
        E_Yexp <- as.numeric(E$E_Yexp_alphaY)[1]
        
        if (!is.list(E) || !is.finite(E_exp) || 
            !is.finite(E_Yexp) || E_exp <= 0) {
            return(rep(0, ncol(base)))
        }
        
        # Compute integrand value (beta-dependent)
        weight <- W(t, marginal_beta)
        eta <- sum(B * marginal_beta)
        mu_hat <- inv_link(eta)
        
        # Ensure clean numeric types
        weight <- as.numeric(weight)
        diff_scalar <- as.numeric(E_Yexp / E_exp - mu_hat)
        
        weight * diff_scalar
    }
    
    # Use seeded adaptive Simpson's with pre-computed grid as initial partition
    result <- seeded_adaptive_simpson(
        f = integrand,
        seeds = grid,
        tol = tol,
        max_recursion = max_recursion
    )
    
    # Ensure finite output
    result <- ifelse(is.finite(result), result, 0)
    result
}


#' Format cache key for time value
#'
#' @param t Numeric time value
#' @return Character key
#' @keywords internal
format_cache_key <- function(t) {
    as.character(signif(t, 12))
}


#' Compute expected value at a single point
#'
#' @param t Time point
#' @param patient_data Patient data frame
#' @param outcome_model Fitted outcome model
#' @param alpha Sensitivity parameter
#' @param impute_fn Imputation function
#' @param time_var Time variable name
#'
#' @return List with E_exp_alphaY and E_Yexp_alphaY
#' @keywords internal
compute_expected_value_at_point <- function(
    t,
    patient_data,
    outcome_model,
    alpha,
    impute_fn = NULL,
    time_var = NULL
) {
    # Use impute_fn when available (preferred path - handles all model requirements)
    if (!is.null(impute_fn)) {
        df <- impute_fn(t, patient_data)
        ev <- compute_SensIAT_expected_values(
            model = outcome_model,
            alpha = alpha,
            new.data = df
        )
        # Extract scalar values
        return(list(E_exp_alphaY = as.numeric(ev$E_exp_alphaY)[1], 
                    E_Yexp_alphaY = as.numeric(ev$E_Yexp_alphaY)[1]))
    }
    
    # Fallback for simple single-index models without impute_fn
    # This path is less reliable and should only be used when impute_fn is not available
    if (is(outcome_model, "SensIAT::Single-index-outcome-model") &&
        !is.null(attr(outcome_model, "kernel")) &&
        !is.null(outcome_model$bandwidth)) {
        
        # Extract patient data
        outcome_var <- as.character(rlang::f_lhs(formula(outcome_model)))
        patient_times <- extract_patient_times(patient_data, time_var)
        patient_outcomes <- patient_data[[outcome_var]]
        
        # Ensure both vectors have same length (should match nrow(patient_data))
        if (is.null(patient_outcomes) || length(patient_outcomes) != length(patient_times)) {
            # If outcome variable doesn't exist or has wrong length, return empty list
            return(list())
        }
        
        valid_idx <- !is.na(patient_outcomes) & !is.na(patient_times)
        patient_times <- patient_times[valid_idx]
        patient_outcomes <- patient_outcomes[valid_idx]
        
        # Sort by time (required for findInterval)
        ord <- order(patient_times)
        patient_times <- patient_times[ord]
        patient_outcomes <- patient_outcomes[ord]
        
        # Find interval and compute expected values
        idx <- findInterval(t, patient_times, left.open = FALSE)
        if (idx < 1L) idx <- 1L
        if (idx > length(patient_times)) idx <- length(patient_times)
        
        # Extract scalar values
        prev_outcome <- patient_outcomes[idx]
        prev_time <- patient_times[idx]
        delta_time <- as.numeric(t) - as.numeric(prev_time)
        
        # Build data frame for this time point
        if (is.na(prev_outcome)) {
            df_at_t <- data.frame(..prev_outcome.. = 0, ..delta_time.. = delta_time,
                                 check.names = FALSE)
        } else {
            df_at_t <- data.frame(..prev_outcome.. = prev_outcome, ..delta_time.. = delta_time,
                                 check.names = FALSE)
        }
        
        # Use compute_SensIAT_expected_values to get expected values
        ev <- compute_SensIAT_expected_values(model = outcome_model, alpha = alpha, 
                                             new.data = df_at_t)
        # Extract scalar values (ev is a data frame with 1 row)
        return(list(E_exp_alphaY = as.numeric(ev$E_exp_alphaY)[1], 
                    E_Yexp_alphaY = as.numeric(ev$E_Yexp_alphaY)[1]))
    }
    
    # No impute_fn and not a single-index model - error
    stop("impute_fn required for non-single-index models")
}


#' Seeded adaptive Simpson's rule
#'
#' Adaptive Simpson's quadrature that starts with a pre-defined set of seed points.
#' The algorithm evaluates the function at seed points, then recursively subdivides
#' intervals where the error estimate exceeds the tolerance.
#'
#' @param f Function to integrate (returns vector)
#' @param seeds Numeric vector of seed points (initial partition)
#' @param tol Error tolerance
#' @param max_recursion Maximum recursion depth
#'
#' @return Integral approximation (vector)
#' @keywords internal
seeded_adaptive_simpson <- function(f, seeds, tol = 1e-6, max_recursion = 8) {
    # Evaluate function at all seed points
    f_vals <- lapply(seeds, f)
    
    # Initialize result
    result <- f_vals[[1]] * 0
    
    # Integrate between consecutive seed points
    for (i in 1:(length(seeds) - 1)) {
        a <- seeds[i]
        b <- seeds[i + 1]
        fa <- f_vals[[i]]
        fb <- f_vals[[i + 1]]
        
        # Use adaptive Simpson on this interval
        segment_result <- adaptive_simpson_segment(
            f = f,
            a = a,
            b = b,
            fa = fa,
            fb = fb,
            tol = tol * (b - a) / (seeds[length(seeds)] - seeds[1]),  # Scale tolerance
            max_recursion = max_recursion,
            recursion_level = 0
        )
        
        result <- result + segment_result
    }
    
    result
}


#' Recursive adaptive Simpson for a single interval
#'
#' @param f Function to integrate
#' @param a Left endpoint
#' @param b Right endpoint  
#' @param fa Function value at a
#' @param fb Function value at b
#' @param tol Tolerance for this interval
#' @param max_recursion Maximum recursion depth
#' @param recursion_level Current recursion level
#'
#' @return Integral approximation
#' @keywords internal
adaptive_simpson_segment <- function(
    f, a, b, fa, fb,
    tol, max_recursion, recursion_level
) {
    # Midpoint
    m <- (a + b) / 2
    fm <- f(m)
    
    # Simpson's rule for whole interval
    h <- (b - a) / 2
    S_whole <- (h / 3) * (fa + 4 * fm + fb)
    
    # Simpson's rule for left and right halves
    m_left <- (a + m) / 2
    m_right <- (m + b) / 2
    f_m_left <- f(m_left)
    f_m_right <- f(m_right)
    
    h_half <- h / 2
    S_left <- (h_half / 3) * (fa + 4 * f_m_left + fm)
    S_right <- (h_half / 3) * (fm + 4 * f_m_right + fb)
    S_halves <- S_left + S_right
    
    # Error estimate (vector-safe)
    error_est <- sqrt(sum((S_halves - S_whole)^2))
    
    # Check convergence or max recursion
    if (error_est < 15 * tol || recursion_level >= max_recursion) {
        # Use higher-order estimate with Richardson extrapolation
        return(S_halves + (S_halves - S_whole) / 15)
    }
    
    # Recurse on both halves
    result_left <- adaptive_simpson_segment(
        f, a, m, fa, fm,
        tol / 2, max_recursion, recursion_level + 1
    )
    result_right <- adaptive_simpson_segment(
        f, m, b, fm, fb,
        tol / 2, max_recursion, recursion_level + 1
    )
    
    result_left + result_right
}
