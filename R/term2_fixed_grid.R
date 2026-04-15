# Fixed-grid integration for term2 influence computation
#
# This method pre-computes expected values at a fixed grid of time points,
# then uses composite Simpson's rule or trapezoidal rule for integration.
# It is significantly faster than adaptive methods when:
# - Multiple beta iterations are needed (expected values computed once per alpha)
# - The integrand is relatively smooth between observation times
# - Grid density is sufficient for desired accuracy
#
# Performance characteristics:
# - Expected values: O(n_grid) per patient per alpha (computed once)
# - Integration per beta iteration: O(n_grid) (just weighted sum)
# - No adaptive subdivision overhead
# - Predictable, constant cost

#' Compute term2 influence using fixed-grid integration
#'
#' This method evaluates the integrand at a pre-determined grid of points
#' and uses composite Simpson's rule for integration. Expected values are
#' pre-computed at grid points and reused across beta iterations.
#'
#' @param patient_data data.frame with patient's observations
#' @param outcome_model The fitted outcome model
#' @param base `SplineBasis` object for marginal mean model
#' @param alpha Sensitivity parameter
#' @param marginal_beta Coefficients (beta) for the marginal mean spline basis
#' @param V_inv Inverse Gram matrix for base
#' @param tmin Lower integration bound
#' @param tmax Upper integration bound
#' @param impute_fn Function to impute data at time t (not used if expected_grid provided)
#' @param inv_link Inverse link function
#' @param W Weight function W(t, beta)
#' @param expected_grid Optional pre-computed expected values at grid points:
#'        list with elements:
#'        - grid: numeric vector of time points
#'        - B_grid: list of basis evaluations at grid points
#'        - E_grid: list of expected values list(E_exp_alphaY, E_Yexp_alphaY)
#' @param n_grid Number of grid points to use if expected_grid not provided
#' @param rule Integration rule: "simpson" (default) or "trapezoid"
#' @param include_obs_times Whether to include observation times in grid (default TRUE)
#' @param time_var Time variable name (optional)
#' @param ... Additional arguments (not used)
#'
#' @return Numeric vector of length `ncol(base)` with term2 influence values
#'
#' @keywords internal
compute_term2_influence_fixed_grid <- function(
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
  n_grid = 100,
  rule = c("simpson", "trapezoid"),
  include_obs_times = TRUE,
  time_var = NULL,
  ...
) {
    rule <- match.arg(rule)
    
    # If expected_grid is provided, use it; otherwise compute grid and expected values
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
        
        # Pre-compute basis evaluations at grid points
        B_grid <- lapply(grid, function(t) pcoriaccel_evaluate_basis(base, t))
        
        # Compute expected values at grid points
        E_grid <- compute_expected_values_at_grid(
            grid = grid,
            patient_data = patient_data,
            outcome_model = outcome_model,
            alpha = alpha,
            impute_fn = impute_fn,
            time_var = time_var
        )
    }
    
    # Compute integrand values at grid points
    # This is the only beta-dependent step
    integrand_values <- mapply(
        function(t, B, E) {
            # Extract scalar values
            E_exp <- as.numeric(E$E_exp_alphaY)[1]
            E_Yexp <- as.numeric(E$E_Yexp_alphaY)[1]
            
            if (!is.list(E) || !is.finite(E_exp) || 
                !is.finite(E_Yexp) || E_exp <= 0) {
                return(rep(0, ncol(base)))
            }
            
            weight <- W(t, marginal_beta)
            eta <- sum(B * marginal_beta)
            mu_hat <- inv_link(eta)
            
            # Ensure weight is a numeric vector (not matrix)
            weight <- as.numeric(weight)
            # Compute scalar difference
            diff_scalar <- as.numeric(E_Yexp / E_exp - mu_hat)
            
            # Element-wise multiplication
            weight * diff_scalar
        },
        grid,
        B_grid,
        E_grid,
        SIMPLIFY = FALSE
    )
    
    # Integrate using specified rule
    if (rule == "simpson") {
        result <- composite_simpson(integrand_values, grid)
    } else {
        result <- composite_trapezoid(integrand_values, grid)
    }
    
    # Ensure finite output
    result <- ifelse(is.finite(result), result, 0)
    result
}


#' Create integration grid
#'
#' @param tmin Lower bound
#' @param tmax Upper bound
#' @param n_grid Number of grid points
#' @param obs_times Optional vector of observation times to include
#'
#' @return Sorted numeric vector of grid points
#' @keywords internal
create_integration_grid <- function(tmin, tmax, n_grid = 100, obs_times = NULL) {
    # Start with uniform grid
    grid <- seq(tmin, tmax, length.out = n_grid)
    
    # Add observation times if provided
    if (!is.null(obs_times)) {
        obs_in_range <- obs_times[obs_times >= tmin & obs_times <= tmax]
        grid <- sort(unique(c(grid, obs_in_range)))
    }
    
    grid
}


#' Extract patient times from data
#'
#' @param patient_data Patient data frame
#' @param time_var Time variable name (optional, can be a quosure or string)
#'
#' @return Numeric vector of times
#' @keywords internal
extract_patient_times <- function(patient_data, time_var = NULL) {
    if (!is.null(time_var)) {
        # Handle both quosures and strings
        if (rlang::is_quosure(time_var)) {
            # Evaluate the quosure in the context of patient_data
            time_values <- rlang::eval_tidy(time_var, patient_data)
            return(time_values)
        } else if (is.character(time_var)) {
            # Already a character, use it directly
            if (!time_var %in% names(patient_data)) {
                stop("Time variable '", time_var, "' not found in patient data")
            }
            return(patient_data[[time_var]])
        } else {
            stop("time_var must be a quosure or character string")
        }
    }
    
    # Auto-detect time variable
    time_candidates <- c("..time..", "Time", "time", "t", "T", "obstime", "obs_time")
    for (candidate in time_candidates) {
        if (candidate %in% names(patient_data)) {
            return(patient_data[[candidate]])
        }
    }
    
    # Fallback: first numeric column
    numeric_cols <- names(patient_data)[sapply(patient_data, is.numeric)]
    if (length(numeric_cols) > 0) {
        return(patient_data[[numeric_cols[1]]])
    }
    
    numeric(0)
}


#' Compute expected values at grid points for a patient
#'
#' @param grid Numeric vector of time points
#' @param patient_data Patient data frame
#' @param outcome_model Fitted outcome model
#' @param alpha Sensitivity parameter
#' @param impute_fn Imputation function
#' @param time_var Time variable name (optional)
#'
#' @return List of expected value lists at each grid point
#' @keywords internal
compute_expected_values_at_grid <- function(
    grid,
    patient_data,
    outcome_model,
    alpha,
    impute_fn = NULL,
    time_var = NULL
) {
    # Use impute_fn when available (preferred path - handles all model requirements)
    if (!is.null(impute_fn)) {
        return(lapply(grid, function(t) {
            df <- impute_fn(t, patient_data)
            ev <- compute_SensIAT_expected_values(
                model = outcome_model,
                alpha = alpha,
                new.data = df
            )
            # Extract scalar values
            list(E_exp_alphaY = as.numeric(ev$E_exp_alphaY)[1], 
                 E_Yexp_alphaY = as.numeric(ev$E_Yexp_alphaY)[1])
        }))
    }
    
    # Fallback for simple single-index models without impute_fn
    # This path is less reliable and should only be used when impute_fn is not available
    if (is(outcome_model, "SensIAT::Single-index-outcome-model") &&
        !is.null(attr(outcome_model, "kernel")) &&
        !is.null(outcome_model$bandwidth)) {
        
        # Extract patient times and outcomes
        outcome_var <- as.character(rlang::f_lhs(formula(outcome_model)))
        patient_times <- extract_patient_times(patient_data, time_var)
        patient_outcomes <- patient_data[[outcome_var]]
        
        # Ensure both vectors have same length (should match nrow(patient_data))
        if (is.null(patient_outcomes) || length(patient_outcomes) != length(patient_times)) {
            # If outcome variable doesn't exist or has wrong length, return empty list
            return(list())
        }
        
        # Remove NA observations
        valid_idx <- !is.na(patient_outcomes) & !is.na(patient_times)
        patient_times <- patient_times[valid_idx]
        patient_outcomes <- patient_outcomes[valid_idx]
        
        # Sort by time (required for findInterval)
        ord <- order(patient_times)
        patient_times <- patient_times[ord]
        patient_outcomes <- patient_outcomes[ord]
        
        # Compute expected values at each grid point
        E_grid <- lapply(grid, function(t) {
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
            # This handles all the model.matrix complexity for us
            ev <- compute_SensIAT_expected_values(model = outcome_model, alpha = alpha, 
                                                 new.data = df_at_t)
            # Extract scalar values (ev is a data frame with 1 row)
            list(E_exp_alphaY = as.numeric(ev$E_exp_alphaY)[1], 
                 E_Yexp_alphaY = as.numeric(ev$E_Yexp_alphaY)[1])
        })
        
        return(E_grid)
    }
    
    # No impute_fn and not a single-index model - error
    stop("impute_fn required for non-single-index models in fixed-grid integration")
}


#' Composite Simpson's rule
#'
#' @param values List of function values at grid points (each is a vector)
#' @param grid Numeric vector of grid points
#'
#' @return Integral approximation (vector of same length as each element in values)
#' @keywords internal
composite_simpson <- function(values, grid) {
    n <- length(grid)
    if (n < 3) {
        # Fall back to trapezoid for small grids
        return(composite_trapezoid(values, grid))
    }
    
    # Ensure odd number of points for Simpson's rule
    if (n %% 2 == 0) {
        # Use trapezoid for last interval
        simpson_n <- n - 1
        values_simpson <- values[1:simpson_n]
        grid_simpson <- grid[1:simpson_n]
        
        result_simpson <- composite_simpson_odd(values_simpson, grid_simpson)
        
        # Add last trapezoid interval
        h_last <- grid[n] - grid[n-1]
        last_interval <- h_last / 2 * (values[[n-1]] + values[[n]])
        
        return(result_simpson + last_interval)
    }
    
    composite_simpson_odd(values, grid)
}


#' Composite Simpson's rule for odd number of points
#'
#' @keywords internal
composite_simpson_odd <- function(values, grid) {
    n <- length(grid)
    stopifnot(n %% 2 == 1)
    
    # Collect all interval contributions
    # This approach avoids accumulation issues with vector arithmetic
    intervals <- list()
    
    for (i in seq(1, n-2, by = 2)) {
        h1 <- grid[i+1] - grid[i]
        h2 <- grid[i+2] - grid[i+1]
        
        if (abs(h1 - h2) < 1e-10) {
            # Uniform spacing: classic Simpson's rule
            h <- h1
            interval <- (h / 3) * (values[[i]] + 4 * values[[i+1]] + values[[i+2]])
        } else {
            # Non-uniform spacing: use weighted Simpson's formula
            # Compute all coefficients as scalars first
            coef1 <- (h1 + h2) / 6 * (2 - h2/h1)
            coef2 <- (h1 + h2) / 6 * ((h1 + h2)^2 / (h1 * h2))
            coef3 <- (h1 + h2) / 6 * (2 - h1/h2)
            
            # Apply coefficients to vectors
            v1 <- coef1 * values[[i]]
            v2 <- coef2 * values[[i+1]]
            v3 <- coef3 * values[[i+2]]
            
            # Sum vectors
            interval <- v1 + v2 + v3
        }
        intervals[[length(intervals) + 1]] <- interval
    }
    
    # Sum all intervals
    if (length(intervals) == 0) {
        return(values[[1]] * 0)
    }
    
    result <- intervals[[1]]
    if (length(intervals) > 1) {
        for (i in 2:length(intervals)) {
            result <- result + intervals[[i]]
        }
    }
    
    result
}


#' Composite trapezoid rule
#'
#' @param values List of function values at grid points
#' @param grid Numeric vector of grid points
#'
#' @return Integral approximation
#' @keywords internal
composite_trapezoid <- function(values, grid) {
    n <- length(grid)
    if (n < 2) return(values[[1]] * 0)
    
    # Collect all intervals
    intervals <- list()
    for (i in 1:(n-1)) {
        h <- grid[i+1] - grid[i]
        interval <- (h / 2) * (values[[i]] + values[[i+1]])
        intervals[[i]] <- interval
    }
    
    # Sum all intervals
    result <- intervals[[1]]
    if (length(intervals) > 1) {
        for (i in 2:length(intervals)) {
            result <- result + intervals[[i]]
        }
    }
    
    result
}
