# Gauss-Legendre quadrature integration for term2 influence computation
#
# This method uses Gauss-Legendre quadrature (via statmod::gauss.quad) for
# numerical integration of the term2 influence function. Gauss-Legendre
# quadrature is highly accurate for smooth integrands and achieves
# high accuracy with relatively few evaluation points.
#
# Performance characteristics:
# - Accuracy: Exact for polynomials up to degree 2n-1 using n points
# - Efficiency: Fewer function evaluations than Simpson's rule for same accuracy
# - Best for: Smooth integrands without discontinuities
# - Limitation: Fixed evaluation points (cannot incorporate observation times directly)

#' Compute term2 influence using Gauss-Legendre quadrature
#'
#' This method uses Gauss-Legendre quadrature for numerical integration.
#' The integrand is evaluated at Gauss-Legendre nodes, which are optimal
#' for polynomial approximation of smooth functions.
#'
#' @param patient_data data.frame with patient's observations
#' @param outcome_model The fitted outcome model
#' @param base SplineBasis for marginal mean model
#' @param alpha Sensitivity parameter
#' @param marginal_beta Coefficients (beta) for the marginal mean spline basis
#' @param V_inv Inverse Gram matrix for base
#' @param tmin Lower integration bound
#' @param tmax Upper integration bound
#' @param impute_fn Function to impute data at time t
#' @param inv_link Inverse link function
#' @param W Weight function W(t, beta)
#' @param expected_grid Optional pre-computed expected values (not typically used
#'        for Gauss-Legendre since nodes are specialized)
#' @param n_nodes Number of Gauss-Legendre nodes (default 50)
#' @param time_var Time variable name (optional)
#' @param ... Additional arguments (not used)
#'
#' @return Numeric vector of length ncol(base) with term2 influence values
#'
#' @details
#' Gauss-Legendre quadrature approximates the integral as a weighted sum:
#' \deqn{\int_a^b f(x) dx \approx \frac{b-a}{2} \sum_{i=1}^n w_i f(x_i)}
#' where \eqn{x_i} are the Gauss-Legendre nodes transformed to the interval
#' \eqn{(a, b)} and \eqn{w_i} are the corresponding weights.
#'
#' The method requires the `statmod` package for computing nodes and weights.
#'
#' @keywords internal
compute_term2_influence_gauss_legendre <- function(
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
  n_nodes = 50,
  time_var = NULL,
  ...
) {
    # Check if we have pre-computed Gauss-Legendre data
    if (!is.null(expected_grid) && !is.null(expected_grid$weights)) {
        # Use pre-computed GL nodes, weights, and expected values
        nodes <- expected_grid$grid
        weights <- expected_grid$weights
        B_nodes <- expected_grid$B_grid
        E_nodes <- expected_grid$E_grid
    } else {
        # Compute our own GL nodes and values
        if (!requireNamespace("statmod", quietly = TRUE)) {
            stop("Package 'statmod' is required for Gauss-Legendre quadrature. ",
                 "Install it with: install.packages('statmod')")
        }
        
        # Get Gauss-Legendre nodes and weights on [-1, 1]
        gl <- statmod::gauss.quad(n = n_nodes, kind = "legendre")
        
        # Transform nodes from [-1, 1] to [tmin, tmax]
        half_range <- (tmax - tmin) / 2
        mid_point <- (tmax + tmin) / 2
        nodes <- half_range * gl$nodes + mid_point
        weights <- gl$weights * half_range  # Include Jacobian
        
        # Pre-compute basis evaluations at nodes
        B_nodes <- lapply(nodes, function(t) pcoriaccel_evaluate_basis(base, t))
        
        # Compute expected values at Gauss-Legendre nodes
        E_nodes <- compute_expected_values_at_nodes(
            nodes = nodes,
            patient_data = patient_data,
            outcome_model = outcome_model,
            alpha = alpha,
            impute_fn = impute_fn,
            time_var = time_var
        )
    }
    
    # Compute weighted integrand values at nodes
    result <- rep(0, ncol(base))
    
    for (i in seq_along(nodes)) {
        t <- nodes[i]
        B <- B_nodes[[i]]
        E <- E_nodes[[i]]
        w_gl <- weights[i]
        
        # Extract scalar values
        E_exp <- as.numeric(E$E_exp_alphaY)[1]
        E_Yexp <- as.numeric(E$E_Yexp_alphaY)[1]
        
        if (!is.list(E) || !is.finite(E_exp) || 
            !is.finite(E_Yexp) || E_exp <= 0) {
            next
        }
        
        # Compute weight from W function
        weight <- W(t, marginal_beta)
        eta <- sum(B * marginal_beta)
        mu_hat <- inv_link(eta)
        
        # Ensure weight is a numeric vector
        weight <- as.numeric(weight)
        # Compute scalar difference
        diff_scalar <- as.numeric(E_Yexp / E_exp - mu_hat)
        
        # Accumulate weighted contribution
        result <- result + w_gl * weight * diff_scalar
    }
    
    # Ensure finite output
    result <- ifelse(is.finite(result), result, 0)
    result
}


#' Compute expected values at Gauss-Legendre nodes
#'
#' @param nodes Numeric vector of Gauss-Legendre nodes
#' @param patient_data Patient data frame
#' @param outcome_model Fitted outcome model
#' @param alpha Sensitivity parameter
#' @param impute_fn Imputation function
#' @param time_var Time variable name (optional)
#'
#' @return List of expected value lists at each node
#' @keywords internal
compute_expected_values_at_nodes <- function(
    nodes,
    patient_data,
    outcome_model,
    alpha,
    impute_fn = NULL,
    time_var = NULL
) {
    # Use impute_fn when available (preferred path)
    if (!is.null(impute_fn)) {
        return(lapply(nodes, function(t) {
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
    
    # Fallback for simple single-index models
    if (is(outcome_model, "SensIAT::Single-index-outcome-model") &&
        !is.null(attr(outcome_model, "kernel")) &&
        !is.null(outcome_model$bandwidth)) {
        
        # Extract patient times and outcomes
        outcome_var <- as.character(rlang::f_lhs(formula(outcome_model)))
        patient_times <- extract_patient_times(patient_data, time_var)
        patient_outcomes <- patient_data[[outcome_var]]
        
        if (length(patient_times) == 0 || length(patient_outcomes) == 0) {
            return(lapply(nodes, function(t) {
                list(E_exp_alphaY = 1, E_Yexp_alphaY = 0)
            }))
        }
        
        # Get the last observation before integration region
        valid_idx <- which(!is.na(patient_outcomes))
        if (length(valid_idx) == 0) {
            return(lapply(nodes, function(t) {
                list(E_exp_alphaY = 1, E_Yexp_alphaY = 0)
            }))
        }
        
        last_obs_idx <- valid_idx[length(valid_idx)]
        prev_outcome <- patient_outcomes[last_obs_idx]
        prev_time <- patient_times[last_obs_idx]
        
        return(lapply(nodes, function(t) {
            delta_time <- t - prev_time
            if (delta_time < 0) delta_time <- 0
            
            # Build template for model matrix
            template_df <- patient_data[1, , drop = FALSE]
            template_df[["..prev_outcome.."]] <- prev_outcome
            template_df[["..delta_time.."]] <- delta_time
            
            ev <- compute_SensIAT_expected_values(
                model = outcome_model,
                alpha = alpha,
                new.data = template_df
            )
            list(E_exp_alphaY = as.numeric(ev$E_exp_alphaY)[1],
                 E_Yexp_alphaY = as.numeric(ev$E_Yexp_alphaY)[1])
        }))
    }
    
    # Last resort fallback
    lapply(nodes, function(t) {
        list(E_exp_alphaY = 1, E_Yexp_alphaY = 0)
    })
}


#' Pre-compute Gauss-Legendre grid for term2 integration
#'
#' Similar to pre_compute_term2_grid but optimized for Gauss-Legendre quadrature.
#' Pre-computes the expected values at Gauss-Legendre nodes for reuse across
#' beta iterations.
#'
#' @param patient_data Patient data frame
#' @param outcome_model Fitted outcome model
#' @param base SplineBasis object
#' @param alpha Sensitivity parameter
#' @param tmin Lower integration bound
#' @param tmax Upper integration bound
#' @param impute_fn Imputation function
#' @param n_nodes Number of Gauss-Legendre nodes
#' @param time_var Time variable name (optional)
#'
#' @return List with:
#'   - nodes: Gauss-Legendre nodes transformed to the interval (tmin, tmax)
#'   - weights: Gauss-Legendre weights (including Jacobian)
#'   - B_nodes: Basis evaluations at nodes
#'   - E_nodes: Expected values at nodes
#'
#' @keywords internal
pre_compute_gauss_legendre_grid <- function(
    patient_data,
    outcome_model,
    base,
    alpha,
    tmin,
    tmax,
    impute_fn = NULL,
    n_nodes = 50,
    time_var = NULL
) {
    if (!requireNamespace("statmod", quietly = TRUE)) {
        stop("Package 'statmod' is required for Gauss-Legendre quadrature")
    }
    
    # Get Gauss-Legendre nodes and weights
    gl <- statmod::gauss.quad(n = n_nodes, kind = "legendre")
    
    # Transform to [tmin, tmax]
    half_range <- (tmax - tmin) / 2
    mid_point <- (tmax + tmin) / 2
    nodes <- half_range * gl$nodes + mid_point
    weights <- gl$weights * half_range
    
    # Pre-compute basis evaluations
    B_nodes <- lapply(nodes, function(t) pcoriaccel_evaluate_basis(base, t))
    
    # Pre-compute expected values
    E_nodes <- compute_expected_values_at_nodes(
        nodes = nodes,
        patient_data = patient_data,
        outcome_model = outcome_model,
        alpha = alpha,
        impute_fn = impute_fn,
        time_var = time_var
    )
    
    list(
        grid = nodes,  # Use 'grid' for consistency with other methods
        nodes = nodes,
        weights = weights,
        B_grid = B_nodes,
        B_nodes = B_nodes,
        E_grid = E_nodes,
        E_nodes = E_nodes
    )
}
