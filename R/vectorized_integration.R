#' Vectorized Term2 Integration for Multiple Alpha Values
#'
#' This is an experimental implementation that integrates term2 influence
#' across multiple alpha values simultaneously, sharing expensive computations
#' for improved performance.
#'
#' @param patient_data data.frame with patient's observations
#' @param outcome_model The fitted outcome model (any type supported by compute_SensIAT_expected_values)
#' @param base SplineBasis for marginal mean model  
#' @param alpha_vec Numeric vector of sensitivity parameters
#' @param marginal_beta Coefficients for the marginal mean spline basis
#' @param V_inv Inverse Gram matrix for base
#' @param tmin Lower integration bound
#' @param tmax Upper integration bound
#' @param impute_fn Function to impute data at time t
#' @param inv_link Inverse link function
#' @param tol Integration tolerance
#'
#' @return List with one element per alpha value, each containing:
#'   \item{Q}{Integration result (vector of length ncol(base))}
#'   \item{fcnt}{Number of function evaluations}
#'   \item{alpha}{Alpha value}
#'   \item{converged}{Convergence status}
#'
#' @export
compute_term2_influence_vectorized <- function(
    patient_data,
    outcome_model,
    base,
    alpha_vec,
    marginal_beta,
    V_inv,
    tmin,
    tmax,
    impute_fn,
    inv_link,
    tol = 1.490116e-08
) {
  stopifnot(is.numeric(alpha_vec), length(alpha_vec) > 0)
  stopifnot(is.function(impute_fn), is.function(inv_link))
  
  # Create wrapper function for compute_SensIAT_expected_values that handles vectorized alpha
  compute_expected_values_fn <- function(alpha, new.data) {
    compute_SensIAT_expected_values(
      model = outcome_model,
      alpha = alpha,
      new.data = new.data
    )
  }
  
  # Weight function - just the spline basis (NOT weighted by V_inv)
  # The original method uses: B(t) * E[Y|X(t)]
  weight_fn <- function(t) {
    as.vector(pcoriaccel_evaluate_basis(base, t))
  }
  
  # Marginal mean function - NOT NEEDED for term2, but keep for API compatibility
  # Set to always return 0 since we don't subtract it
  marginal_mean_fn <- function(t) {
    0.0
  }
  
  # PIECEWISE INTEGRATION: Match the original method's approach
  # Extract observation times and clip to integration bounds [tmin, tmax]
  observation_times <- sort(unique(patient_data$Time))
  times <- unique(c(tmin, observation_times[observation_times > tmin & observation_times < tmax], tmax))
  
  # Integrate piecewise between consecutive observation times
  # For each piece, impute at boundaries and linearly interpolate
  period_results <- purrr::map2(
    head(times, -1), 
    tail(times, -1), 
    function(lower, upper) {
      
      # Impute at piece boundaries (like original method)
      lower_data <- impute_fn(lower, patient_data)
      upper_data <- impute_fn(upper, patient_data)
      
      # Get Xbeta at boundaries using outcome model
      # Extract model formula and compute linear predictors
      model_frame_lower <- model.frame(outcome_model, data = lower_data)
      model_frame_upper <- model.frame(outcome_model, data = upper_data)
      
      xb_lower <- as.numeric(model.matrix(terms(outcome_model), data = model_frame_lower) %*% outcome_model$coef)
      xb_upper <- as.numeric(model.matrix(terms(outcome_model), data = model_frame_upper) %*% outcome_model$coef)
      
      # Create interpolating impute function for this piece
      # This matches the original: a = (time - lower)/(upper - lower); xb_time = (1-a)*xb_lower + a*xb_upper
      piece_impute_fn <- function(t, data) {
        # Linear interpolation of Xbeta
        a <- (t - lower) / (upper - lower)
        xb_interpolated <- (1 - a) * xb_lower + a * xb_upper
        
        # Return a data frame with the interpolated linear predictor
        # The compute_expected_values_fn will use this
        data.frame(xb = xb_interpolated)
      }
      
      # For this piece, we need a special expected values function that uses xb directly
      # instead of computing it from covariates
      piece_compute_expected_fn <- function(alpha, new.data) {
        # Extract xb from the interpolated data
        xb <- new.data$xb
        
        # Use the outcome model's parameters to compute E[Y|xb]
        # This replicates the pmf calculation from the original method
        mf <- model.frame(outcome_model)
        Xi <- model.matrix(terms(outcome_model), data = mf)
        beta <- outcome_model$coef
        Yi <- model.response(mf)
        y <- sort(unique(Yi))
        Xb <- Xi %*% beta
        
        pmf <- pcoriaccel_estimate_pmf(Xb, Yi, xb, y, outcome_model$bandwidth)
        
        # Handle vectorized alpha: compute for each alpha value
        if (length(alpha) > 1) {
          results <- purrr::map_dfr(alpha, function(a) {
            E_exp_alphaY  <- sum(exp(a * y) * pmf)
            E_Yexp_alphaY <- sum(y * exp(a * y) * pmf)
            E_Y_past <- E_Yexp_alphaY / E_exp_alphaY
            
            data.frame(
              alpha = a,
              E_Y_past = E_Y_past,
              E_exp_alphaY = E_exp_alphaY,
              E_Yexp_alphaY = E_Yexp_alphaY
            )
          })
          return(results)
        }
        
        # Single alpha case
        E_exp_alphaY  <- sum(exp(alpha * y) * pmf)
        E_Yexp_alphaY <- sum(y * exp(alpha * y) * pmf)
        E_Y_past <- E_Yexp_alphaY / E_exp_alphaY
        
        # Return data frame matching compute_SensIAT_expected_values format
        data.frame(
          alpha = alpha,
          E_Y_past = E_Y_past,
          E_exp_alphaY = E_exp_alphaY,
          E_Yexp_alphaY = E_Yexp_alphaY
        )
      }
      
      integrate_term2_vectorized_alpha(
        compute_expected_values_fn = piece_compute_expected_fn,
        impute_fn = piece_impute_fn,
        weight_fn = weight_fn,
        marginal_mean_fn = marginal_mean_fn,
        alpha_vec = alpha_vec,
        tmin = lower,
        tmax = upper,
        patient_data = patient_data,  # Not actually used now since piece_impute_fn doesn't need it
        tol = tol
      )
    }
  )
  
  # Sum up the Q values across all periods for each alpha
  # period_results is a list of length (# pieces), each containing results for all alphas
  # We need to transpose this to get one result per alpha
  n_alphas <- length(alpha_vec)
  results <- vector("list", n_alphas)
  
  for (i in seq_along(alpha_vec)) {
    # Sum Q values across all periods for this alpha
    Q_total <- Reduce(`+`, lapply(period_results, function(pr) pr[[i]]$Q))
    
    # Sum function evaluation counts
    fcnt_total <- sum(sapply(period_results, function(pr) pr[[i]]$fcnt))
    
    # Check if all pieces converged
    all_converged <- all(sapply(period_results, function(pr) pr[[i]]$converged))
    
    results[[i]] <- list(
      Q = Q_total,
      fcnt = fcnt_total,
      alpha = alpha_vec[i],
      converged = all_converged
    )
  }
  
  # Set names for easy access
  names(results) <- paste0("alpha_", alpha_vec)
  
  return(results)
}


#' Compare Vectorized vs Original Term2 Integration Methods
#'
#' Test function to verify that the new vectorized method produces 
#' the same results as the original method.
#'
#' @inheritParams compute_term2_influence_vectorized
#' @param method_original Function implementing original integration method
#' @param tolerance Numerical tolerance for comparison
#'
#' @return List with comparison results and diagnostics
#'
#' @export
compare_term2_methods <- function(
    patient_data,
    outcome_model,
    base,
    alpha_vec,
    marginal_beta,
    V_inv,
    tmin,
    tmax,
    impute_fn,
    inv_link,
    method_original = compute_term2_influence_original,
    tolerance = 1e-10
) {
  
  # Time the vectorized method
  time_vectorized <- system.time({
    results_vectorized <- compute_term2_influence_vectorized(
      patient_data, outcome_model, base, alpha_vec, marginal_beta, 
      V_inv, tmin, tmax, impute_fn, inv_link
    )
  })
  
  # Time the original method for each alpha
  time_original <- system.time({
    results_original <- purrr::map(alpha_vec, function(alpha) {
      list(
        Q = method_original(
          patient_data, outcome_model, base, alpha, marginal_beta,
          V_inv, tmin, tmax, impute_fn, inv_link
        ),
        alpha = alpha
      )
    })
  })
  
  # Compare results
  comparisons <- purrr::map2(results_vectorized, results_original, function(vec_res, orig_res) {
    diff <- abs(vec_res$Q - orig_res$Q)
    max_diff <- max(diff)
    mean_diff <- mean(diff)
    
    list(
      alpha = vec_res$alpha,
      max_absolute_diff = max_diff,
      mean_absolute_diff = mean_diff,
      all_close = all(diff < tolerance),
      vectorized_fcnt = vec_res$fcnt,
      vectorized_converged = vec_res$converged
    )
  })
  
  names(comparisons) <- paste0("alpha_", alpha_vec)
  
  return(list(
    comparisons = comparisons,
    timing = list(
      vectorized_elapsed = time_vectorized[["elapsed"]],
      original_elapsed = time_original[["elapsed"]],
      speedup = time_original[["elapsed"]] / time_vectorized[["elapsed"]]
    ),
    summary = list(
      all_alphas_close = all(purrr::map_lgl(comparisons, ~ .x$all_close)),
      max_diff_overall = max(purrr::map_dbl(comparisons, ~ .x$max_absolute_diff)),
      mean_diff_overall = mean(purrr::map_dbl(comparisons, ~ .x$mean_absolute_diff))
    )
  ))
}