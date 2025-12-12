# Debug output control - set SENSIAT_TEST_DEBUG=1 to enable
test_debug <- Sys.getenv("SENSIAT_TEST_DEBUG", "0") == "1"

test_that("piecewise integration matches original method exactly", {
  # Use exact setup from test-vectorized-integration-simple.R
  model.data <- SensIAT_example_data |>
    dplyr::group_by(Subject_ID) |>
    dplyr::arrange(Time) |>
    dplyr::mutate(
      prev_time = dplyr::lag(Time),
      prev_outcome = dplyr::lag(Outcome),
      delta_time = Time - dplyr::lag(Time),
      visit.number = seq_along(Time)
    ) |>
    dplyr::filter(!is.na(Outcome))

  followup.data <- model.data |>
    dplyr::filter(Time > 0)

  intensity.model <-
    rlang::inject(coxph(
      Surv(prev_time,Time,!is.na(Outcome)) ~
        prev_outcome+strata(visit.number),
      id = Subject_ID,
      data = followup.data
    ))

  outcome.model <- fit_SensIAT_single_index_fixed_coef_model(
    Outcome ~
      splines::ns(prev_outcome, df=3) +
      prev_outcome +
      delta_time - 1,
    id = Subject_ID,
    data = followup.data
  )

  base <- orthogonalsplinebasis::SplineBasis(c(60,60,60,60,260,460,460,460,460))

  centering.statistics <-
    dplyr::summarize(
      dplyr::ungroup(dplyr::filter(model.data, Time > 0, !is.na(Outcome))),
      dplyr::across(c(Time, delta_time),
                    list(mean = ~ mean(.x, na.rm = TRUE),
                         sd = ~ sd(.x, na.rm = TRUE))
      )
    ) |>
    as.numeric() |>
    matrix(ncol = 2, byrow = TRUE) |>
    `dimnames<-`(list(c("time", "delta_time"), c("mean", "sd")))

  # Get first patient
  df_i <- model.data |> dplyr::filter(Subject_ID == 1) |> as.data.frame()

  variables <- list(
    time = rlang::sym("Time"),
    id = rlang::sym("Subject_ID"),
    outcome = rlang::sym("Outcome"),
    prev_time = rlang::sym("prev_time"),
    prev_outcome = rlang::sym("prev_outcome"),
    delta_time = rlang::sym("delta_time")
  )

  # Test alpha = 0 only
  alpha_test <- 0

  # Get integration bounds
  tmin <- base@knots[base@order]
  tmax <- base@knots[length(base@knots) - base@order + 1]

  # Get piecewise intervals
  observation_times <- sort(unique(df_i$Time))
  times <- unique(c(tmin, observation_times[observation_times > tmin & observation_times < tmax], tmax))

  if (test_debug) {
    cat("\nIntegration setup:\n")
    cat("  tmin:", tmin, "tmax:", tmax, "\n")
    cat("  Observation times:", observation_times, "\n")
    cat("  Piece boundaries:", times, "\n\n")
  }

  # For FIRST piece [60, 214], manually compute what original method does
  lower <- times[1]
  upper <- times[2]

  if (test_debug) {
    cat(sprintf("Testing piece [%.0f, %.0f]:\n", lower, upper))
  }

  # Impute at boundaries using impute_patient_df (available in testthat)
  lower_df <- impute_patient_df(
    eval.times = lower,
    df_i = df_i,
    variables = variables,
    centering = centering.statistics,
    right = FALSE
  )
  upper_df <- impute_patient_df(
    eval.times = upper,
    df_i = df_i,
    variables = variables,
    centering = centering.statistics,
    right = TRUE
  )

  # Get xb at boundaries
  xb_lower <- as.numeric(model.matrix(terms(outcome.model), data = lower_df) %*% outcome.model$coef)
  xb_upper <- as.numeric(model.matrix(terms(outcome.model), data = upper_df) %*% outcome.model$coef)

  if (test_debug) {
    cat(sprintf("  xb_lower = %.4f, xb_upper = %.4f\n", xb_lower, xb_upper))
  }

  # Test integrand at middle of piece
  t_mid <- (lower + upper) / 2
  a_mid <- (t_mid - lower) / (upper - lower)
  xb_mid_interp <- (1 - a_mid) * xb_lower + a_mid * xb_upper

  if (test_debug) {
    cat(sprintf("  At t=%.1f (midpoint): a=%.3f, xb_interpolated=%.4f\n", t_mid, a_mid, xb_mid_interp))
  }

  # Compute pmf and expected value at midpoint
  mf <- model.frame(outcome.model)
  Xi <- model.matrix(terms(outcome.model), data = mf)
  Yi <- model.response(mf)
  y <- sort(unique(Yi))
  Xb <- Xi %*% outcome.model$coef

  pmf_mid <- pcoriaccel_estimate_pmf(Xb, Yi, xb_mid_interp, y, outcome.model$bandwidth)
  E_exp_alphaY_mid <- sum(exp(alpha_test * y) * pmf_mid)
  E_Yexp_alphaY_mid <- sum(y * exp(alpha_test * y) * pmf_mid)
  E_Y_mid <- E_Yexp_alphaY_mid / E_exp_alphaY_mid

  cat(sprintf("  E[Y|X(%.1f)] = %.4f\n", t_mid, E_Y_mid))

  # Get weight at midpoint
  B_mid <- pcoriaccel_evaluate_basis(base, t_mid)
  cat(sprintf("  B(%.1f) = [%.4f, %.4f, %.4f, %.4f, %.4f]\n", t_mid, B_mid[1], B_mid[2], B_mid[3], B_mid[4], B_mid[5]))

  # Integrand at midpoint (for alpha=0, marginal_mean=0)
  integrand_mid <- B_mid * E_Y_mid
  if (test_debug) {
    cat(sprintf("  Integrand(%.1f) = [%.4f, %.4f, %.4f, %.4f, %.4f]\n\n", t_mid,
                integrand_mid[1], integrand_mid[2], integrand_mid[3], integrand_mid[4], integrand_mid[5]))
  }

  # Now test what the original pracma::quadv does for this piece
  integrand_fn_original <- function(time) {
    a <- (time - lower) / (upper - lower)
    xb_time <- (1 - a) * xb_lower + a * xb_upper

    pmf <- pcoriaccel_estimate_pmf(Xb, Yi, xb_time, y, outcome.model$bandwidth)
    E_exp_alphaY <- sum(exp(alpha_test * y) * pmf)
    E_Yexp_alphaY <- sum(y * exp(alpha_test * y) * pmf)
    ev <- E_Yexp_alphaY / E_exp_alphaY

    B_t <- pcoriaccel_evaluate_basis(base, time)
    return(as.vector(B_t * ev))
  }

  # Integrate using pracma::quadv (original method)
  original_integral <- pracma::quadv(
    integrand_fn_original,
    a = lower,
    b = upper,
    tol = 1.490116e-08
  )

  cat(sprintf("Original (pracma::quadv) integral over [%.0f, %.0f]:\n", lower, upper))
  cat(sprintf("  Q = [%.4f, %.4f, %.4f, %.4f, %.4f]\n\n", 
              original_integral$Q[1], original_integral$Q[2], original_integral$Q[3], 
              original_integral$Q[4], original_integral$Q[5]))

  # Now test the vectorized C++ integrator for this piece
  # Set up the impute and compute functions as in vectorized_integration.R

  piece_impute_fn <- function(t, data) {
    a <- (t - lower) / (upper - lower)
    xb_interpolated <- (1 - a) * xb_lower + a * xb_upper
    data.frame(xb = xb_interpolated)
  }

  piece_compute_expected_fn <- function(alpha, new.data) {
    xb <- new.data$xb

    pmf <- pcoriaccel_estimate_pmf(Xb, Yi, xb, y, outcome.model$bandwidth)

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

    E_exp_alphaY  <- sum(exp(alpha * y) * pmf)
    E_Yexp_alphaY <- sum(y * exp(alpha * y) * pmf)
    E_Y_past <- E_Yexp_alphaY / E_exp_alphaY

    data.frame(
      alpha = alpha,
      E_Y_past = E_Y_past,
      E_exp_alphaY = E_exp_alphaY,
      E_Yexp_alphaY = E_Yexp_alphaY
    )
  }

  weight_fn <- function(t) {
    as.vector(pcoriaccel_evaluate_basis(base, t))
  }

  marginal_mean_fn <- function(t) {
    0.0
  }

  # Test that the integrand functions match at midpoint
  cat("\nTesting integrand at midpoint t=137:\n")
  
  # Original integrand
  original_integrand_mid <- integrand_fn_original(t_mid)
  cat(sprintf("  Original: [%.4f, %.4f, %.4f, %.4f, %.4f]\n",
              original_integrand_mid[1], original_integrand_mid[2], original_integrand_mid[3],
              original_integrand_mid[4], original_integrand_mid[5]))
  
  # Vectorized integrand - need to manually call the integrand computation
  imputed_mid <- piece_impute_fn(t_mid, df_i)
  expected_mid <- piece_compute_expected_fn(alpha_test, imputed_mid)
  weight_mid <- weight_fn(t_mid)
  marginal_mid <- marginal_mean_fn(t_mid)
  
  # Match C++ computation
  conditional_mean_mid <- expected_mid$E_Yexp_alphaY / expected_mid$E_exp_alphaY
  integrand_scalar_mid <- conditional_mean_mid - marginal_mid
  vectorized_integrand_mid <- weight_mid * integrand_scalar_mid
  
  cat(sprintf("  Vectorized: [%.4f, %.4f, %.4f, %.4f, %.4f]\n",
              vectorized_integrand_mid[1], vectorized_integrand_mid[2], vectorized_integrand_mid[3],
              vectorized_integrand_mid[4], vectorized_integrand_mid[5]))
  
  cat(sprintf("  Ratio: [%.4f, %.4f, %.4f, %.4f, %.4f]\n\n",
              vectorized_integrand_mid[1] / original_integrand_mid[1],
              vectorized_integrand_mid[2] / original_integrand_mid[2],
              vectorized_integrand_mid[3] / original_integrand_mid[3],
              vectorized_integrand_mid[4] / original_integrand_mid[4],
              vectorized_integrand_mid[5] / original_integrand_mid[5]))

  # Sanity check: compute integral using simple trapezoidal rule with fine grid
  cat("Sanity check - trapezoidal rule with 1000 points:\n")
  t_grid <- seq(lower, upper, length.out = 1000)
  integrand_grid <- t(sapply(t_grid, integrand_fn_original))
  trap_integral <- colSums((integrand_grid[-1,] + integrand_grid[-nrow(integrand_grid),]) / 2) * (upper - lower) / (length(t_grid) - 1)
  cat(sprintf("  Q_trap = [%.4f, %.4f, %.4f, %.4f, %.4f]\n\n",
              trap_integral[1], trap_integral[2], trap_integral[3],
              trap_integral[4], trap_integral[5]))

  vectorized_integral <- integrate_term2_vectorized_alpha(
    compute_expected_values_fn = piece_compute_expected_fn,
    impute_fn = piece_impute_fn,
    weight_fn = weight_fn,
    marginal_mean_fn = marginal_mean_fn,
    alpha_vec = c(alpha_test),
    tmin = lower,
    tmax = upper,
    patient_data = df_i,
    tol = 1.490116e-08
  )

  vectorized_Q <- vectorized_integral[[1]]$Q

  cat(sprintf("Vectorized (C++ adaptive Simpson) integral over [%.0f, %.0f]:\n", lower, upper))
  cat(sprintf("  Q = [%.4f, %.4f, %.4f, %.4f, %.4f]\n", 
              vectorized_Q[1], vectorized_Q[2], vectorized_Q[3], 
              vectorized_Q[4], vectorized_Q[5]))
  cat(sprintf("  Function evaluations: %d\n", vectorized_integral[[1]]$fcnt))
  cat(sprintf("  Converged: %s\n\n", vectorized_integral[[1]]$converged))

  # Compare
  diff <- original_integral$Q - vectorized_Q
  cat(sprintf("Difference (Original - Vectorized):\n"))
  cat(sprintf("  [%.4f, %.4f, %.4f, %.4f, %.4f]\n\n", diff[1], diff[2], diff[3], diff[4], diff[5]))

  # They should match within tolerance
  expect_true(max(abs(diff)) < 1e-5, 
              info = sprintf("Max difference: %.8f", max(abs(diff))))
})
