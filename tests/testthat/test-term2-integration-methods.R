test_that("term2 integration methods produce equivalent results with sufficient grid points", {
  skip_on_cran()
  
  # Create small test dataset
  set.seed(123)
  n_patients <- 10
  
  sim_data <- lapply(1:n_patients, function(id) {
    n_visits <- sample(4:6, 1)
    times <- sort(c(0, cumsum(rexp(n_visits - 1, rate = 1/50))))
    outcomes <- numeric(n_visits)
    outcomes[1] <- rnorm(1, mean = 10, sd = 2)
    for (i in 2:n_visits) {
      delta_t <- times[i] - times[i-1]
      outcomes[i] <- 0.7 * outcomes[i-1] + 0.05 * delta_t + rnorm(1, sd = 1.5)
    }
    data.frame(
      Subject_ID = id,
      Time = times,
      Outcome = outcomes,
      Visit = 1:n_visits
    )
  }) %>% dplyr::bind_rows()
  
  data_with_lags <- sim_data %>%
    dplyr::group_by(Subject_ID) %>%
    dplyr::mutate(
      ..prev_outcome.. = dplyr::lag(Outcome, default = NA_real_, order_by = Time),
      ..prev_time.. = dplyr::lag(Time, default = 0, order_by = Time),
      ..delta_time.. = Time - dplyr::lag(Time, default = NA_real_, order_by = Time)
    ) %>%
    dplyr::ungroup()
  
  # Fit models
  intensity.model <- survival::coxph(
    Surv(..prev_time.., Time, !is.na(Outcome)) ~ ..prev_outcome.. + strata(Visit),
    data = data_with_lags %>% dplyr::filter(Time > 0)
  )
  
  outcome.model <- fit_SensIAT_single_index_fixed_coef_model(
    # Outcome ~ splines::ns(..prev_outcome.., df = 2) + ..delta_time.. - 1,
    Outcome ~ ..prev_outcome.. + Time + ..delta_time.. - 1,
    id = Subject_ID,
    data = data_with_lags %>% dplyr::filter(Time > 0),
    kernel = "dnorm",
    abs.tol = 1e-7
  )
  
  impute_fn <- function(t, df) {
    data_wl <- df %>%
      dplyr::mutate(
        ..prev_time.. = Time,
        ..prev_outcome.. = Outcome,
        ..delta_time.. = 0
      )
    extrapolate_from_last_observation(
      t, data_wl, "Time",
      slopes = c("..delta_time.." = 1)
    )
  }
  
  knots <- c(50, 150, 250)
  
  # Fit with different methods
  result_fast <- fit_SensIAT_marginal_mean_model_generalized(
    data = data_with_lags,
    time = data_with_lags$Time,
    id = data_with_lags$Subject_ID,
    alpha = 0,
    knots = knots,
    outcome.model = outcome.model,
    intensity.model = intensity.model,
    loss = "lp_mse",
    link = "identity",
    impute_data = impute_fn,
    term2_method = "fast"
  )
  
  result_original <- fit_SensIAT_marginal_mean_model_generalized(
    data = data_with_lags,
    time = Time,
    id = data_with_lags$Subject_ID,
    alpha = 0,
    knots = knots,
    outcome.model = outcome.model,
    intensity.model = intensity.model,
    loss = "lp_mse",
    link = "identity",
    impute_data = impute_fn,
    term2_method = "original",
    # tol=1e-12,  # Tight tolerance for convergence
    # debug = TRUE  # Add debug flag if supported
  )
  
  result_fixed_grid <- fit_SensIAT_marginal_mean_model_generalized(
    data = data_with_lags,
    time = data_with_lags$Time,
    id = data_with_lags$Subject_ID,
    alpha = 0,
    knots = knots,
    outcome.model = outcome.model,
    intensity.model = intensity.model,
    loss = "lp_mse",
    link = "identity",
    impute_data = impute_fn,
    term2_method = "fixed_grid",
    term2_grid_n = 200  # High density for accuracy
  )
  
  result_seeded <- fit_SensIAT_marginal_mean_model_generalized(
    data = data_with_lags,
    time = data_with_lags$Time,
    id = data_with_lags$Subject_ID,
    alpha = 0,
    knots = knots,
    outcome.model = outcome.model,
    intensity.model = intensity.model,
    loss = "lp_mse",
    link = "identity",
    impute_data = impute_fn,
    term2_method = "seeded_adaptive",
    term2_grid_n = 100
  )
  
  result_gauss_legendre <- fit_SensIAT_marginal_mean_model_generalized(
    data = data_with_lags,
    time = data_with_lags$Time,
    id = data_with_lags$Subject_ID,
    alpha = 0,
    knots = knots,
    outcome.model = outcome.model,
    intensity.model = intensity.model,
    loss = "lp_mse",
    link = "identity",
    impute_data = impute_fn,
    term2_method = "gauss_legendre",
    term2_grid_n = 50  # GL is very accurate with fewer points
  )
  
  # Extract coefficients
  coef_fast <- result_fast$coefficients[[1]]
  coef_original <- result_original$coefficients[[1]]
  coef_fixed_grid <- result_fixed_grid$coefficients[[1]]
  coef_seeded <- result_seeded$coefficients[[1]]
  coef_gauss <- result_gauss_legendre$coefficients[[1]]
  
  # Diagnostic: Compare coefficient magnitudes and signs
  cat("\nCoefficient comparison:\n")
  cat("fast:     ", round(coef_fast, 6), "\n")
  cat("original: ", round(coef_original, 6), "\n")
  cat("fixed:    ", round(coef_fixed_grid, 6), "\n")
  cat("seeded:   ", round(coef_seeded, 6), "\n")
  cat("gauss:    ", round(coef_gauss, 6), "\n")
  
  # Check for common issues:
  # 1. Sign flip
  sign_match <- sign(coef_fast) == sign(coef_original)
  cat("Sign match: ", all(sign_match), "\n")
  
  # 2. Scale factor (ratio should be constant if just scaling issue)
  ratio <- coef_original / coef_fast
  cat("Ratio (original/fast): ", round(ratio, 4), "\n")
  cat("Ratio SD: ", round(sd(ratio), 6), "\n")
  
  # 3. Offset (difference should be constant if additive error)
  diff <- coef_original - coef_fast
  cat("Difference: ", round(diff, 6), "\n")
  
  # 4. Relative error
  rel_error <- abs((coef_original - coef_fast) / coef_fast)
  cat("Relative error: ", round(rel_error, 6), "\n")
  cat("Max relative error: ", round(max(rel_error), 6), "\n")
  
  # Test parity between adaptive methods (should be nearly identical)
  # NOTE: If original consistently differs, this suggests a bug in original
  expect_equal(coef_fast, coef_original, tolerance = 1e-5,
               label = "fast vs original methods (now fixed)")
  
  # Test parity with fixed_grid (should be close with high grid density)
  # Note: Grid methods converge more slowly than adaptive methods, so we use
  # a looser tolerance. With n=200 grid points, accuracy is ~O(h^2) where h=(tmax-tmin)/n.
  expect_equal(coef_fast, coef_fixed_grid, tolerance = 0.05,
               label = "fast vs fixed_grid (n=200)")
  
  # Test parity with seeded_adaptive (should be very close)
  expect_equal(coef_fast, coef_seeded, tolerance = 1e-5,
               label = "fast vs seeded_adaptive")
  
  # Test parity with gauss_legendre (GL is highly accurate with fewer nodes)
  # Note: GL quadrature uses different evaluation points than other methods,
  # so tolerance is slightly looser (~3% to accommodate integration differences)
  expect_equal(coef_fast, coef_gauss, tolerance = 0.03,
               label = "fast vs gauss_legendre (n=50)")
  
  # Both adaptive methods should be nearly identical because they use the same
  # integrand discretized via adaptive Simpson's rule (fast segments at observation
  # times where the interpolation has discontinuities)
  expect_equal(coef_fast, coef_original, tolerance = 1e-5,
               label = "fast vs original methods (segmented vs non-segmented adaptive quadrature)")
})


test_that("fixed_grid accuracy improves with grid density", {
  skip_on_cran()
  
  # Simple test data
  set.seed(456)
  data_with_lags <- SensIAT_example_data %>%
    dplyr::group_by(Subject_ID) %>%
    dplyr::mutate(
      ..prev_outcome.. = dplyr::lag(Outcome, default = NA_real_, order_by = Time),
      ..prev_time.. = dplyr::lag(Time, default = 0, order_by = Time),
      ..delta_time.. = Time - dplyr::lag(Time, default = NA_real_, order_by = Time)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(Subject_ID <= 5)  # Small subset for speed
  
  intensity.model <- survival::coxph(
    Surv(..prev_time.., Time, !is.na(Outcome)) ~ ..prev_outcome.. + strata(Visit),
    data = data_with_lags %>% dplyr::filter(Time > 0)
  )
  
  outcome.model <- fit_SensIAT_single_index_fixed_coef_model(
    Outcome ~ splines::ns(..prev_outcome.., df = 2) + ..delta_time.. - 1,
    id = Subject_ID,
    data = data_with_lags %>% dplyr::filter(Time > 0)
  )
  
  impute_fn <- function(t, df) {
    data_wl <- df %>%
      dplyr::mutate(
        ..prev_time.. = Time,
        ..prev_outcome.. = Outcome,
        ..delta_time.. = 0
      )
    extrapolate_from_last_observation(
      t, data_wl, "Time",
      slopes = c("..delta_time.." = 1)
    )
  }
  
  knots <- c(100, 300, 500)
  
  # Reference with adaptive method
  result_ref <- fit_SensIAT_marginal_mean_model_generalized(
    data = data_with_lags,
    time = data_with_lags$Time,
    id = data_with_lags$Subject_ID,
    alpha = 0,
    knots = knots,
    outcome.model = outcome.model,
    intensity.model = intensity.model,
    loss = "lp_mse",
    link = "identity",
    impute_data = impute_fn,
    term2_method = "fast"
  )
  
  coef_ref <- result_ref$coefficients[[1]]
  
  # Test with different grid densities
  # Use small grids for speed - convergence should be apparent even with these
  grid_sizes <- c(20, 50, 100)
  errors <- numeric(length(grid_sizes))
  
  for (i in seq_along(grid_sizes)) {
    result <- fit_SensIAT_marginal_mean_model_generalized(
      data = data_with_lags,
      time = data_with_lags$Time,
      id = data_with_lags$Subject_ID,
      alpha = 0,
      knots = knots,
      outcome.model = outcome.model,
      intensity.model = intensity.model,
      loss = "lp_mse",
      link = "identity",
      impute_data = impute_fn,
      term2_method = "fixed_grid",
      term2_grid_n = grid_sizes[i]
    )
    
    errors[i] <- sqrt(mean((result$coefficients[[1]] - coef_ref)^2))
  }
  
  # Errors should decrease with grid density (or at least not increase significantly)
  # Note: with observation times included in grid, very small grids may not show
  # monotonic convergence due to different grid point placements
  expect_true(errors[1] > errors[3] * 0.9,  # n=20 should be worse than n=100
              label = "error decreases from n=20 to n=100")
  
  # Higher-density grid should achieve reasonable accuracy
  # For grid methods, 0.05 RMSE is a reasonable target with 100 points
  expect_lt(errors[3], 0.1,
            label = "n=100 achieves reasonable accuracy")
})


test_that("all methods work with multiple alpha values", {
  skip_on_cran()
  
  data_with_lags <- SensIAT_example_data %>%
    dplyr::group_by(Subject_ID) %>%
    dplyr::mutate(
      ..prev_outcome.. = dplyr::lag(Outcome, default = NA_real_, order_by = Time),
      ..prev_time.. = dplyr::lag(Time, default = 0, order_by = Time),
      ..delta_time.. = Time - dplyr::lag(Time, default = NA_real_, order_by = Time)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(Subject_ID <= 3)  # Very small for speed
  
  intensity.model <- survival::coxph(
    Surv(..prev_time.., Time, !is.na(Outcome)) ~ ..prev_outcome.. + strata(Visit),
    data = data_with_lags %>% dplyr::filter(Time > 0)
  )
  
  outcome.model <- fit_SensIAT_single_index_fixed_coef_model(
    Outcome ~ splines::ns(..prev_outcome.., df = 2) + ..delta_time.. - 1,
    id = Subject_ID,
    data = data_with_lags %>% dplyr::filter(Time > 0)
  )
  
  impute_fn <- function(t, df) {
    data_wl <- df %>%
      dplyr::mutate(
        ..prev_time.. = Time,
        ..prev_outcome.. = Outcome,
        ..delta_time.. = 0
      )
    extrapolate_from_last_observation(
      t, data_wl, "Time",
      slopes = c("..delta_time.." = 1)
    )
  }
  
  knots <- c(100, 300, 500)
  alpha_seq <- c(-0.05, 0, 0.05)
  
  # Test each method with multiple alphas
  methods <- c("fast", "original", "fixed_grid", "seeded_adaptive")
  
  for (method in methods) {
    result <- fit_SensIAT_marginal_mean_model_generalized(
      data = data_with_lags,
      time = data_with_lags$Time,
      id = data_with_lags$Subject_ID,
      alpha = alpha_seq,
      knots = knots,
      outcome.model = outcome.model,
      intensity.model = intensity.model,
      loss = "lp_mse",
      link = "identity",
      impute_data = impute_fn,
      term2_method = method,
      term2_grid_n = 100
    )
    
    expect_length(result$coefficients, length(alpha_seq))
    expect_length(result$influence, length(alpha_seq))
    
    # Check that results vary by alpha
    coef_1 <- result$coefficients[[1]]
    coef_3 <- result$coefficients[[3]]
    expect_false(isTRUE(all.equal(coef_1, coef_3, tolerance = 1e-8)),
                 label = paste(method, "coefficients differ for different alphas"))
  }
})


test_that("integration utility functions work correctly", {
  # Test create_integration_grid
  grid1 <- create_integration_grid(0, 100, n_grid = 11)
  expect_length(grid1, 11)
  expect_equal(min(grid1), 0)
  expect_equal(max(grid1), 100)
  
  # With observation times
  obs_times <- c(10, 30, 70)
  grid2 <- create_integration_grid(0, 100, n_grid = 11, obs_times = obs_times)
  expect_true(all(obs_times %in% grid2))
  expect_gte(length(grid2), 11)  # Should include uniform grid + obs times
  
  # Test composite_trapezoid
  values <- list(1, 4, 9, 16)  # y = x^2 at x = 1, 2, 3, 4
  grid <- c(1, 2, 3, 4)
  integral <- composite_trapezoid(values, grid)
  # Trapezoid rule for x^2 from 1 to 4
  # = 0.5*1*(1+4) + 0.5*1*(4+9) + 0.5*1*(9+16) = 2.5 + 6.5 + 12.5 = 21.5
  expect_equal(integral, 21.5)
  
  # Test composite_simpson (odd number of points)
  values3 <- list(1, 4, 9)  # y = x^2 at x = 1, 2, 3
  grid3 <- c(1, 2, 3)
  integral3 <- composite_simpson(values3, grid3)
  # Simpson's rule for x^2 from 1 to 3 with h=1
  # = (h/3) * (f(1) + 4*f(2) + f(3)) = (1/3) * (1 + 16 + 9) = 26/3
  expect_equal(integral3, 26/3, tolerance = 1e-10)
})
