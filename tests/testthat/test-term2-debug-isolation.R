test_that("isolate and debug term2 computation only", {
  skip_on_cran()
  
  # Create minimal reproducible data
  set.seed(12345)
  n_patients <- 2  # Minimal
  
  sim_data <- lapply(1:n_patients, function(id) {
    times <- c(0, 100, 200)  # Simple regular times
    outcomes <- c(10, 11, 12)  # Simple linear outcomes
    data.frame(
      Subject_ID = id,
      Time = times,
      Outcome = outcomes,
      Visit = 1:3
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
  
  cat("\n=== Input Data ===\n")
  print(data_with_lags)
  
  # Fit models
  intensity.model <- survival::coxph(
    Surv(..prev_time.., Time, !is.na(Outcome)) ~ 1,  # Null for simplicity
    data = data_with_lags %>% dplyr::filter(Time > 0)
  )
  
  outcome.model <- fit_SensIAT_single_index_fixed_coef_model(
    Outcome ~ ..prev_outcome.. + ..delta_time.. - 1,
    id = Subject_ID, 
    initial = c(0,0),  # Start at zero for simplicity
    data = data_with_lags %>% dplyr::filter(Time > 0)
  )
  
  cat("\n=== Outcome Model Summary ===\n")
  print(summary(outcome.model))
  
  cat("\n=== Intensity Model Summary ===\n")
  print(summary(intensity.model))
  
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
  
  knots <- c(10, 100, 200)
  
  cat("\n=== Fitting with term2_method = 'fast' ===\n")
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
  
  cat("\n=== Fitting with term2_method = 'original' ===\n")
  result_original <- fit_SensIAT_marginal_mean_model_generalized(
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
    term2_method = "original"
  )
  
  cat("\n=== Fitting with term2_method = 'fixed_grid' (n=500) ===\n")
  result_fixed <- fit_SensIAT_marginal_mean_model_generalized(
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
    term2_grid_n = 500
  )
  
  # Extract and compare results
  coef_fast <- result_fast$coefficients[[1]]
  coef_original <- result_original$coefficients[[1]]
  coef_fixed <- result_fixed$coefficients[[1]]
  
  cat("\n=== Results Comparison ===\n")
  cat("fast:     ", coef_fast, "\n")
  cat("original: ", coef_original, "\n")
  cat("fixed:    ", coef_fixed, "\n")
  
  cat("\n=== Error Analysis ===\n")
  error_original <- coef_original - coef_fast
  error_fixed <- coef_fixed - coef_fast
  rel_error_original <- abs(error_original / coef_fast)
  rel_error_fixed <- abs(error_fixed / coef_fast)
  
  cat("Absolute error (original - fast): ", error_original, "\n")
  cat("Absolute error (fixed - fast): ", error_fixed, "\n")
  cat("Relative error (original): ", rel_error_original, "\n")
  cat("Relative error (fixed): ", rel_error_fixed, "\n")
  
  cat("\n=== Objective Values ===\n")
  cat("fast loss:     ", result_fast$loss, "\n")
  cat("original loss: ", result_original$loss, "\n")
  cat("fixed loss:    ", result_fixed$loss, "\n")
  
  # Check if original matches fixed with coarse grid
  result_fixed_coarse <- fit_SensIAT_marginal_mean_model_generalized(
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
    term2_grid_n = 3  # Match number of observation points
  )
  
  coef_fixed_coarse <- result_fixed_coarse$coefficients[[1]]
  
  cat("\n=== Comparison with Coarse Fixed Grid (n=3) ===\n")
  cat("fixed(n=3):  ", coef_fixed_coarse, "\n")
  cat("original:    ", coef_original, "\n")
  cat("Difference:  ", coef_original - coef_fixed_coarse, "\n")
  
  # Try with observed times only
  result_fixed_obs <- fit_SensIAT_marginal_mean_model_generalized(
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
    term2_grid_n = 2  # Just start and end
  )
  
  coef_fixed_obs <- result_fixed_obs$coefficients[[1]]
  
  cat("\n=== Comparison with Fixed Grid (n=2, endpoints only) ===\n")
  cat("fixed(n=2):  ", coef_fixed_obs, "\n")
  cat("original:    ", coef_original, "\n")
  cat("Difference:  ", coef_original - coef_fixed_obs, "\n")
})


test_that("check parameterization differences in term2 methods", {
  skip_on_cran()
  
  # Check if original method uses different parameterization
  # For example: different link function, different scale, different order of operations
  
  set.seed(54321)
  data_with_lags <- SensIAT_example_data %>%
    dplyr::group_by(Subject_ID) %>%
    dplyr::mutate(
      ..prev_outcome.. = dplyr::lag(Outcome, default = NA_real_, order_by = Time),
      ..prev_time.. = dplyr::lag(Time, default = 0, order_by = Time),
      ..delta_time.. = Time - dplyr::lag(Time, default = NA_real_, order_by = Time)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(Subject_ID <= 2)
  
  intensity.model <- survival::coxph(
    Surv(..prev_time.., Time, !is.na(Outcome)) ~ 1,
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
  
  knots <- c(100, 300)
  
  # Test 1: Different link functions with same method
  cat("\n=== Test different link functions with fast method ===\n")
  for (link_fn in c("identity", "log")) {
    result <- fit_SensIAT_marginal_mean_model_generalized(
      data = data_with_lags,
      time = data_with_lags$Time,
      id = data_with_lags$Subject_ID,
      alpha = 0,
      knots = knots,
      outcome.model = outcome.model,
      intensity.model = intensity.model,
      loss = "lp_mse",
      link = link_fn,
      impute_data = impute_fn,
      term2_method = "fast"
    )
    cat("Link =", link_fn, ": ", result$coefficients[[1]], "\n")
  }
  
  # Test 2: Different alpha values with original vs fast
  cat("\n=== Test different alpha values ===\n")
  alpha_values <- c(-0.1, 0, 0.1)
  
  for (alpha_val in alpha_values) {
    result_fast <- fit_SensIAT_marginal_mean_model_generalized(
      data = data_with_lags,
      time = data_with_lags$Time,
      id = data_with_lags$Subject_ID,
      alpha = alpha_val,
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
      time = data_with_lags$Time,
      id = data_with_lags$Subject_ID,
      alpha = alpha_val,
      knots = knots,
      outcome.model = outcome.model,
      intensity.model = intensity.model,
      loss = "lp_mse",
      link = "identity",
      impute_data = impute_fn,
      term2_method = "original"
    )
    
    cat("Alpha =", alpha_val, "\n")
    cat("  fast:     ", result_fast$coefficients[[1]], "\n")
    cat("  original: ", result_original$coefficients[[1]], "\n")
    cat("  diff:     ", result_original$coefficients[[1]] - result_fast$coefficients[[1]], "\n")
  }
})
