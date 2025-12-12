test_that("vectorized integration matches original for linear model", {
  skip("Deprecated test - API changed to use variables/centering instead of impute_fn/inv_link. See test-vectorized-integration-simple.R for current tests.")
  skip_on_cran()
  
  # Generate synthetic data for testing
  set.seed(42)
  n_obs <- 20
  patient_data <- data.frame(
    Time = sort(runif(n_obs, 0, 100)),
    Outcome = rnorm(n_obs, mean = 10, sd = 2),
    prev_outcome_lag1 = c(NA, rnorm(n_obs-1, mean = 9.5, sd = 1.8)),
    delta_time_lag1 = c(NA, diff(sort(runif(n_obs, 0, 100))))
  )
  
  # Remove first row (has NAs)
  patient_data <- patient_data[-1, ]
  
  # Create a simple linear outcome model for testing
  outcome_model <- lm(Outcome ~ prev_outcome_lag1 + delta_time_lag1, data = patient_data)
  
  # Create spline basis for marginal model
  knots <- c(25, 25, 25, 25, 50, 75, 75, 75, 75)
  base <- orthogonalsplinebasis::SplineBasis(knots = knots, order = 4L)
  
  # Generate synthetic marginal beta coefficients
  n_basis <- ncol(base) 
  marginal_beta <- rnorm(n_basis, mean = 0, sd = 0.5)
  
  # Compute Gram matrix inverse
  V_inv <- solve(GramMatrix(base))
  
  # Test parameters
  alpha_vec <- c(-0.3, 0, 0.3, 0.6)
  tmin <- 10
  tmax <- 90
  
  # Imputation function for testing
  impute_fn <- function(t, data) {
    # Find the most recent observation
    recent_idx <- max(which(data$Time <= t))
    if (length(recent_idx) == 0 || recent_idx == 0) {
      # Use first observation if t is before all observations
      recent_idx <- 1
    }
    
    data.frame(
      Time = t,
      prev_outcome_lag1 = data$Outcome[recent_idx],
      delta_time_lag1 = t - data$Time[recent_idx]
    )
  }
  
  # Inverse link function (identity for linear model)
  inv_link <- function(x) x
  
  # Run comparison test
  comparison <- compare_term2_methods(
    patient_data = patient_data,
    outcome_model = outcome_model,
    base = base,
    alpha_vec = alpha_vec,
    marginal_beta = marginal_beta,
    V_inv = V_inv,
    tmin = tmin,
    tmax = tmax,
    impute_fn = impute_fn,
    inv_link = inv_link,
    tolerance = 1e-8
  )
  
  # Verify results match
  expect_true(comparison$summary$all_alphas_close,
              info = paste("Max diff:", comparison$summary$max_diff_overall))
  
  expect_true(comparison$summary$max_diff_overall < 1e-6,
              info = paste("Maximum difference too large:", comparison$summary$max_diff_overall))
  
  # Check that vectorized method is reasonably efficient
  expect_true(comparison$timing$speedup > 0.1,
              info = paste("Vectorized method should not be much slower. Speedup:", comparison$timing$speedup))
  
  # Verify structure of results
  expect_equal(length(comparison$comparisons), length(alpha_vec))
  expect_true(all(purrr::map_lgl(comparison$comparisons, ~ .x$vectorized_converged)))
})

test_that("vectorized integration matches original for GLM Poisson", {
  skip("Deprecated test - API changed to use variables/centering instead of impute_fn/inv_link. See test-vectorized-integration-simple.R for current tests.")
  skip_on_cran()
  
  # Generate synthetic data for Poisson GLM
  set.seed(123)
  n_obs <- 15
  patient_data <- data.frame(
    Time = sort(runif(n_obs, 0, 100)),
    Outcome = rpois(n_obs, lambda = 3),
    prev_outcome_lag1 = c(NA, rpois(n_obs-1, lambda = 2.8)),
    delta_time_lag1 = c(NA, diff(sort(runif(n_obs, 0, 100))))
  )
  
  # Remove first row
  patient_data <- patient_data[-1, ]
  
  # Fit Poisson GLM
  outcome_model <- glm(Outcome ~ prev_outcome_lag1 + delta_time_lag1, 
                       data = patient_data, family = poisson())
  
  # Create spline basis
  knots <- c(30, 70)
  base <- orthogonalsplinebasis::SplineBasis(knots = knots, order = 3L)
  
  n_basis <- ncol(base@basis)
  marginal_beta <- rnorm(n_basis, mean = 0, sd = 0.3)
  V_inv <- solve(GramMatrix(base))
  
  # Test with fewer alphas for speed
  alpha_vec <- c(-0.2, 0.2)
  tmin <- 15
  tmax <- 85
  
  impute_fn <- function(t, data) {
    recent_idx <- max(which(data$Time <= t))
    if (length(recent_idx) == 0 || recent_idx == 0) {
      recent_idx <- 1
    }
    
    data.frame(
      prev_outcome_lag1 = data$Outcome[recent_idx],
      delta_time_lag1 = pmax(0, t - data$Time[recent_idx])  # Ensure non-negative
    )
  }
  
  inv_link <- exp  # Log link inverse
  
  # Run comparison
  comparison <- compare_term2_methods(
    patient_data = patient_data,
    outcome_model = outcome_model,
    base = base,
    alpha_vec = alpha_vec,
    marginal_beta = marginal_beta,
    V_inv = V_inv,
    tmin = tmin,
    tmax = tmax,
    impute_fn = impute_fn,
    inv_link = inv_link,
    tolerance = 1e-8
  )
  
  # Verify results
  expect_true(comparison$summary$all_alphas_close,
              info = paste("Max diff:", comparison$summary$max_diff_overall))
  
  expect_true(comparison$summary$max_diff_overall < 1e-5)
  
  # Check convergence
  expect_true(all(purrr::map_lgl(comparison$comparisons, ~ .x$vectorized_converged)))
})

test_that("vectorized integration handles edge cases correctly", {
  skip("Deprecated test - API changed to use variables/centering instead of impute_fn/inv_link. See test-vectorized-integration-simple.R for current tests.")
  skip_on_cran()
  
  # Test with single alpha value
  set.seed(456)
  patient_data <- data.frame(
    Time = c(10, 20, 30),
    Outcome = c(5, 6, 7),
    prev_outcome_lag1 = c(4.5, 5.5, 6.5),
    delta_time_lag1 = c(10, 10, 10)
  )
  
  outcome_model <- lm(Outcome ~ prev_outcome_lag1, data = patient_data)
  
  # Simple basis
  base <- orthogonalsplinebasis::SplineBasis(knots = c(25), order = 2L)
  
  marginal_beta <- rep(0.1, ncol(base@basis))
  V_inv <- solve(GramMatrix(base))
  
  impute_fn <- function(t, data) {
    data.frame(prev_outcome_lag1 = data$Outcome[nrow(data)])
  }
  
  inv_link <- function(x) x
  
  # Test single alpha
  single_alpha_result <- compute_term2_influence_vectorized(
    patient_data = patient_data,
    outcome_model = outcome_model,
    base = base,
    alpha_vec = c(0),
    marginal_beta = marginal_beta,
    V_inv = V_inv,
    tmin = 5,
    tmax = 35,
    impute_fn = impute_fn,
    inv_link = inv_link
  )
  
  expect_equal(length(single_alpha_result), 1)
  expect_equal(single_alpha_result[[1]]$alpha, 0)
  expect_true(single_alpha_result[[1]]$converged)
  expect_equal(length(single_alpha_result[[1]]$Q), ncol(base@basis))
  
  # Test error handling
  expect_error(
    compute_term2_influence_vectorized(
      patient_data, outcome_model, base, numeric(0),
      marginal_beta, V_inv, 5, 35, impute_fn, inv_link
    ),
    "alpha_vec must contain at least one element"
  )
  
  expect_error(
    compute_term2_influence_vectorized(
      patient_data, outcome_model, base, c(0), 
      marginal_beta, V_inv, 35, 5, impute_fn, inv_link  # tmax < tmin
    ),
    "tmax must be greater than tmin"
  )
})

test_that("vectorized integration performance scales well with number of alphas", {
  skip("Deprecated test - API changed to use variables/centering instead of impute_fn/inv_link. See test-vectorized-integration-simple.R for current tests.")
  skip_on_cran()
  skip_if(Sys.getenv("SKIP_PERFORMANCE_TESTS") == "true", "Performance tests skipped")
  
  # Generate test data
  set.seed(789)
  patient_data <- data.frame(
    Time = sort(runif(10, 0, 50)),
    Outcome = rnorm(10, mean = 8, sd = 1.5),
    prev_outcome_lag1 = c(NA, rnorm(9, mean = 7.8, sd = 1.2)),
    delta_time_lag1 = c(NA, diff(sort(runif(10, 0, 50))))
  )
  patient_data <- patient_data[-1, ]
  
  outcome_model <- lm(Outcome ~ prev_outcome_lag1 + delta_time_lag1, data = patient_data)
  
  # Simple basis for speed
  base <- orthogonalsplinebasis::SplineBasis(knots = c(25), order = 2L)
  
  marginal_beta <- rnorm(ncol(base@basis), sd = 0.2)
  V_inv <- solve(GramMatrix(base))
  
  impute_fn <- function(t, data) {
    recent_idx <- max(which(data$Time <= t))
    if (length(recent_idx) == 0) recent_idx <- 1
    data.frame(
      prev_outcome_lag1 = data$Outcome[recent_idx],
      delta_time_lag1 = t - data$Time[recent_idx]
    )
  }
  
  inv_link <- function(x) x
  
  # Test with different numbers of alphas
  alpha_counts <- c(2, 5, 10)
  times <- numeric(length(alpha_counts))
  
  for (i in seq_along(alpha_counts)) {
    n_alpha <- alpha_counts[i]
    alpha_vec <- seq(-0.5, 0.5, length.out = n_alpha)
    
    time_taken <- system.time({
      result <- compute_term2_influence_vectorized(
        patient_data, outcome_model, base, alpha_vec,
        marginal_beta, V_inv, 5, 45, impute_fn, inv_link
      )
    })
    
    times[i] <- time_taken[["elapsed"]]
    
    # Verify all alphas converged
    expect_true(all(purrr::map_lgl(result, ~ .x$converged)))
  }
  
  # Performance should not degrade too severely with more alphas
  # (This is more of a smoke test than a strict requirement)
  expect_true(times[3] < times[1] * 10,  # 10x alphas shouldn't take >10x time
              info = paste("Times:", paste(times, collapse = ", ")))
  
  cat("Performance test - Alpha counts:", alpha_counts, "Times:", times, "\n")
})