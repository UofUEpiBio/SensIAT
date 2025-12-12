test_that("diagnostic: integrand values match at sample points", {
  skip("Diagnostic test - run manually when debugging. All integration issues resolved.")
  # This diagnostic test was used to identify the weight function bug (now fixed)
  # It's kept for historical reference but skipped in normal test runs
  skip_on_cran()
  
  # Same setup as main test
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

  df_i <- model.data |>
    dplyr::filter(Subject_ID == 1)

  variables <- list(
    time = rlang::sym("Time"),
    id = rlang::sym("Subject_ID"),
    outcome = rlang::sym("Outcome"),
    prev_time = rlang::sym("prev_time"),
    prev_outcome = rlang::sym("prev_outcome"),
    delta_time = rlang::sym("delta_time")
  )

  alpha_test <- 0
  n_basis <- dim(base@Matrices)[2]
  marginal_beta <- rep(0, n_basis)
  V_inv <- solve(orthogonalsplinebasis::GramMatrix(base))

  # Test at several time points
  test_times <- c(100, 200, 300, 400)
  
  cat("\n=== DIAGNOSTIC: Integrand Component Analysis ===\n")
  
  for (t in test_times) {
    cat(sprintf("\n--- Time = %.1f ---\n", t))
    
    # 1. Compute what ORIGINAL method uses as integrand
    imputed <- SensIAT:::impute_patient_df(t, df_i, variables, centering.statistics, TRUE)
    expected <- compute_SensIAT_expected_values(outcome.model, alpha_test, imputed)
    
    conditional_mean <- expected$E_Yexp_alphaY / expected$E_exp_alphaY
    B <- orthogonalsplinebasis::evaluate(base, t)
    original_integrand <- as.vector(B * conditional_mean)
    
    cat("Original integrand = B(t) * conditional_mean:\n")
    cat("  B(t):", as.vector(B), "\n")
    cat("  conditional_mean:", conditional_mean, "\n")
    cat("  E_Yexp_alphaY:", expected$E_Yexp_alphaY, "\n")
    cat("  E_exp_alphaY:", expected$E_exp_alphaY, "\n")
    cat("  integrand:", original_integrand, "\n")
    
    # 2. Test what the ACTUAL vectorized implementation produces now
    # Create the actual weight function as used in the implementation
    weight_fn_actual <- function(t) {
      as.vector(pcoriaccel_evaluate_basis(base, t))
    }
    marginal_mean_fn_actual <- function(t) { 0.0 }
    
    actual_weight <- weight_fn_actual(t)
    actual_marginal <- marginal_mean_fn_actual(t)
    actual_integrand <- actual_weight * (conditional_mean - actual_marginal)
    
    cat("\nACTUAL Vectorized implementation NOW:\n")
    cat("  weight = B(t):", actual_weight, "\n")
    cat("  marginal_mean:", actual_marginal, "\n")
    cat("  integrand = weight * conditional_mean:", actual_integrand, "\n")
    
    # Also compute OLD (wrong) version for comparison
    B_vec <- as.vector(B)
    mu <- sum(B_vec * marginal_beta)
    old_weight <- as.vector((V_inv %*% B_vec) * exp(-mu))
    marginal_mean <- sum(B_vec * marginal_beta)
    old_integrand <- old_weight * (conditional_mean - marginal_mean)
    
    cat("\nOLD (WRONG) Vectorized integrand:\n")
    cat("  weight = (V_inv %*% B) * exp(-mu):", old_weight, "\n")
    cat("  integrand:", old_integrand, "\n")
    
    cat("\nComparison:\n")
    cat("  Original:", original_integrand, "\n")
    cat("  NEW Vectorized:", actual_integrand, "\n")
    cat("  Ratio (should be ~1):", original_integrand / (actual_integrand + 1e-10), "\n")
  }
  
  cat("\n=== KEY FINDINGS ===\n")
  cat("1. Original uses: B(t) * E[Y|X(t)]\n")
  cat("2. Vectorized uses: (V_inv %*% B) * exp(-mu) * E[Y|X(t)]\n")
  cat("   (since marginal_mean = mu = 0 when marginal_beta = 0)\n")
  cat("3. The weight transformation changes magnitude dramatically\n")
  cat("\nHYPOTHESIS: Vectorized should use B(t) directly, not weighted!\n")
})
