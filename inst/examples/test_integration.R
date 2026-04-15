#!/usr/bin/env Rscript
# Test script to verify term2_method integration in fit_SensIAT_marginal_mean_model_generalized

# Load package from source when running inside the repo; otherwise try installed package
if (file.exists("DESCRIPTION")) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools is required to load the package from source. Please install it.")
  }
  devtools::load_all(quiet = TRUE)
} else {
  library(SensIAT)
}

library(dplyr)
library(survival)
library(splines)  # for ns() function

cat("Setting up test data...\n")

# Setup test data (subset to first 10 patients to keep it fast but stable)
patient_ids <- unique(SensIAT_example_data$Subject_ID)[1:3]
data_with_lags <- SensIAT_example_data |>
  dplyr::filter(Subject_ID %in% patient_ids) |>
  group_by(Subject_ID) |>
  mutate(
    ..prev_outcome.. = lag(Outcome, default = NA_real_, order_by = Time),
    ..prev_time.. = lag(Time, default = 0, order_by = Time),
    ..delta_time.. = Time - lag(.data$Time, default = NA_real_, order_by = Time)
  )

# Create models
intensity.model <-
  coxph(Surv(..prev_time.., Time, !is.na(Outcome)) ~ ..prev_outcome.. + strata(Visit),
        data = data_with_lags |> filter(.data$Time > 0))

outcome.model <-
  fit_SensIAT_single_index_fixed_coef_model(
    Outcome ~ ns(..prev_outcome.., df=3) + ..delta_time.. - 1,
    id = Subject_ID,
    data = data_with_lags |> filter(Time > 0))

impute_data <- function(t, df) {
  data_wl <- df |>
    mutate(..prev_time.. = Time,
           ..prev_outcome.. = Outcome,
           ..delta_time.. = 0)
  extrapolate_from_last_observation(t, data_wl, 'Time', slopes = c('..delta_time..' = 1))
}

cat("\nTesting FAST method (default)...\n")
time_fast <- system.time({
  fit_fast <- fit_SensIAT_marginal_mean_model_generalized(
    data = data_with_lags,
    time = data_with_lags$Time,
    id = data_with_lags$Subject_ID,
    alpha = 0,
    knots = c(60, 260, 460),
    outcome.model = outcome.model,
    intensity.model = intensity.model,
    loss = 'lp_mse',
    link = 'log',
    term2_method = 'fast',
    impute_data = impute_data
  )
})

cat("FAST method completed in", time_fast["elapsed"], "seconds\n")
cat("Convergence:", fit_fast$convergence, "\n")
cat("Coefficients:", fit_fast$par, "\n")

cat("\nTesting ORIGINAL method...\n")
time_original <- system.time({
  fit_original <- fit_SensIAT_marginal_mean_model_generalized(
    data = data_with_lags,
    time = data_with_lags$Time,
    id = data_with_lags$Subject_ID,
    alpha = 0,
    knots = c(60, 260, 460),
    outcome.model = outcome.model,
    intensity.model = intensity.model,
    loss = 'lp_mse',
    link = 'log',
    term2_method = 'original',
    impute_data = impute_data
  )
})

cat("ORIGINAL method completed in", time_original["elapsed"], "seconds\n")
cat("Convergence:", fit_original$convergence, "\n")
cat("Coefficients:", fit_original$par, "\n")

# Compare results
cat("\n=== COMPARISON ===\n")
cat("Speedup:", time_original["elapsed"] / time_fast["elapsed"], "x\n")
cat("Max coefficient difference:", max(abs(fit_fast$par - fit_original$par)), "\n")
cat("Convergence match:", fit_fast$convergence == fit_original$convergence, "\n")

if (max(abs(fit_fast$par - fit_original$par)) < 1e-4) {
  cat("\n✓ PASS: Methods produce nearly identical results\n")
} else {
  cat("\n✗ FAIL: Methods differ significantly\n")
}

if (time_fast["elapsed"] < time_original["elapsed"]) {
  cat("✓ PASS: FAST method is faster\n")
} else {
  cat("✗ FAIL: FAST method is not faster\n")
}
