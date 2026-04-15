#!/usr/bin/env Rscript
# Quick test of the updated fast method with segmented integration

if (file.exists("DESCRIPTION")) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools is required to load the package from source. Please install it.")
  }
  devtools::load_all(quiet = TRUE)
} else {
  library(SensIAT)
}

library(dplyr)
library(splines)
library(survival)

cat("Setting up test data (first 3 patients)...\n")

# Setup test data
patient_ids <- unique(SensIAT_example_data$Subject_ID)[1:3]
data_with_lags <- SensIAT_example_data |>
  filter(Subject_ID %in% patient_ids) |>
  group_by(Subject_ID) |>
  mutate(
    ..prev_outcome.. = lag(Outcome, default = NA_real_, order_by = Time),
    ..prev_time.. = lag(Time, default = 0, order_by = Time),
    ..delta_time.. = Time - lag(.data$Time, default = NA_real_, order_by = Time)
  ) |>
  ungroup()

cat("Fitting outcome model...\n")

outcome.model <-
  fit_SensIAT_single_index_fixed_coef_model(
    Outcome ~ ns(..prev_outcome.., df=3) + ..delta_time.. - 1,
    id = Subject_ID,
    data = data_with_lags |> filter(Time > 0)
  )

cat("Setting up marginal mean parameters...\n")

alpha <- 0
knots <- c(60, 260, 460)
spline.degree <- 3L
knots_extended <- c(
  rep(head(knots, 1), spline.degree),
  knots,
  rep(tail(knots, 1), spline.degree)
)
base <- orthogonalsplinebasis::SplineBasis(knots_extended, order = spline.degree + 1L)
V.inv <- solve(orthogonalsplinebasis::GramMatrix(base))

tmin <- base@knots[base@order]
tmax <- base@knots[length(base@knots) - base@order + 1]

set.seed(123)
beta_test <- rep(0.5, ncol(base))

impute_data_fn <- function(t, df) {
  data_wl <- df |>
    mutate(
      ..prev_time.. = .data$Time,
      ..prev_outcome.. = .data$Outcome,
      ..delta_time.. = 0
    )
  extrapolate_from_last_observation(t, data_wl, 'Time', slopes = c('..delta_time..' = 1))
}

inv_link <- exp

cat("\nTesting FAST method on first patient...\n")

patient_id <- patient_ids[1]
patient_data <- data_with_lags[data_with_lags$Subject_ID == patient_id, ]

cat("Patient", patient_id, "has", nrow(patient_data), "observations\n")
cat("Observation times:", paste(patient_data$Time[!is.na(patient_data$Outcome)], collapse=", "), "\n")

tryCatch({
  result_fast <- SensIAT:::compute_term2_influence_fast(
    patient_data = patient_data,
    outcome_model = outcome.model,
    base = base,
    alpha = alpha,
    marginal_beta = beta_test,
    V_inv = V.inv,
    tmin = tmin,
    tmax = tmax,
    impute_fn = impute_data_fn,
    inv_link = inv_link
  )
  
  cat("✓ FAST method completed successfully\n")
  cat("Result length:", length(result_fast), "\n")
  cat("Result range: [", min(result_fast), ",", max(result_fast), "]\n")
  cat("All finite:", all(is.finite(result_fast)), "\n")
  
}, error = function(e) {
  cat("✗ FAST method failed with error:\n")
  cat(conditionMessage(e), "\n")
  traceback()
})

cat("\nTesting ORIGINAL method on same patient...\n")

tryCatch({
  result_original <- SensIAT:::compute_term2_influence_original(
    patient_data = patient_data,
    outcome_model = outcome.model,
    base = base,
    alpha = alpha,
    marginal_beta = beta_test,
    V_inv = V.inv,
    tmin = tmin,
    tmax = tmax,
    impute_fn = impute_data_fn,
    inv_link = inv_link
  )
  
  cat("✓ ORIGINAL method completed successfully\n")
  cat("Result length:", length(result_original), "\n")
  cat("Result range: [", min(result_original), ",", max(result_original), "]\n")
  
  if (exists("result_fast")) {
    cat("\n=== COMPARISON ===\n")
    cat("Max difference:", max(abs(result_fast - result_original)), "\n")
    if (max(abs(result_fast - result_original)) < 1e-4) {
      cat("✓ Methods match within tolerance\n")
    } else {
      cat("✗ Methods differ significantly\n")
      cat("Differences:", result_fast - result_original, "\n")
    }
  }
  
}, error = function(e) {
  cat("✗ ORIGINAL method failed with error:\n")
  cat(conditionMessage(e), "\n")
})

cat("\n=== TEST COMPLETE ===\n")
