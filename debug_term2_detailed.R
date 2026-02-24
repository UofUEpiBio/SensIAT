# Debug script for term2 zero values - detailed investigation
library(dplyr)
devtools::load_all("/workspaces/pcoriRPackage")

set.seed(456)
data_with_lags <- SensIAT_example_data %>%
  group_by(Subject_ID) %>%
  mutate(
    ..prev_outcome.. = lag(Outcome, default = NA_real_, order_by = Time),
    ..prev_time.. = lag(Time, default = 0, order_by = Time),
    ..delta_time.. = Time - lag(Time, default = NA_real_, order_by = Time)
  ) %>%
  ungroup() %>%
  filter(Subject_ID <= 5)

intensity.model <- survival::coxph(
  Surv(..prev_time.., Time, !is.na(Outcome)) ~ ..prev_outcome.. + strata(Visit),
  data = data_with_lags %>% filter(Time > 0)
)

outcome.model <- fit_SensIAT_single_index_fixed_coef_model(
  Outcome ~ splines::ns(..prev_outcome.., df = 2) + ..delta_time.. - 1,
  id = Subject_ID,
  data = data_with_lags %>% filter(Time > 0)
)

impute_fn <- function(t, df) {
  data_wl <- df %>%
    mutate(
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
alpha <- 0.5

# Test expected values computation at different time points for patient 1
patient_data <- data_with_lags %>% filter(Subject_ID == 1)
cat("=== Patient 1 data ===\n")
print(patient_data)

cat("\n=== Testing imputation and expected values ===\n")
test_times <- c(50, 100, 150, 200, 300, 400, 500)
for (t in test_times) {
  imputed <- impute_fn(t, patient_data)
  cat(sprintf("\nt=%d imputed data:\n", t))
  print(imputed)
  
  expected <- compute_SensIAT_expected_values(
    model = outcome.model,
    alpha = alpha,
    new.data = imputed
  )
  cat(sprintf("Expected: E_exp_alphaY=%.6f, E_Yexp_alphaY=%.6f, ratio=%.6f\n", 
              expected$E_exp_alphaY, expected$E_Yexp_alphaY, 
              expected$E_Yexp_alphaY / expected$E_exp_alphaY))
}

# Now run full fit and check term2
cat("\n=== Running full fit ===\n")
result <- fit_SensIAT_marginal_mean_model_generalized(
  data = data_with_lags,
  time = data_with_lags$Time,
  id = data_with_lags$Subject_ID,
  alpha = alpha,
  knots = knots,
  outcome.model = outcome.model,
  intensity.model = intensity.model,
  loss = "lp_mse",
  link = "identity",
  impute_data = impute_fn,
  term2_method = "fast"
)

cat("Coefficients:\n")
print(result$coefficients[[1]])

cat("\nTerm2 for all patients:\n")
for (i in seq_along(result$influence[[1]]$term2)) {
  cat(sprintf("  Patient %d: ", i))
  cat(paste(result$influence[[1]]$term2[[i]], collapse = " "), "\n")
}

cat("\nTerm1 for patient 1 (first 3):\n")
for (i in 1:min(3, length(result$influence[[1]]$term1))) {
  cat(sprintf("  Obs %d: ", i))
  cat(paste(round(result$influence[[1]]$term1[[i]], 4), collapse = " "), "\n")
}
