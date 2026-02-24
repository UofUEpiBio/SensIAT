# Debug: Check what expected values are returning
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
  filter(Subject_ID == 1)  # Test single patient

outcome.model <- fit_SensIAT_single_index_fixed_coef_model(
  Outcome ~ splines::ns(..prev_outcome.., df = 2) + ..delta_time.. - 1,
  id = Subject_ID,
  data = data_with_lags %>% filter(Time > 0)
)

cat("Model class:", class(outcome.model), "\n")
cat("Is single-index:", is(outcome.model, "SensIAT::Single-index-outcome-model"), "\n")
cat("Has kernel attr:", !is.null(attr(outcome.model, "kernel")), "\n")
cat("Has bandwidth:", !is.null(outcome.model$bandwidth), "\n")

# Test compute_expected_values_at_grid
grid <- c(100, 200, 300, 400, 500)
cat("\nTesting compute_expected_values_at_grid:\n")
E_grid <- SensIAT:::compute_expected_values_at_grid(
  grid = grid,
  patient_data = data_with_lags,
  outcome_model = outcome.model,
  alpha = 0,
  time_var = "Time"
)

cat("\nE_grid results:\n")
for (i in seq_along(grid)) {
  cat(sprintf("t=%.0f: E_exp=%.4f, E_Yexp=%.4f\n", 
              grid[i], E_grid[[i]]$E_exp_alphaY, E_grid[[i]]$E_Yexp_alphaY))
}

# Also check what fast method uses for expected values
cat("\n=== What fast method computes ===\n")
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

for (t in grid) {
  df_at_t <- impute_fn(t, data_with_lags)
  ev <- compute_SensIAT_expected_values(model = outcome.model, alpha = 0, new.data = df_at_t)
  cat(sprintf("t=%.0f: E_exp=%.4f, E_Yexp=%.4f\n", 
              t, ev$E_exp_alphaY, ev$E_Yexp_alphaY))
}
