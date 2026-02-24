# Debug term2 with proper tracing
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

# Monkey-patch compute_term2_influence_fast to add debugging
original_fn <- SensIAT:::compute_term2_influence_fast
traced_fn <- function(...) {
  cat("compute_term2_influence_fast called\n")
  args <- list(...)
  cat("  patient_data nrow:", nrow(args$patient_data), "\n")
  cat("  patient_data Time range:", range(args$patient_data[[args$time_var %||% "Time"]]), "\n")
  cat("  tmin:", args$tmin, "tmax:", args$tmax, "\n")
  cat("  alpha:", args$alpha, "\n")
  cat("  ncol(base):", ncol(args$base), "\n")
  
  result <- tryCatch({
    original_fn(...)
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    traceback()
    rep(0, ncol(args$base))
  })
  
  cat("  result:", head(result), "\n")
  result
}

# Replace the function in the package namespace
assignInNamespace("compute_term2_influence_fast", traced_fn, "SensIAT")

# Now run main function
cat("\n=== Running main function with tracing ===\n")
result <- fit_SensIAT_marginal_mean_model_generalized(
  data = data_with_lags,
  time = data_with_lags$Time,
  id = data_with_lags$Subject_ID,
  alpha = 0.5,
  knots = knots,
  outcome.model = outcome.model,
  intensity.model = intensity.model,
  loss = "lp_mse",
  link = "identity",
  impute_data = impute_fn,
  term2_method = "fast",
  use_expected_cache = FALSE
)

cat("\n=== Final result ===\n")
cat("Coefficients:", result$coefficients[[1]], "\n")
cat("Term2 for patient 1:", result$influence[[1]]$term2[[1]], "\n")
