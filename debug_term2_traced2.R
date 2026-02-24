# Debug term2 with in-function tracing
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

# Create a traced version that shows actual execution
traced_compute_term2 <- function(
  patient_data,
  outcome_model,
  base,
  alpha,
  marginal_beta,
  V_inv,
  tmin,
  tmax,
  impute_fn,
  inv_link,
  W,
  expected_get = NULL,
  time_var = NULL,
  ...
) {
  cat("TRACED: compute_term2_influence_fast called\n")
  cat("  time_var:", deparse(time_var), "\n")
  cat("  class(time_var):", class(time_var), "\n")
  cat("  patient_data ncol:", ncol(patient_data), "nrow:", nrow(patient_data), "\n")
  cat("  colnames:", paste(names(patient_data), collapse=", "), "\n")
  cat("  tmin:", tmin, "tmax:", tmax, "\n")
  
  # Test what time_var access gives
  if (!is.null(time_var)) {
    if (is.character(time_var)) {
      cat("  time_var is character\n")
      cat("  patient_data[[time_var]]:", head(patient_data[[time_var]]), "\n")
    } else if (rlang::is_quosure(time_var)) {
      cat("  time_var is quosure\n")
      cat("  quo_name:", tryCatch(rlang::quo_name(time_var), error = function(e) paste("ERROR:", e$message)), "\n")
    } else {
      cat("  time_var is other type:", class(time_var), "\n")
    }
  }
  
  # Now call the original function
  result <- SensIAT:::compute_term2_influence_fast(
    patient_data = patient_data,
    outcome_model = outcome_model,
    base = base,
    alpha = alpha,
    marginal_beta = marginal_beta,
    V_inv = V_inv,
    tmin = tmin,
    tmax = tmax,
    impute_fn = impute_fn,
    inv_link = inv_link,
    W = W,
    expected_get = expected_get,
    time_var = time_var
  )
  
  cat("  result:", head(result), "\n\n")
  result
}

# Temporarily replace the function
assignInNamespace("compute_term2_influence_fast", traced_compute_term2, "SensIAT")

# Now run main function
cat("=== Running main function with tracing ===\n\n")
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
cat("Term2 for patient 1:", result$influence[[1]]$term2[[1]], "\n")
