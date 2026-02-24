# Debug: enable term2 tracing in main function
library(dplyr)
devtools::load_all("/workspaces/pcoriRPackage")

# Enable debug output for first call
assign("DEBUG_TERM2_ONCE", TRUE, envir = .GlobalEnv)

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

cat("\nFinal term2 for patient 1:", result$influence[[1]]$term2[[1]], "\n")
