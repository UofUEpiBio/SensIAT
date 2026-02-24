# Quick test of all four term2 methods
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
  filter(Subject_ID <= 3)  # Small subset for speed

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

# Test all four methods
methods <- c("fast", "original", "fixed_grid", "seeded_adaptive")

results <- list()
for (m in methods) {
  cat("\n=== Testing method:", m, "===\n")
  t1 <- Sys.time()
  result <- fit_SensIAT_marginal_mean_model_generalized(
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
    term2_method = m,
    term2_grid_n = if(m %in% c("fixed_grid", "seeded_adaptive")) 200 else 100
  )
  t2 <- Sys.time()
  cat("Time:", round(difftime(t2, t1, units="secs"), 2), "secs\n")
  cat("Coefficients:", round(result$coefficients[[1]], 6), "\n")
  results[[m]] <- result$coefficients[[1]]
}

cat("\n\n=== SUMMARY ===\n")
cat("Method            Coefficients\n")
cat(sprintf("%-17s %s\n", "fast:", paste(round(results$fast, 6), collapse=" ")))
cat(sprintf("%-17s %s\n", "original:", paste(round(results$original, 6), collapse=" ")))
cat(sprintf("%-17s %s\n", "fixed_grid:", paste(round(results$fixed_grid, 6), collapse=" ")))
cat(sprintf("%-17s %s\n", "seeded_adaptive:", paste(round(results$seeded_adaptive, 6), collapse=" ")))

cat("\nMax absolute differences:\n")
cat("fast vs original:", max(abs(results$fast - results$original)), "\n")
cat("fast vs fixed_grid:", max(abs(results$fast - results$fixed_grid)), "\n")
cat("fast vs seeded_adaptive:", max(abs(results$fast - results$seeded_adaptive)), "\n")
