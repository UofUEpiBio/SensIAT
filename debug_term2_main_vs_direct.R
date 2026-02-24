# Debug to find why term2 is zero when called from main function
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

# Run main function
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
  term2_method = "fast"
)

cat("Coefficients:", result$coefficients[[1]], "\n")
cat("Term2 for patient 1:", result$influence[[1]]$term2[[1]], "\n")
cat("All term2 zero?:", all(sapply(result$influence[[1]]$term2, function(x) all(x == 0))), "\n")

# Now let's manually replicate what the main function does
spline.degree <- 3
full_knots <- c(
    rep(head(knots, 1), spline.degree),
    knots,
    rep(tail(knots, 1), spline.degree)
)
cat("\nFull knots:", full_knots, "\n")

base <- orthogonalsplinebasis::SplineBasis(full_knots, order = spline.degree + 1L)
tmin <- base@knots[base@order]
tmax <- base@knots[length(base@knots) - base@order + 1]
cat("tmin:", tmin, "tmax:", tmax, "\n")

V_inv <- solve(orthogonalsplinebasis::GramMatrix(base))

W <- function(t, beta) {
  B <- pcoriaccel_evaluate_basis(base, t)
  as.numeric(B %*% V_inv %*% B)
}

marginal_beta <- rep(0, ncol(base))
alpha <- 0.5

# Patient data
patient_data <- data_with_lags %>% filter(Subject_ID == 1)

# Direct call
cat("\n=== Direct call to compute_term2_influence_fast ===\n")
direct_result <- SensIAT:::compute_term2_influence_fast(
  patient_data = patient_data,
  outcome_model = outcome.model,
  base = base,
  alpha = alpha,
  marginal_beta = marginal_beta,
  V_inv = V_inv,
  tmin = tmin,
  tmax = tmax,
  impute_fn = impute_fn,
  inv_link = identity,
  W = W,
  expected_get = NULL,
  time_var = "Time"
)
cat("Direct result:", direct_result, "\n")
