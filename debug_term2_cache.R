# Debug to find why term2 is zero - with error capture
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

# Create base the same way main function does
spline.degree <- 3
full_knots <- c(
    rep(head(knots, 1), spline.degree),
    knots,
    rep(tail(knots, 1), spline.degree)
)
base <- orthogonalsplinebasis::SplineBasis(full_knots, order = spline.degree + 1L)
tmin <- base@knots[base@order]
tmax <- base@knots[length(base@knots) - base@order + 1]
cat("tmin:", tmin, "tmax:", tmax, "\n")

V_inv <- solve(orthogonalsplinebasis::GramMatrix(base))
W <- function(t, beta) {
  B <- pcoriaccel_evaluate_basis(base, t)
  as.numeric(B %*% V_inv %*% B)
}

# Set up parameters as main function does
marginal_beta <- rep(0, ncol(base))
alpha <- 0.5

# Check outcome model type
cat("\nOutcome model class:", class(outcome.model), "\n")
cat("Has kernel attr?:", !is.null(attr(outcome.model, "kernel")), "\n")
cat("Has bandwidth?:", !is.null(outcome.model$bandwidth), "\n")

# Patient data
patient_data <- data_with_lags %>% filter(Subject_ID == 1)
cat("\nPatient data time range:", range(patient_data$Time), "\n")
cat("Patient data observations in [tmin, tmax]:", sum(patient_data$Time >= tmin & patient_data$Time <= tmax), "\n")

# Test expected_get function behavior (simulating the cache)
cat("\n=== Testing expected_get at different time points ===\n")
# Test the impute_fn directly
test_times <- c(tmin, 200, 300, tmax)
for (t in test_times) {
  df <- impute_fn(t, patient_data)
  cat(sprintf("\nt=%f:\n", t))
  print(df)
  
  ev <- compute_SensIAT_expected_values(
    model = outcome.model,
    alpha = alpha,
    new.data = df
  )
  cat("E_exp_alphaY:", as.numeric(ev$E_exp_alphaY)[1], "\n")
  cat("E_Yexp_alphaY:", as.numeric(ev$E_Yexp_alphaY)[1], "\n")
}

# Test compute_term2_influence_fast with expected_get = NULL
cat("\n=== Direct call without cache ===\n")
result_no_cache <- SensIAT:::compute_term2_influence_fast(
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
cat("Result:", result_no_cache, "\n")

# Now test with cache (simulating what main function does)
cat("\n=== Testing with cache function ===\n")
cache_env <- new.env(parent = emptyenv())
expected_get <- function(t) {
  k <- as.character(signif(t, 12))
  if (!exists(k, cache_env, inherits = FALSE)) {
    df <- impute_fn(t, patient_data)
    ev <- compute_SensIAT_expected_values(
      model = outcome.model,
      alpha = alpha,
      new.data = df
    )
    cache_env[[k]] <- list(
      E_exp_alphaY = as.numeric(ev$E_exp_alphaY)[1],
      E_Yexp_alphaY = as.numeric(ev$E_Yexp_alphaY)[1]
    )
  }
  cache_env[[k]]
}

result_with_cache <- SensIAT:::compute_term2_influence_fast(
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
  expected_get = expected_get,
  time_var = "Time"
)
cat("Result with cache:", result_with_cache, "\n")

# Now compare with what main function actually returns
cat("\n=== Main function result ===\n")
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
  term2_method = "fast",
  use_expected_cache = FALSE  # Try without cache
)
cat("Term2 without cache:", result$influence[[1]]$term2[[1]], "\n")

result2 <- fit_SensIAT_marginal_mean_model_generalized(
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
  term2_method = "fast",
  use_expected_cache = TRUE  # With cache (default)
)
cat("Term2 with cache:", result2$influence[[1]]$term2[[1]], "\n")
