# Debug term2 - check observation times in integration range
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

knots <- c(100, 300, 500)
spline.degree <- 3
full_knots <- c(
    rep(head(knots, 1), spline.degree),
    knots,
    rep(tail(knots, 1), spline.degree)
)

base <- orthogonalsplinebasis::SplineBasis(full_knots, order = spline.degree + 1L)
tmin <- base@knots[base@order]
tmax <- base@knots[length(base@knots) - base@order + 1]

cat("Integration range: [", tmin, ", ", tmax, "]\n\n")

# Check each patient's observation times
for (pid in 1:5) {
  patient_data <- data_with_lags %>% filter(Subject_ID == pid)
  times <- patient_data$Time
  in_range <- times[times >= tmin & times <= tmax]
  cat("Patient", pid, "times:", paste(times, collapse=", "), "\n")
  cat("  in range [", tmin, ",", tmax, "]:", paste(in_range, collapse=", "), "\n")
  cat("  count:", length(in_range), "\n\n")
}

# The problem: if there are no observation times in range,
# the breakpoints will only be [tmin, tmax] which should still work
# Let's test with patient 1 directly

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

V_inv <- solve(orthogonalsplinebasis::GramMatrix(base))
W <- function(t, beta) {
  B <- pcoriaccel_evaluate_basis(base, t)
  as.numeric(B %*% V_inv %*% B)
}

patient_data <- data_with_lags %>% filter(Subject_ID == 1)
marginal_beta <- rep(0, ncol(base))
alpha <- 0.5

cat("\n=== Direct call to compute_term2_influence_fast ===\n")
result <- tryCatch({
  SensIAT:::compute_term2_influence_fast(
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
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  traceback()
  NULL
})
cat("Result:", result, "\n")
