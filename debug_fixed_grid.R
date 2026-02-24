# Debug: Compare fixed_grid vs fast at term2 level
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

cat("Patient 1 data:\n")
print(data_with_lags %>% select(Subject_ID, Time, Outcome, ..prev_outcome.., ..delta_time..))

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
spline.degree <- 3L
full_knots <- c(rep(knots[1], spline.degree), knots, rep(knots[length(knots)], spline.degree))
base <- orthogonalsplinebasis::SplineBasis(full_knots, order = spline.degree + 1L)
V.inv <- solve(orthogonalsplinebasis::GramMatrix(base))

# Setup parameters
tmin <- 100
tmax <- 500
alpha <- 0
marginal_beta <- rep(0, ncol(base))
inv_link <- function(x) x  # identity link
W <- function(t, beta) SensIAT:::pcoriaccel_evaluate_basis(base, t) %*% V.inv

patient_data <- data_with_lags

# Call fast method directly
cat("\n=== FAST method ===\n")
term2_fast <- SensIAT:::compute_term2_influence_fast(
  patient_data = patient_data,
  outcome_model = outcome.model,
  base = base,
  alpha = alpha,
  marginal_beta = marginal_beta,
  V_inv = V.inv,
  tmin = tmin,
  tmax = tmax,
  impute_fn = impute_fn,
  inv_link = inv_link,
  W = W,
  time_var = "Time"
)
cat("term2_fast:", term2_fast, "\n")

# Call fixed_grid method directly
cat("\n=== FIXED_GRID method (n=200) ===\n")
term2_fixed <- SensIAT:::compute_term2_influence_fixed_grid(
  patient_data = patient_data,
  outcome_model = outcome.model,
  base = base,
  alpha = alpha,
  marginal_beta = marginal_beta,
  V_inv = V.inv,
  tmin = tmin,
  tmax = tmax,
  impute_fn = impute_fn,
  inv_link = inv_link,
  W = W,
  n_grid = 200,
  time_var = "Time"
)
cat("term2_fixed:", term2_fixed, "\n")

cat("\n=== Comparison ===\n")
cat("Difference:", term2_fast - term2_fixed, "\n")
cat("Ratio:", term2_fast / term2_fixed, "\n")
