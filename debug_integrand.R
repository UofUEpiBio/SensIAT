library(dplyr)
library(SensIAT)

set.seed(123)
n_patients <- 1

sim_data <- lapply(1:n_patients, function(id) {
  n_visits <- 4
  times <- c(0, 42.2, 71.0, 137)
  outcomes <- c(12.4, 8.24, 9.07, 9.51)
  data.frame(Subject_ID = id, Time = times, Outcome = outcomes, Visit = 1:n_visits)
}) |> bind_rows()

data_with_lags <- sim_data |>
  group_by(Subject_ID) |>
  mutate(
    ..prev_outcome.. = lag(Outcome, default = NA_real_, order_by = Time),
    ..prev_time.. = lag(Time, default = 0, order_by = Time),
    ..delta_time.. = Time - lag(Time, default = NA_real_, order_by = Time)
  ) |>
  ungroup()

intensity.model <- survival::coxph(
  Surv(..prev_time.., Time, !is.na(Outcome)) ~ ..prev_outcome.. + strata(Visit),
  data = data_with_lags |> filter(Time > 0)
)

outcome.model <- fit_SensIAT_single_index_fixed_coef_model(
  Outcome ~ ..prev_outcome.. + Time + ..delta_time.. - 1,
  id = Subject_ID,
  data = data_with_lags |> filter(Time > 0),
  kernel = "dnorm",
  abs.tol = 1e-7
)

impute_fn <- function(t, df) {
  data_wl <- df |>
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

knots <- c(50, 150, 250)
# Create a simple basis for testing
base <- list()  # We'll just skip the basis evaluation for now

patient_data <- data_with_lags |> filter(Subject_ID == 1)
tmin <- min(data_with_lags$Time)
tmax <- max(data_with_lags$Time)

# Compute W function
beta <- c(1, 1, 1)  # Dummy coefficients
V_inv <- diag(3)
inv_link <- function(x) x
W <- function(t, marginal_beta) {
  rep(1, ncol(V_inv))  # Dummy weight for testing
}

# Compare the integrand at a few points
cat("Comparing integrand evaluations:\n")
for (t in c(50, 75, 100)) {
  cat("\n=== At t =", t, "===\n")
  
  # Get imputed data
  imputed <- impute_fn(t, patient_data)
  cat("Imputed data:\n")
  print(imputed)
  
  # Compute expected values
  expected <- compute_SensIAT_expected_values(
    model = outcome.model,
    alpha = 0,
    new.data = imputed
  )
  cat("Expected values:\n")
  cat("  E_Yexp_alphaY:", expected$E_Yexp_alphaY, "\n")
  cat("  E_exp_alphaY:", expected$E_exp_alphaY, "\n")
  cat("  Ratio:", expected$E_Yexp_alphaY / expected$E_exp_alphaY, "\n")
}
