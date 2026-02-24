library(dplyr)
library(SensIAT)

set.seed(123)
n_patients <- 2

sim_data <- lapply(1:n_patients, function(id) {
  n_visits <- 4
  times <- sort(c(0, cumsum(rexp(n_visits - 1, rate = 1/50))))
  outcomes <- numeric(n_visits)
  outcomes[1] <- rnorm(1, mean = 10, sd = 2)
  for (i in 2:n_visits) {
    delta_t <- times[i] - times[i-1]
    outcomes[i] <- 0.7 * outcomes[i-1] + 0.05 * delta_t + rnorm(1, sd = 1.5)
  }
  data.frame(
    Subject_ID = id,
    Time = times,
    Outcome = outcomes,
    Visit = 1:n_visits
  )
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

# Check what the outcome model looks like
cat("Outcome model formula:\n")
print(formula(outcome.model))
cat("\nOutcome model data columns:\n")
print(names(outcome.model$data))
cat("\nOutcome model coefficients:\n")
print(outcome.model$coef)

# Check what happens when we compute expected values
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

# Get the first patient's data
patient_1_data <- data_with_lags |> filter(Subject_ID == 1)
cat("\n\nPatient 1 data:\n")
print(patient_1_data)

# Impute at a time point
imputed_100 <- impute_fn(100, patient_1_data)
cat("\n\nImputed data at t=100:\n")
print(imputed_100)

# Compute expected values using the outcome model
expected_100 <- compute_SensIAT_expected_values(
  model = outcome.model,
  alpha = 0,
  new.data = imputed_100
)
cat("\n\nExpected values at t=100:\n")
print(expected_100)

# Print the expected values with more digits
cat("\n\nDetailed expected values:\n")
cat("E_Yexp_alphaY:", expected_100$E_Yexp_alphaY, "\n")
cat("E_exp_alphaY:", expected_100$E_exp_alphaY, "\n")
cat("Ratio:", expected_100$E_Yexp_alphaY / expected_100$E_exp_alphaY, "\n")
