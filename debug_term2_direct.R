# Debug script for term2 zero values - direct call to term2 function
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

knots <- c(150, 300, 450)  # interior knots between tmin and tmax
alpha <- 0.5
spline.degree <- 3

# Set up the same way the main function does
tmin <- min(data_with_lags$Time[data_with_lags$Time > 0])
tmax <- max(data_with_lags$Time)
cat("tmin:", tmin, "tmax:", tmax, "\n")

# Make sure knots are within [tmin, tmax]
knots <- knots[knots > tmin & knots < tmax]
cat("Effective knots:", knots, "\n")

# Pad knots the same way the main function does
full_knots <- c(
    rep(head(knots, 1), spline.degree),
    knots,
    rep(tail(knots, 1), spline.degree)
)
cat("Full knots:", full_knots, "\n")
base <- orthogonalsplinebasis::SplineBasis(full_knots, order = spline.degree + 1L)
V_inv <- solve(orthogonalsplinebasis::GramMatrix(base))
cat("Number of basis functions:", ncol(base), "\n")

# Weight function W
W <- function(t, beta) {
  B <- pcoriaccel_evaluate_basis(base, t)
  as.numeric(B %*% V_inv %*% B)
}

# Initial beta
marginal_beta <- rep(0, ncol(base))

# Patient data
patient_data <- data_with_lags %>% filter(Subject_ID == 1)
cat("\nPatient 1 data:\n")
print(patient_data)

# Call compute_term2_influence_fast directly
cat("\n=== Calling compute_term2_influence_fast directly ===\n")
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
  cat("Error:", conditionMessage(e), "\n")
  cat("Traceback:\n")
  traceback()
  NULL
})

cat("\nResult:\n")
print(result)

# Now let's manually compute the integrand at a few points
cat("\n=== Manual integrand computation ===\n")
test_times <- c(tmin, (tmin + tmax)/2, tmax)
for (t in test_times) {
  imputed <- impute_fn(t, patient_data)
  expected <- compute_SensIAT_expected_values(
    model = outcome.model,
    alpha = alpha,
    new.data = imputed
  )
  
  B <- pcoriaccel_evaluate_basis(base, t)
  weight <- W(t, marginal_beta)
  
  term <- expected$E_Yexp_alphaY / expected$E_exp_alphaY - crossprod(B, marginal_beta)
  integrand_value <- weight * as.numeric(term)
  
  cat(sprintf("\nt=%f:\n", t))
  cat("  E_exp_alphaY:", expected$E_exp_alphaY, "\n")
  cat("  E_Yexp_alphaY:", expected$E_Yexp_alphaY, "\n")
  cat("  ratio:", expected$E_Yexp_alphaY / expected$E_exp_alphaY, "\n")
  cat("  B'beta:", as.numeric(crossprod(B, marginal_beta)), "\n")
  cat("  term:", as.numeric(term), "\n")
  cat("  weight:", weight, "\n")
  cat("  integrand:", integrand_value, "\n")
}
