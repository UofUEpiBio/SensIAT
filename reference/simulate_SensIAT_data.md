# Simulate SensIAT Data

Generates simulated longitudinal data following the SensIAT model
structure, alternating between observation time generation (Cox
proportional hazards model) and outcome generation.

## Usage

``` r
simulate_SensIAT_data(
  n_subjects,
  End,
  intensity_coef = -0.05,
  outcome_coef = list(prev_outcome = c(0.7, -0.1, 0.05), time = -0.001, delta_time =
    -0.002, intercept = 2),
  baseline_hazard = 0.005,
  outcome_sd = 1.5,
  baseline_outcome_fn = NULL,
  initial_outcome_mean = 5,
  initial_outcome_sd = 2,
  max_visits = 50,
  seed = NULL,
  link = "identity",
  intensity_fn = NULL,
  intensity_bound = NULL,
  outcome_model = NULL,
  outcome_simulator = NULL
)
```

## Arguments

- n_subjects:

  Number of subjects to simulate.

- End:

  Maximum follow-up time.

- intensity_coef:

  Coefficient for the effect of previous outcome on observation
  intensity. Can be a scalar or vector (one per visit number stratum).

- outcome_coef:

  Named list of coefficients for outcome model including:

  - `prev_outcome` - coefficients for natural spline of previous outcome
    (length 3)

  - `time` - coefficient for time since baseline

  - `delta_time` - coefficient for time since last observation

  - `intercept` - intercept term

- baseline_hazard:

  Baseline hazard function. Either a function of time and visit number,
  or a numeric value for constant baseline hazard.

- outcome_sd:

  Standard deviation of the outcome residuals.

- baseline_outcome_fn:

  Optional function to generate baseline outcome value for each subject.

- initial_outcome_mean:

  Mean of the initial (baseline) outcome.

- initial_outcome_sd:

  Standard deviation of the initial outcome.

- max_visits:

  Maximum number of visits per subject (to prevent infinite loops).

- seed:

  Random seed for reproducibility.

- link:

  Link function for outcome model. One of "identity", "log", or "logit".
  Determines the scale on which the outcome model operates.

- intensity_fn:

  Optional function to compute intensity (hazard) of observation. If
  provided, should take arguments (`time`, `prev_outcome`, `visit_num`)
  and return a scalar intensity value. If `NULL` (default), intensity is
  computed from `intensity_coef` and `baseline_hazard`.

- intensity_bound:

  Upper bound on intensity for rejection sampling. Required if
  intensity_fn is provided. Represents the supremum of the intensity
  function on the interval of interest.

- outcome_model:

  Optional fitted single-index outcome model. If provided, outcomes for
  follow-up visits are generated from the fitted model via
  [`make_single_index_simulator()`](https://uofuepibio.github.io/SensIAT/reference/make_single_index_simulator.md).

- outcome_simulator:

  Optional simulator function for follow-up outcomes. When provided, it
  overrides the internal outcome generation function. This function
  should accept `prev_outcome`, `time`, `delta_time`, and optionally
  `newdata`.

## Value

A tibble with columns:

- `Subject_ID` - Subject identifier

- `Time` - Observation time

- `Outcome` - Observed outcome value

- Additional columns may be added for internal use

## Examples

``` r
# \donttest{
# Default usage (uses exponential gaps derived from intensity_coef and baseline_hazard)
sim_data <- simulate_SensIAT_data(
    n_subjects = 100,
    End = 830,
    intensity_coef = -0.05,
    outcome_coef = list(
        prev_outcome = c(0.7, -0.1, 0.05),
        time = -0.001,
        delta_time = -0.002,
        intercept = 2
    ),
    baseline_hazard = 0.005,
    outcome_sd = 1.5
)

# Example with custom intensity function and thinning (Exp(lambda_star))
intensity_fn <- function(t, prev_outcome, visit_num) {
  lambda0 <- 0.005
  gamma  <- -0.05
  lambda0 * exp(gamma * prev_outcome)
}
sim_data2 <- simulate_SensIAT_data(
    n_subjects = 50,
    End = 200,
    seed = 123,
    intensity_fn = intensity_fn,
    intensity_bound = 0.05,
    max_visits = 20
)

# Example using a fitted single-index outcome model to generate outcomes
# via the fitted conditional distribution.
#
# Note: this example uses a fitted model object and is intended for
# parametric bootstrap-style simulation.
#
# 
# 
# outcome_model <- fit_SensIAT_single_index_fixed_coef_model(
#     Outcome ~ prev_outcome + time + delta_time,
#     data = training_data,
#     id = Subject_ID
# )
# sim_data3 <- simulate_SensIAT_data(
#     n_subjects = 50,
#     End = 200,
#     seed = 123,
#     outcome_model = outcome_model,
#     intensity_fn = intensity_fn,
#     intensity_bound = 0.05,
#     max_visits = 20
# )
# }
```
