# Simulate Treatment and Control Groups

Generate simulated data for both treatment and control groups with
potentially different parameters.

## Usage

``` r
simulate_SensIAT_two_groups(
  n_subjects,
  End,
  intensity_coef = -0.05,
  outcome_coef = list(prev_outcome = c(0.7, -0.1, 0.05), time = -0.001, delta_time =
    -0.002, intercept = 2),
  baseline_hazard = 0.005,
  outcome_sd = 1.5,
  initial_outcome_mean = 5,
  initial_outcome_sd = 2,
  max_visits = 50,
  treatment_effect = 0,
  treatment_intensity_effect = 1,
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

- initial_outcome_mean:

  Mean of the initial (baseline) outcome.

- initial_outcome_sd:

  Standard deviation of the initial outcome.

- max_visits:

  Maximum number of visits per subject (to prevent infinite loops).

- treatment_effect:

  Additive treatment effect on outcomes (added to intercept on link
  scale).

- treatment_intensity_effect:

  Multiplicative effect on observation intensity (values \< 1 mean fewer
  observations in treatment group).

- seed:

  Random seed for reproducibility.

- link:

  Link function for outcome model. One of "identity", "log", or "logit".

- intensity_fn:

  Optional function to compute intensity (hazard) of observation. If
  provided, should take arguments (`time`, `prev_outcome`, `visit_num`)
  and return a scalar intensity value. If `NULL` (default), intensity is
  computed from `intensity_coef` and `baseline_hazard`.

- intensity_bound:

  Upper bound on intensity for rejection sampling. Required if
  `intensity_fn` is provided. Represents the supremum of the intensity
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

A tibble with an additional `Treatment` column indicating group
assignment.

## Examples

``` r
# \donttest{
# Default treatment/control simulation (uses exponential gaps derived from coefficients)
sim_data <- simulate_SensIAT_two_groups(
    n_subjects = 100,
    End = 830,
    treatment_effect = 1.5,
    treatment_intensity_effect = 0.9
)

# Example using custom intensity with thinning
intensity_fn <- function(t, prev_outcome, visit_num) {
  lambda0 <- 0.005
  gamma  <- -0.05
  lambda0 * exp(gamma * prev_outcome)
}
sim_data2 <- simulate_SensIAT_two_groups(
    n_subjects = 100,
    End = 200,
    seed = 123,
    intensity_fn = intensity_fn,
    intensity_bound = 0.05,
    treatment_effect = 1.0,
    treatment_intensity_effect = 0.9
)
# }
```
