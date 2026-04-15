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
  link = "identity"
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

## Value

A tibble with an additional `Treatment` column indicating group
assignment.

## Examples

``` r
# \donttest{
# Simulate data with treatment effect
sim_data <- simulate_SensIAT_two_groups(
    n_subjects = 100,
    End = 830,
    treatment_effect = 1.5,
    treatment_intensity_effect = 0.9
)
# }
```
