# Simulate Data for a Single Subject

Internal function to simulate observation times and outcomes for one
subject, alternating between the two processes.

## Usage

``` r
simulate_single_subject(
  subject_id,
  End,
  intensity_coef,
  outcome_coef,
  baseline_hazard,
  outcome_sd,
  initial_outcome_mean,
  initial_outcome_sd,
  max_visits,
  link = "identity"
)
```

## Arguments

- subject_id:

  ID for this subject.

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

- link:

  Link function for outcome model. One of "identity", "log", or "logit".
  Determines the scale on which the outcome model operates.

## Value

A tibble with the subject's longitudinal data.
