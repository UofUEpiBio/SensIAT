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
  max_visits,
  baseline_outcome_fn = NULL,
  initial_outcome_mean,
  initial_outcome_sd,
  link = "identity",
  intensity_fn = NULL,
  intensity_bound = NULL,
  outcome_simulator = NULL
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

- max_visits:

  Maximum number of visits per subject (to prevent infinite loops).

- baseline_outcome_fn:

  Optional function to generate baseline outcome value for each subject.

- initial_outcome_mean:

  Mean of the initial (baseline) outcome.

- initial_outcome_sd:

  Standard deviation of the initial outcome.

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

- outcome_simulator:

  Optional simulator function for follow-up outcomes. When provided, it
  overrides the internal outcome generation function. This function
  should accept `prev_outcome`, `time`, `delta_time`, and optionally
  `newdata`.

## Value

A tibble with the subject's longitudinal data.
