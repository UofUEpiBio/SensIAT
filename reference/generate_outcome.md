# Generate Outcome Value

Generate an outcome value based on previous outcome, time since
baseline, and time since last observation using a non-linear model.

## Usage

``` r
generate_outcome(
  prev_outcome,
  time,
  delta_time,
  outcome_coef,
  outcome_sd,
  link = "identity"
)
```

## Arguments

- prev_outcome:

  Previous outcome value.

- time:

  Current time (time since baseline).

- delta_time:

  Time since last observation.

- outcome_coef:

  Named list of coefficients for outcome model including:

  - `prev_outcome` - coefficients for natural spline of previous outcome
    (length 3)

  - `time` - coefficient for time since baseline

  - `delta_time` - coefficient for time since last observation

  - `intercept` - intercept term

- outcome_sd:

  Standard deviation of the outcome residuals.

- link:

  Link function for outcome model. One of "identity", "log", or "logit".
  Determines the scale on which the outcome model operates.

## Value

Generated outcome value.

## Details

The non-linear transformation of previous outcome approximates the
natural spline basis (df=3) used in fitting functions. While fitting
uses [`splines::ns()`](https://rdrr.io/r/splines/ns.html), simulation
uses a polynomial approximation for computational efficiency during
one-at-a-time data generation.
