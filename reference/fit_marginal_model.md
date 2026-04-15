# Fit SensIAT Marginal Mean Model (Unified)

This function fits a marginal mean model for sensitivity analysis,
supporting both single-index (linear) and generalized linear outcome
models.

## Usage

``` r
fit_marginal_model(
  data,
  outcome.model,
  intensity.model,
  id,
  time,
  alpha = 0,
  knots,
  ...
)
```

## Arguments

- data:

  Data frame containing all required variables.

- outcome.model:

  Outcome model object or formula (GLM, single-index, or LM).

- intensity.model:

  Intensity model object (e.g., `coxph`).

- id:

  Subject identifier variable.

- time:

  Time variable.

- alpha:

  Sensitivity parameter(s).

- knots:

  Spline knot locations.

- ...:

  Additional arguments passed to underlying methods.

## Value

A fitted SensIAT marginal mean model object.
