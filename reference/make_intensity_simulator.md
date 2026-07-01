# Create an intensity simulator from a fitted intensity model

This factory builds a closure that computes an observation intensity
(hazard) at time t given covariates like `prev_outcome` and `visit_num`.
Supports user-supplied functions and `coxph` fitted models (approximate
baseline hazard via
[`survival::basehaz`](https://rdrr.io/pkg/survival/man/basehaz.html)).

## Usage

``` r
make_intensity_simulator(intensity_model, covariate_mapping = NULL)
```

## Arguments

- intensity_model:

  A fitted model (e.g., `coxph`) or a
  `function(t, prev_outcome, visit_num)`.

- covariate_mapping:

  Optional named character vector mapping expected covariate names to
  model variable names.

## Value

A function of signature `function(t, prev_outcome, visit_num)` returning
non-negative numeric intensity.
