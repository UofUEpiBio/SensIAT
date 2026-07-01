# Create a simulator function from a fitted Single-index outcome model

This factory builds a closure that can be used to simulate outcome
values from a fitted `SensIAT::Single-index-outcome-model`. The returned
function samples from the estimated conditional distribution of the
outcome given covariates using the Nadaraya-Watson estimator implemented
in `pcoriaccel_NW`.

## Usage

``` r
make_single_index_simulator(outcome_model, covariate_mapping = NULL)
```

## Arguments

- outcome_model:

  A fitted object of class `SensIAT::Single-index-outcome-model`.

- covariate_mapping:

  Optional named character vector mapping expected covariate names (e.g.
  `prev_outcome`, `time`, `delta_time`) to the names used in the
  original model formula. If `NULL`, the factory will attempt to use
  `prev_outcome`, `time`, and `delta_time` directly.

## Value

A function with signature
`function(prev_outcome, time, delta_time, newdata = NULL)` which returns
a sampled outcome value consistent with the fitted model.
