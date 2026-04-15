# Compute expected values at grid points for a patient

Compute expected values at grid points for a patient

## Usage

``` r
compute_expected_values_at_grid(
  grid,
  patient_data,
  outcome_model,
  alpha,
  impute_fn = NULL,
  time_var = NULL
)
```

## Arguments

- grid:

  Numeric vector of time points

- patient_data:

  Patient data frame

- outcome_model:

  Fitted outcome model

- alpha:

  Sensitivity parameter

- impute_fn:

  Imputation function

- time_var:

  Time variable name (optional)

## Value

List of expected value lists at each grid point
