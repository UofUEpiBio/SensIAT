# Compute expected values at Gauss-Legendre nodes

Compute expected values at Gauss-Legendre nodes

## Usage

``` r
compute_expected_values_at_nodes(
  nodes,
  patient_data,
  outcome_model,
  alpha,
  impute_fn = NULL,
  time_var = NULL
)
```

## Arguments

- nodes:

  Numeric vector of Gauss-Legendre nodes

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

List of expected value lists at each node
