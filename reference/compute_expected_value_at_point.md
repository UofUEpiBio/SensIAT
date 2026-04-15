# Compute expected value at a single point

Compute expected value at a single point

## Usage

``` r
compute_expected_value_at_point(
  t,
  patient_data,
  outcome_model,
  alpha,
  impute_fn = NULL,
  time_var = NULL
)
```

## Arguments

- t:

  Time point

- patient_data:

  Patient data frame

- outcome_model:

  Fitted outcome model

- alpha:

  Sensitivity parameter

- impute_fn:

  Imputation function

- time_var:

  Time variable name

## Value

List with E_exp_alphaY and E_Yexp_alphaY
