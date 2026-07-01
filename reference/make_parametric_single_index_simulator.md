# Build a parametric outcome simulator from a fitted single-index outcome model using sampled coefficients for the single-index projection.

Build a parametric outcome simulator from a fitted single-index outcome
model using sampled coefficients for the single-index projection.

## Usage

``` r
make_parametric_single_index_simulator(
  outcome_model,
  sampled_coef,
  covariate_mapping = NULL
)
```

## Arguments

- outcome_model:

  A fitted `SensIAT::Single-index-outcome-model` object.

- sampled_coef:

  Named numeric vector of sampled coefficients (from
  [`get_bootstrap_coeffs()`](https://uofuepibio.github.io/SensIAT/reference/get_bootstrap_coeffs.md)).

- covariate_mapping:

  Optional named character vector mapping expected covariate names to
  model variable names. **\[experimental\]**
