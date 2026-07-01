# Build a parametric intensity simulator using sampled coefficients

Supports `coxph` objects: sampled coefficients are used to compute
linear predictor while baseline hazard is taken from
[`survival::basehaz()`](https://rdrr.io/pkg/survival/man/basehaz.html).

## Usage

``` r
make_parametric_intensity_simulator(
  intensity_model,
  sampled_coef,
  covariate_mapping = NULL
)
```

## Arguments

- intensity_model:

  A fitted model (e.g., `coxph`) or a
  `function(t, prev_outcome, visit_num)`.

- sampled_coef:

  Named numeric vector of sampled coefficients (from
  [`get_bootstrap_coeffs()`](https://uofuepibio.github.io/SensIAT/reference/get_bootstrap_coeffs.md)).

- covariate_mapping:

  Optional named character vector mapping expected covariate names to
  model variable names.

## Value

A function of signature `function(t, prev_outcome, visit_num)` returning
non-negative numeric intensity.

**\[experimental\]**
