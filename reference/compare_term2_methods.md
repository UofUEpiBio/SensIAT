# Compare Vectorized vs Original Term2 Integration Methods

Test function to verify that the new vectorized method produces the same
results as the original method.

## Usage

``` r
compare_term2_methods(
  patient_data,
  outcome_model,
  base,
  alpha_vec,
  marginal_beta,
  V_inv,
  tmin,
  tmax,
  impute_fn,
  inv_link,
  method_original = compute_term2_influence_original,
  tolerance = 1e-10
)
```

## Arguments

- patient_data:

  data.frame with patient's observations

- outcome_model:

  The fitted outcome model (any type supported by
  compute_SensIAT_expected_values)

- base:

  `SplineBasis` for marginal mean model

- alpha_vec:

  Numeric vector of sensitivity parameters

- marginal_beta:

  Coefficients for the marginal mean spline basis

- V_inv:

  Inverse Gram matrix for base

- tmin:

  Lower integration bound

- tmax:

  Upper integration bound

- inv_link:

  Inverse link function

- method_original:

  Function implementing original integration method

- tolerance:

  Numerical tolerance for comparison

## Value

List with comparison results and diagnostics
