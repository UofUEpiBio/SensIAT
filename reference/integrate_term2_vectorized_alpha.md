# Vectorized Term2 Integration Across Multiple Alpha Values

Integrates the term2 influence function across multiple alpha values
simultaneously, sharing expensive computations while maintaining
separate convergence criteria.

## Usage

``` r
integrate_term2_vectorized_alpha(
  compute_expected_values_fn,
  impute_fn,
  weight_fn,
  marginal_mean_fn,
  alpha_vec,
  tmin,
  tmax,
  patient_data,
  tol = 1.490116e-08
)
```

## Arguments

- compute_expected_values_fn:

  R function that computes expected values for all alphas

- impute_fn:

  R function for data imputation: impute_fn(t, patient_data)

- weight_fn:

  R function for weight computation: weight_fn(t)

- marginal_mean_fn:

  R function for marginal mean: marginal_mean_fn(t)

- alpha_vec:

  Numeric vector of sensitivity parameters

- tmin:

  Lower integration bound

- tmax:

  Upper integration bound

- patient_data:

  Patient data object (passed to R functions)

- tol:

  Convergence tolerance

## Value

List of integration results, one per alpha value
