# Vectorized Term2 Integration for Multiple Alpha Values

This is an experimental implementation that integrates term2 influence
across multiple alpha values simultaneously, sharing expensive
computations for improved performance.

## Usage

``` r
compute_term2_influence_vectorized(
  patient_data,
  outcome_model,
  base,
  alpha_vec,
  marginal_beta,
  V_inv,
  tmin,
  tmax,
  variables,
  centering,
  inv_link,
  tol = 1.490116e-08
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

- variables:

  List of variable names (time, outcome, prev_outcome, etc.)

- centering:

  Centering statistics for standardization

- inv_link:

  Inverse link function

- tol:

  Integration tolerance

## Value

List with one element per alpha value, each containing:

- Q:

  Integration result (vector of length ncol(base))

- fcnt:

  Number of function evaluations

- alpha:

  Alpha value

- converged:

  Convergence status
