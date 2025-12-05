# Runs an optimized implementation of the `compute_influence_term_2_quadv_sim_via_matrix` function.

Runs an optimized implementation of the
`compute_influence_term_2_quadv_sim_via_matrix` function.

## Usage

``` r
pcoriaccel_compute_influence_term_2_quadv_sim_via_matrix(
  X,
  Y,
  times,
  individual_X,
  x_slope,
  alpha,
  beta,
  spline_basis,
  bandwidth,
  tol = 0.0001220703,
  kernel = "K2_Biweight"
)
```

## Arguments

- X:

  Matrix of all covariates, transformed as necessary by model

- Y:

  Vector of all outcomes (same length as a column of `X`)

- times:

  Vector of observation times for individual

- individual_X:

  Matrix of covariates for individual rows correspond to times prepared
  for inferences for integration.

- x_slope:

  Vector of numeric(length(beta)) indicating how

- alpha:

  Vector of sensitivity parameters

- beta:

  Vector of coefficients of the outcome model

- spline_basis:

  Spline basis object
  ([`orthogonalsplinebasis::SplineBasis`](https://rdrr.io/pkg/orthogonalsplinebasis/man/SplineBasis.html))

- bandwidth:

  Bandwidth for the kernel density estimate of the outcome model.

- tol:

  Tolerance for integration

- kernel:

  Kernel function to use for the kernel density estimate

## Value

integration result
