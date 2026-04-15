# Pre-compute Gauss-Legendre grid for term2 integration

Similar to pre_compute_term2_grid but optimized for Gauss-Legendre
quadrature. Pre-computes the expected values at Gauss-Legendre nodes for
reuse across beta iterations.

## Usage

``` r
pre_compute_gauss_legendre_grid(
  patient_data,
  outcome_model,
  base,
  alpha,
  tmin,
  tmax,
  impute_fn = NULL,
  n_nodes = 50,
  time_var = NULL
)
```

## Arguments

- patient_data:

  Patient data frame

- outcome_model:

  Fitted outcome model

- base:

  `SplineBasis` object

- alpha:

  Sensitivity parameter

- tmin:

  Lower integration bound

- tmax:

  Upper integration bound

- impute_fn:

  Imputation function

- n_nodes:

  Number of Gauss-Legendre nodes

- time_var:

  Time variable name (optional)

## Value

List with:

- nodes: Gauss-Legendre nodes transformed to the interval (`tmin`,
  `tmax`)

- weights: Gauss-Legendre weights (including Jacobian)

- B_nodes: Basis evaluations at nodes

- E_nodes: Expected values at nodes
