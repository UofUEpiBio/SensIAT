# Compute term2 influence using Gauss-Legendre quadrature

This method uses Gauss-Legendre quadrature for numerical integration.
The integrand is evaluated at Gauss-Legendre nodes, which are optimal
for polynomial approximation of smooth functions.

## Usage

``` r
compute_term2_influence_gauss_legendre(
  patient_data,
  outcome_model,
  base,
  alpha,
  marginal_beta,
  V_inv,
  tmin,
  tmax,
  impute_fn = NULL,
  inv_link,
  W,
  expected_grid = NULL,
  n_nodes = 50,
  time_var = NULL,
  ...
)
```

## Arguments

- patient_data:

  data.frame with patient's observations

- outcome_model:

  The fitted outcome model

- base:

  `SplineBasis` object for marginal mean model

- alpha:

  Sensitivity parameter

- marginal_beta:

  Coefficients (beta) for the marginal mean spline basis

- V_inv:

  Inverse Gram matrix for base

- tmin:

  Lower integration bound

- tmax:

  Upper integration bound

- impute_fn:

  Function to impute data at time t

- inv_link:

  Inverse link function

- W:

  Weight function W(t, beta)

- expected_grid:

  Optional pre-computed expected values (not typically used for
  Gauss-Legendre since nodes are specialized)

- n_nodes:

  Number of Gauss-Legendre nodes (default 50)

- time_var:

  Time variable name (optional)

- ...:

  Additional arguments (not used)

## Value

Numeric vector of length `ncol(base)` with term2 influence values

## Details

Gauss-Legendre quadrature approximates the integral as a weighted sum:
\$\$\int_a^b f(x) dx \approx \frac{b-a}{2} \sum\_{i=1}^n w_i f(x_i)\$\$
where \\x_i\\ are the Gauss-Legendre nodes transformed to the interval
\\(a, b)\\ and \\w_i\\ are the corresponding weights.

The method requires the statmod package for computing nodes and weights.
