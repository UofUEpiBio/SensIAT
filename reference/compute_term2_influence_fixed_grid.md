# Compute term2 influence using fixed-grid integration

This method evaluates the integrand at a pre-determined grid of points
and uses composite Simpson's rule for integration. Expected values are
pre-computed at grid points and reused across beta iterations.

## Usage

``` r
compute_term2_influence_fixed_grid(
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
  n_grid = 100,
  rule = c("simpson", "trapezoid"),
  include_obs_times = TRUE,
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

  Function to impute data at time t (not used if expected_grid provided)

- inv_link:

  Inverse link function

- W:

  Weight function W(t, beta)

- expected_grid:

  Optional pre-computed expected values at grid points: list with
  elements: - grid: numeric vector of time points - B_grid: list of
  basis evaluations at grid points - E_grid: list of expected values
  list(E_exp_alphaY, E_Yexp_alphaY)

- n_grid:

  Number of grid points to use if expected_grid not provided

- rule:

  Integration rule: "simpson" (default) or "trapezoid"

- include_obs_times:

  Whether to include observation times in grid (default TRUE)

- time_var:

  Time variable name (optional)

- ...:

  Additional arguments (not used)

## Value

Numeric vector of length `ncol(base)` with term2 influence values
