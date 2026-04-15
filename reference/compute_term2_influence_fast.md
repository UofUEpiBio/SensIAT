# Compute term2 influence for a patient (fast method with segmentation)

Segments integration at observation times to handle discontinuities in
the interpolation function, then uses adaptive Simpson's quadrature on
each segment.

## Usage

``` r
compute_term2_influence_fast(
  patient_data,
  outcome_model,
  base,
  alpha,
  marginal_beta,
  V_inv,
  tmin,
  tmax,
  impute_fn,
  inv_link,
  W,
  expected_get = NULL,
  time_var = NULL,
  ...
)
```

## Arguments

- patient_data:

  data.frame with patient's observations (Time, Outcome, and lag
  variables)

- outcome_model:

  The fitted Single-index outcome model

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

  Function to impute data at time t: impute_fn(t, patient_data) -\>
  data.frame

- inv_link:

  Inverse link function (e.g., exp for log link)

- W:

  Weight function W(t, beta)

- expected_get:

  Optional caching function: `expected_get(t)` -\>
  `list(E_exp_alphaY, E_Yexp_alphaY)`

- time_var:

  Name of the time variable in patient_data (if NULL, auto-detected)

- ...:

  Additional arguments (not used)

## Value

Numeric vector of length `ncol(base)` with term2 influence values
