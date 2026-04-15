# Build fast term2 integrand using closure optimization

Build fast term2 integrand using closure optimization

Build a fast term2 integrand for a single patient (closure)

## Usage

``` r
make_term2_integrand_fast(
  outcome.model,
  base,
  alpha,
  patient_times,
  patient_outcomes,
  marginal_beta,
  V_inv = NULL,
  W = NULL,
  expected_get = NULL
)

make_term2_integrand_fast(
  outcome.model,
  base,
  alpha,
  patient_times,
  patient_outcomes,
  marginal_beta,
  V_inv = NULL,
  W = NULL,
  expected_get = NULL
)
```

## Arguments

- outcome.model:

  The fitted Single-index outcome model

- base:

  `SplineBasis` object for marginal mean model

- alpha:

  Sensitivity parameter

- patient_times:

  Numeric vector of observation times for the patient (sorted asc)

- patient_outcomes:

  Numeric vector of outcomes aligned with patient_times

- marginal_beta:

  Coefficients (beta) for the marginal mean spline basis

- V_inv:

  Precomputed inverse Gram matrix for base (optional; computed if NULL)

- W:

  Weight function W(t, beta)

- expected_get:

  Optional caching function for expected values: expected_get(t) -\>
  list(E_exp_alphaY, E_Yexp_alphaY)

## Value

A function f(t) computing the term2 integrand at scalar t
