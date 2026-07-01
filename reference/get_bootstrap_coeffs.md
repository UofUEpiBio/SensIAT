# Get coefficients for bootstrap simulation.

By default this returns original fitted coefficients. If sampling is
requested, coefficients are sampled from an asymptotic multivariate
normal distribution only when
[`vcov()`](https://rdrr.io/r/stats/vcov.html) is available.

## Usage

``` r
get_bootstrap_coeffs(model, sample_coefficients = FALSE, model_label = "model")
```

## Arguments

- model:

  A fitted model with [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) methods.

- sample_coefficients:

  Logical; if `TRUE`, sample coefficients from an asymptotic
  multivariate normal distribution when
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) is available. If `FALSE`
  (default), use original fitted coefficients.

- model_label:

  Character; label for the model used in warning messages.
