# Parametric bootstrap orchestration **\[experimental\]**

Parametric bootstrap orchestration **\[experimental\]**

## Usage

``` r
parametric_bootstrap(
  nboot = 100,
  intensity_model = NULL,
  outcome_model = NULL,
  simulate_args = list(),
  seed = NULL,
  progress = interactive(),
  sample_coefficients = FALSE,
  verbosity = c("none", "basic", "detailed"),
  verbose = NULL
)
```

## Arguments

- nboot:

  Number of bootstrap replicates

- intensity_model:

  Fitted intensity model
  ([survival::coxph](https://rdrr.io/pkg/survival/man/coxph.html)) or
  function

- outcome_model:

  Fitted outcome model (single-index) or NULL

- simulate_args:

  List of arguments to pass to
  [`simulate_SensIAT_data()`](https://uofuepibio.github.io/SensIAT/reference/simulate_SensIAT_data.md)
  (e.g., n_subjects, End, intensity_bound)

- seed:

  Optional seed for reproducibility

- progress:

  Logical; show progress bar when available.

- sample_coefficients:

  Logical; if `TRUE`, sample coefficients from an asymptotic
  multivariate normal distribution when
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) is available. If `FALSE`
  (default), use original fitted coefficients.

- verbosity:

  Logging verbosity for bootstrap internals: one of `"none"`, `"basic"`,
  or `"detailed"`.

- verbose:

  Deprecated shortcut; if `TRUE`, equivalent to
  `verbosity = "detailed"`.

## Value

A list of simulated datasets (length `nboot`)
