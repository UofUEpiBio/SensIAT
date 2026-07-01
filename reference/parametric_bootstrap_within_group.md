# Parametric bootstrap for a within-group SensIAT model **\[experimental\]**

Parametric bootstrap for a within-group SensIAT model
**\[experimental\]**

## Usage

``` r
parametric_bootstrap_within_group(
  within_group_model,
  nboot = 100,
  simulate_args = list(),
  seed = NULL,
  progress = interactive(),
  sample_coefficients = FALSE,
  refit = TRUE,
  return = c("coefficients", "models", "data"),
  prune_models = FALSE,
  gc_every = 10L,
  verbosity = c("none", "basic", "detailed"),
  verbose = NULL
)
```

## Arguments

- within_group_model:

  A fitted `SensIAT_within_group_model` object.

- nboot:

  Number of bootstrap replicates.

- simulate_args:

  List of arguments to pass to
  [`simulate_SensIAT_data()`](https://uofuepibio.github.io/SensIAT/reference/simulate_SensIAT_data.md).
  If not specified, `End`, `n_subjects`, `initial_outcome_mean`, and
  `initial_outcome_sd` are inferred from the fitted model.

- seed:

  Optional seed for reproducibility.

- progress:

  Logical; show progress bar when available.

- sample_coefficients:

  Logical; if `TRUE`, sample coefficients from an asymptotic
  multivariate normal distribution when
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) is available. If `FALSE`
  (default), use original fitted coefficients.

- refit:

  Logical; if `TRUE` (default), fit a `SensIAT_within_group_model` on
  each simulated replicate using the original model's settings.

- return:

  One of `"coefficients"` (default), `"models"`, or `"data"`.
  `"coefficients"` is memory-efficient and stores only replicated
  marginal mean coefficients.

- prune_models:

  Logical; when `return = "models"`, prune each replicated model before
  returning.

- gc_every:

  Integer. Run `gc(FALSE)` every `gc_every` replications. Use `NULL` to
  disable explicit garbage collection.

- verbosity:

  Logging verbosity for bootstrap internals: one of `"none"`, `"basic"`,
  or `"detailed"`.

- verbose:

  Deprecated shortcut; if `TRUE`, equivalent to
  `verbosity = "detailed"`.

## Value

If `return = "coefficients"`, a `SensIAT_withingroup_bootstrap_results`
object. Otherwise returns a list of replicated fitted models
(`"models"`) or simulated datasets (`"data"`).

## Examples

``` r
if (FALSE) { # \dontrun{
data("SensIAT_example_data", package = "SensIAT")

# Fit a single-index outcome model on a small subset of the example data.
small_data <- dplyr::filter(
  SensIAT_example_data,
  Subject_ID %in% head(unique(SensIAT_example_data$Subject_ID), 8)
)

model <- fit_SensIAT_within_group_model(
  group.data = small_data,
  outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
  alpha = 0,
  id = Subject_ID,
  outcome = Outcome,
  time = Time,
  End = 830,
  knots = c(60, 260, 460)
)

# This example may take a long time because it fits a single-index outcome model
# and generates bootstrap replicates.
res <- parametric_bootstrap_within_group(
  nboot = 2,
  within_group_model = model,
  simulate_args = list(
    n_subjects = 3,
    End = 5,
    max_visits = 5
  ),
  seed = 123,
  refit = TRUE
)

print(res[[1]])
} # }
```
