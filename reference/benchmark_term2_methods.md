# Benchmark Term2 Integration Methods

Efficiently compares different term2 integration methods by isolating
just the term2 computation from the full model fitting. This is useful
for performance analysis and method selection without the overhead of
repeated model fitting.

## Usage

``` r
benchmark_term2_methods(
  data,
  id,
  time,
  outcome.model,
  knots,
  spline_degree = 3L,
  alpha = 0,
  impute_data,
  link = c("identity", "log"),
  methods = c("fast", "original", "fixed_grid", "seeded_adaptive"),
  grid_sizes = c(50, 100, 200),
  n_patients = NULL,
  n_iterations = 3,
  reference_method = "fast",
  seed = 42
)
```

## Arguments

- data:

  Data frame with longitudinal observations

- id:

  Unquoted column name for subject identifier

- time:

  Numeric vector of observation times

- outcome.model:

  A fitted outcome model (e.g., from
  fit_SensIAT_single_index_fixed_coef_model)

- knots:

  Knots for marginal mean model spline basis

- spline_degree:

  Degree of B-spline basis (default: 3)

- alpha:

  Numeric vector of sensitivity parameters to test

- impute_data:

  Function to impute data at arbitrary times:
  `function(t, patient_data) -> data.frame`

- link:

  Link function: "identity" or "log"

- methods:

  Character vector of methods to benchmark. Options: "fast", "original",
  "fixed_grid", "seeded_adaptive", "gauss_legendre"

- grid_sizes:

  For grid-based methods, vector of grid sizes to test (default: c(50,
  100, 200))

- n_patients:

  Number of patients to use for benchmark (NULL = all patients)

- n_iterations:

  Number of timing iterations per method

- reference_method:

  Method to use as accuracy reference (default: "fast")

- seed:

  Random seed for reproducibility

## Value

A list with components:

- timing:

  Data frame with timing results per method

- accuracy:

  Data frame with accuracy metrics vs reference

- reference_results:

  List of reference results per patient

- setup_info:

  List with setup parameters and dimensions

## Examples

``` r
if (FALSE) { # \dontrun{
# Setup data with lag variables
data_with_lags <- SensIAT_example_data |>
  dplyr::group_by(Subject_ID) |>
  dplyr::mutate(
    ..prev_outcome.. = dplyr::lag(Outcome, order_by = Time),
    ..prev_time.. = dplyr::lag(Time, default = 0, order_by = Time),
    ..delta_time.. = Time - dplyr::lag(Time, order_by = Time)
  ) |>
  dplyr::ungroup()

# Fit outcome model
outcome.model <- fit_SensIAT_single_index_fixed_coef_model(
  Outcome ~ splines::ns(..prev_outcome.., df = 3) + ..delta_time.. - 1,
  id = Subject_ID,
  data = data_with_lags |> dplyr::filter(Time > 0)
)

# Imputation function
impute_fn <- function(t, df) {
  extrapolate_from_last_observation(t, df, "Time",
    slopes = c("..delta_time.." = 1))
}

# Run benchmark
results <- benchmark_term2_methods(
  data = data_with_lags,
  id = Subject_ID,
  time = data_with_lags$Time,
  outcome.model = outcome.model,
  knots = c(60, 260, 460),
  alpha = c(-0.1, 0, 0.1),
  impute_data = impute_fn,
  methods = c("fast", "fixed_grid", "seeded_adaptive"),
  grid_sizes = c(50, 100),
  n_patients = 10
)

# View timing results
print(results$timing)
} # }
```
