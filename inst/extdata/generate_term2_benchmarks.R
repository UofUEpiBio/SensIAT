# Generate precomputed benchmark data for term2-integration-methods vignette
#
# This script generates the benchmark results stored in term2_benchmark_results.rda
# Run from package root: Rscript inst/extdata/generate_term2_benchmarks.R
#
# Output: inst/extdata/term2_benchmark_results.rda containing:
#   - term2_benchmark_results: Main method comparison benchmark
#   - term2_grid_benchmark_results: Grid density analysis benchmark

library(SensIAT)
library(dplyr)

cat("Generating term2 benchmark data for vignette...\n\n")

# Prepare example data with lag variables
data_with_lags <- SensIAT_example_data |>
  group_by(Subject_ID) |>
  mutate(
    ..prev_outcome.. = lag(Outcome, default = NA_real_, order_by = Time),
    ..prev_time.. = lag(Time, default = 0, order_by = Time),
    ..delta_time.. = Time - lag(Time, default = NA_real_, order_by = Time)
  ) |>
  ungroup()

# Fit outcome model
outcome.model <- fit_SensIAT_single_index_fixed_coef_model(
  Outcome ~ splines::ns(..prev_outcome.., df = 3) + ..delta_time.. - 1,
  id = Subject_ID,
  data = data_with_lags |> filter(Time > 0)
)

# Imputation function
impute_fn <- function(t, df) {
  extrapolate_from_last_observation(t, df, "Time", slopes = c("..delta_time.." = 1))
}

knots <- c(60, 260, 460)

# -----------------------------------------------------------------------------
# Main benchmark: Compare methods
# -----------------------------------------------------------------------------
cat("Running main benchmark (method comparison)...\n")

term2_benchmark_results <- benchmark_term2_methods(
  data = data_with_lags,
  id = Subject_ID,
  time = data_with_lags$Time,
  outcome.model = outcome.model,
  knots = knots,
  alpha = 0,
  impute_data = impute_fn,
  link = "identity",
  methods = c("fast", "fixed_grid", "seeded_adaptive"),
  grid_sizes = c(50, 100),
  n_patients = 5,
  n_iterations = 2,
  reference_method = "fast"
)

cat("\nMain benchmark timing:\n")
print(term2_benchmark_results$timing)

# -----------------------------------------------------------------------------
# Grid analysis benchmark: Varying grid sizes
# -----------------------------------------------------------------------------
cat("\nRunning grid density analysis benchmark...\n")

term2_grid_benchmark_results <- benchmark_term2_methods(
  data = data_with_lags,
  id = Subject_ID,
  time = data_with_lags$Time,
  outcome.model = outcome.model,
  knots = knots,
  alpha = 0,
  impute_data = impute_fn,
  link = "identity",
  methods = c("fast", "fixed_grid"),
  grid_sizes = c(50, 100, 200),
  n_patients = 5,
  n_iterations = 2,
  reference_method = "fast"
)

cat("\nGrid benchmark timing:\n")
print(term2_grid_benchmark_results$timing)

# -----------------------------------------------------------------------------
# Save results
# -----------------------------------------------------------------------------
output_file <- file.path("inst", "extdata", "term2_benchmark_results.rda")
save(term2_benchmark_results, term2_grid_benchmark_results, file = output_file)
cat("\nSaved to:", output_file, "\n")

# Record metadata
cat("\nBenchmark metadata:\n")
cat("  Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n
")
cat("  R version:", R.version.string, "\n")
cat("  Platform:", R.version$platform, "\n")
