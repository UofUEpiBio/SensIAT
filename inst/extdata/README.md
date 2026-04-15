# Precomputed Data Files

This directory contains precomputed data used by package vignettes to avoid
lengthy computation times during rendering.

## Files

### term2_benchmark_results.rda

Precomputed benchmark results for the `term2-integration-methods` vignette.
Contains:

- **`term2_benchmark_results`**: List with timing and accuracy comparison 
  of term2 integration methods (`fast`, `fixed_grid`, `seeded_adaptive`)
  
- **`term2_grid_benchmark_results`**: List with timing and accuracy for
  fixed-grid method at varying grid densities (50, 100, 200 points)

Both objects include:
- `$timing`: Data frame with mean/sd/min/max times per method
- `$accuracy`: Data frame with max_abs_diff, mean_abs_diff, rmse vs reference
- `$setup_info`: List with benchmark parameters (n_patients, n_alphas, etc.)

## Regenerating Data

To regenerate the benchmark data:

```r
source("inst/extdata/generate_term2_benchmarks.R")
```

Or from the command line:

```bash
Rscript inst/extdata/generate_term2_benchmarks.R
```

## Usage in Vignettes

```r
# Load precomputed results
load(system.file("extdata", "term2_benchmark_results.rda", package = "SensIAT"))

# Access benchmark data
term2_benchmark_results$timing
term2_benchmark_results$accuracy
```
