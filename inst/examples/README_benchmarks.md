# Term 2 Influence Computation Benchmark

This directory contains benchmarking scripts for optimizing the term 2 influence computation in the generalized marginal mean model fitting.

## Files

- `benchmark_term2_caching.R` - Compares original vs cached/memoized implementation

## Running the Benchmark

```r
source("inst/examples/benchmark_term2_caching.R")
```

## What It Tests

The benchmark compares two implementations of `compute_term2_for_patient`:

1. **Original Method**: Computes expected values from scratch at each integration point
2. **Cached Method**: Uses memoization to avoid recomputing at repeated integration points

## Output

The script provides:

- **Timing comparison**: Total execution time for both methods
- **Validation**: Verifies that both methods produce identical results
- **Per-patient breakdown**: Shows computation time per patient
- **Speedup metrics**: Calculates performance improvement
- **Extrapolation**: Estimates time needed for full dataset

## Customization

You can modify the test parameters:

- `test_ids`: Change which patients to test (default: first 5)
- `alpha`: Sensitivity parameter (default: 0)
- `knots`: Spline knot locations
- `beta_test`: Initial beta vector for testing

## Expected Results

The cached method should show significant speedup (2-5x) since:
- Integration uses adaptive Simpson's rule with repeated evaluations at same points
- Expected value computation involves expensive kernel density estimation
- Caching eliminates redundant calculations

## Next Steps

If caching shows improvement:
1. Apply to full `fit_SensIAT_marginal_mean_model_generalized()` function
2. Consider additional optimizations (pre-computation, parallelization)
3. Profile to identify remaining bottlenecks
