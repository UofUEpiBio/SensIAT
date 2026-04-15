# BENCHMARK RUNNING INSTRUCTIONS
================================================================================

Since Rscript is not in your PATH, run the benchmark from an R console.

## Steps to Run:

### 1. Open R (or RStudio)

### 2. Load the package
```r
setwd("c:/Users/u0092104/OneDrive - University of Utah/Projects/Dan S - PCORI/pcoriRPackage")
devtools::load_all()  # Or use Ctrl+Shift+L in RStudio
```

### 3. Quick test (verify fast integrand works)
```r
source("inst/examples/test_fast_integrand.R")
```

This should output:
- ✓ Setup complete
- ✓ Fast integrand created  
- ✓ Integrand evaluated
- ✓ Integration successful
- SUCCESS message

### 4. Run full benchmark (if test passes)
```r
source("inst/examples/benchmark_term2_caching.R")
```

This will compare 3 methods:
- **Original**: Uses impute_data() + compute_SensIAT_expected_values()
- **Cached**: Memoizes results (but still calls the expensive functions)
- **FAST**: Low-allocation closure, avoids data.frame creation

## Expected Output:

```
============================================================
BENCHMARK: Term 2 Influence Computation
============================================================

Setup:
  - Number of patients: 5
  - Alpha: 0
  - Integration range: [ 60 , 460 ]
  - Spline basis dimension: 7

Warming up...
Warm-up complete.

Running ORIGINAL method...
Original method completed in X.XX seconds

Running CACHED method...
Cached method completed in X.XX seconds

Running FAST method...
FAST method completed in X.XX seconds

============================================================
VALIDATION: Comparing Results
============================================================

Results match: TRUE

Per-patient comparison:
# A tibble: 5 × 5
  patient_id   original    cached difference rel_diff
       <int>      <dbl>     <dbl>      <dbl>    <dbl>
...

============================================================
PERFORMANCE SUMMARY
============================================================

Original method: X.XX seconds
Cached method:   X.XX seconds
FAST method:     X.XX seconds
Speedup (Cached vs Original): X.XXx
Speedup (FAST vs Original):   X.XXx
Time saved (FAST):            X.XX seconds (XX.X%)

Per-patient average:
Original: X.XXX seconds/patient
Cached:   X.XXX seconds/patient
FAST:     X.XXX seconds/patient

Estimated time for all NN patients:
Original: X.X minutes
Cached:   X.X minutes
FAST:     X.X minutes

============================================================
```

## What to Report:

Please paste the output showing:
1. Whether the test_fast_integrand.R succeeded
2. The timing results from the full benchmark (especially "PERFORMANCE SUMMARY" section)
3. Whether validation passed (Results match: TRUE)
4. The speedup values for FAST vs Original

## Troubleshooting:

If you get errors, paste:
- The error message
- Which step failed (setup, original, cached, or fast)
- The traceback if available

## Files Created:

- R/term2_fast.R - Fast integrand builder
- inst/examples/benchmark_term2_caching.R - Full 3-way benchmark
- inst/examples/test_fast_integrand.R - Quick validation test
- inst/examples/check_setup.R - Dependency checker
- inst/examples/run_benchmark_simple.R - Alternative launcher
- inst/examples/INSTRUCTIONS.md - This file

================================================================================
