# Test Suite README

## Debug Output Control

By default, tests run with minimal output. To enable detailed diagnostic output for debugging, set the environment variable:

```bash
export SENSIAT_TEST_DEBUG=1
```

Or in R:

```r
Sys.setenv(SENSIAT_TEST_DEBUG = "1")
devtools::test()  # or test_active_file(), etc.
```

This will show detailed output in tests that support it, including:
- Integration bounds and observation times
- Intermediate calculation values
- Comparison results  
- Diagnostic information

## Test Organization

### Active Tests
- `test-vectorized-integration-simple.R` - Main integration validation (minimal output)
- `test-jackknife*.R` - Jackknife variance estimation tests
- `test-SensIAT_sim_outcome_modeler_mave.R` - MAVE outcome model tests (optimized for speed)

### Debug/Diagnostic Tests (require SENSIAT_TEST_DEBUG=1)
- `test-integration-detailed.R` - Detailed single-piece integration with verbose output
- `test-vectorized-integration-diagnostic.R` - Historical integrand diagnostics (skipped)

### Skipped/Manual Tests  
- `test-vectorized-integration.R` - Skipped, outdated API tests
- `test-jackknife-comprehensive.R` - Manual tests for full dataset validation

To run manual tests:
```r
testthat::test_file("tests/testthat/test-jackknife-comprehensive.R")
```
