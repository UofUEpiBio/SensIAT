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
- `test-vectorized-integration-simple.R` - Main integration validation (clean output by default)
- `test-integration-detailed.R` - Detailed single-piece integration test (verbose when debug enabled)
- `test-jackknife*.R` - Jackknife variance estimation tests
- `test-SensIAT_sim_outcome_modeler_mave.R` - MAVE outcome model tests (optimized for speed)

### Diagnostic/Skipped Tests  
- `test-vectorized-integration-diagnostic.R` - Skipped, used for historical bug diagnosis
- `test-vectorized-integration.R` - Skipped, outdated API tests
- `test-jackknife-comprehensive.R` - Manual tests for full dataset validation

## Performance

The test suite has been optimized for speed:
- MAVE tests use 10 subjects instead of 200 (~12 seconds)
- Jackknife tests use minimal subjects and simplified models
- Comprehensive tests are skipped by default

To run comprehensive tests:
```r
testthat::test_file("tests/testthat/test-jackknife-comprehensive.R")
```
