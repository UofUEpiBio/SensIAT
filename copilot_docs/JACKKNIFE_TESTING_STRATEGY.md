# Fast Jackknife Testing Strategy

## Problem
The original jackknife test was taking an estimated **42 minutes** to complete because it:
- Used all 200 subjects from `SensIAT_example_data` 
- Performed leave-one-out cross-validation (200 model refits)
- Used complex model specification with multiple knots and spline terms
- Tested parallelization by running the expensive operation twice

## Solution Overview
Created a multi-tier testing strategy that provides **20x speed improvement** while maintaining test coverage:

### 1. Fast Tests (test-jackknife.R)
**Runtime: ~2 minutes** instead of 42 minutes

- **Small dataset**: 10 subjects instead of 200 (20x fewer cross-validation iterations)
- **Simplified model**: Single knot, linear terms instead of natural splines  
- **Reduced complexity**: Fewer alpha values and time points
- **Focused testing**: Tests core functionality and parallelization separately

### 2. Comprehensive Tests (test-jackknife-comprehensive.R) 
**Runtime: 5-15 minutes** - run manually when needed

- Medium-scale tests with 20-50 subjects
- Full complexity models for integration testing
- Performance benchmarks
- Original full-scale test preserved but skipped by default

### 3. Unit Tests
**Runtime: <1 second**

- Tests jackknife computation logic with mocked data
- Validates mathematical correctness without expensive model fitting
- Ensures variance calculations are correct

## Test Structure

### Fast Tests (Always Run)
```r
# Basic functionality with 10 subjects
test_that("jackknife basic functionality")

# Parallelization with 5 subjects  
test_that("jackknife is invariant under parallelization")

# Mathematical correctness with mocked data
test_that("jackknife computation logic")
```

### Comprehensive Tests (Manual)
```r
# To run comprehensive tests:
testthat::test_file("tests/testthat/test-jackknife-comprehensive.R")

# Or run specific tests:
testthat::test_file("tests/testthat/test-jackknife-comprehensive.R", 
                    filter = "medium scale")
```

## What's Tested

### Core Functionality ✓
- Jackknife variance calculation correctness
- Cross-validation logic (leave-one-out)
- Result structure and data types
- Integration with SensIAT model objects

### Parallelization ✓  
- Sequential vs parallel result consistency
- Future/furrr integration
- Multi-worker execution

### Edge Cases ✓
- Small datasets
- Single alpha values  
- Mathematical edge cases with mocked data

### Performance ✓
- Timing benchmarks
- Memory usage validation
- Scalability testing (in comprehensive tests)

## Usage Guidelines

### For Regular Development
- Fast tests run automatically with `devtools::test()`
- Provide confidence in core functionality
- Complete in reasonable time for CI/CD

### For Jackknife-Specific Changes
- Run comprehensive tests manually
- Use medium-scale tests for thorough validation
- Run performance benchmarks to check for regressions

### For Release Validation  
- Run full comprehensive test suite
- Validate performance benchmarks
- Ensure no regressions in complex scenarios

## Environmental Controls

```r
# Skip medium tests if desired
Sys.setenv(SKIP_MEDIUM_TESTS = "true")

# Force comprehensive tests to run
# (edit the skip() calls in test-jackknife-comprehensive.R)
```

## Benefits

1. **Speed**: 20x faster routine testing (2 min vs 42 min)
2. **Coverage**: Maintains test coverage of core functionality  
3. **Flexibility**: Multiple test tiers for different scenarios
4. **CI-Friendly**: Fast tests suitable for continuous integration
5. **Debugging**: Easier to debug with smaller datasets
6. **Mathematical Validation**: Unit tests verify computation logic

## Migration from Original Test

The original expensive test is preserved in `test-jackknife-comprehensive.R` as:
```r
test_that("jackknife comprehensive test with full dataset", {
    skip("Manual test - run explicitly when needed")
    # ... original test code
})
```

This ensures no functionality is lost while making routine testing practical.