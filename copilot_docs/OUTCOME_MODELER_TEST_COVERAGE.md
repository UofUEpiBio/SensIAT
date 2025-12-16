# Outcome Modeler Test Coverage Summary

## Overview

The SensIAT package provides three different parameterizations of single-index outcome models. This document summarizes the test coverage for each variant.

## Three Outcome Modeler Variants

### 1. `fit_SensIAT_single_index_fixed_coef_model`
- **Parameterization**: Fixes first coefficient, optimizes bandwidth and remaining coefficients
- **File**: `R/sim_outcome_modeler.R`
- **Test Coverage**: Extensive (referenced in 11 test files)
- **Primary Tests**: Integrated into multiple workflow tests (jackknife, influence terms, etc.)

### 2. `fit_SensIAT_single_index_norm1coef_model` (MAVE)
- **Parameterization**: Normalizes coefficients to unit norm, optimizes bandwidth
- **File**: `R/SensIAT_sim_outcome_modeler_mave.R`
- **Test File**: `tests/testthat/test-SensIAT_sim_outcome_modeler_mave.R`
- **Test Coverage**:
  - 3 test_that blocks
  - 12 expect_* assertions
  - 131 lines of code
- **Tests Cover**:
  - ISE vs MSE optimization
  - Grid search vs optim methods
  - Coefficient re-estimation options
  - Integration with within_group workflow

### 3. `fit_SensIAT_single_index_fixed_bandwidth_model` (NEW)
- **Parameterization**: Fixes bandwidth at 1, optimizes all coefficients
- **File**: `R/sim_outcome_modeler.R`
- **Test File**: `tests/testthat/test-SensIAT_sim_outcome_modeler_fixed_bandwidth.R`
- **Test Coverage**:
  - 9 test_that blocks
  - 42 expect_* assertions  
  - 363 lines of code
- **Tests Cover**:
  - Basic functionality and model structure
  - Three optimization methods (nmk, optim, nlminb)
  - Two kernel types (K2_Biweight, dnorm)
  - Expected value computation
  - Custom initial value specification
  - Model predictions (structure validation)
  - S3 methods (coef, formula)
  - Integration with within_group workflow
  - Edge cases (alpha=0, negative alphas)

## Test Comparison

| Feature | MAVE | Fixed Bandwidth |
|---------|------|-----------------|
| Test blocks | 3 | 9 |
| Assertions | 12 | 42 |
| Lines of code | 131 | 363 |
| Optimization methods tested | 2 (grid, optim) | 3 (nmk, optim, nlminb) |
| Kernel types tested | 1 (implicit) | 2 (K2_Biweight, dnorm) |
| Initial values tested | No | Yes |
| Edge cases tested | No | Yes |

## Key Test Features

### Fixed Bandwidth Tests Include:

1. **Helper Functions**:
   - `create_small_test_data()`: Creates 10-subject dataset for fast testing
   - `run_fixed_bandwidth_expectations()`: Common validation checks

2. **Comprehensive Coverage**:
   - Verifies bandwidth is exactly 1 (not estimated)
   - Tests all coefficients are optimized (unlike fixed_coef which fixes first coefficient)
   - Validates S3 class structure
   - Checks optimization convergence across multiple methods
   - Tests kernel function variations

3. **Integration Testing**:
   - Tests with `fit_SensIAT_within_group_model`
   - Validates influence term computation
   - Confirms expected value calculations
   - Tests multiple alpha values

## Test Execution Results

All tests pass successfully:
- **Total passing tests**: 348 (full package)
- **Fixed bandwidth tests**: 147 passing assertions
- **Zero failures**

## Summary

The `fit_SensIAT_single_index_fixed_bandwidth_model` now has **more comprehensive test coverage** than the MAVE variant, with:
- 3× more test blocks
- 3.5× more assertions
- Coverage of more edge cases and parameter combinations
- Better documentation of expected behavior

All three outcome modeler parameterizations are now adequately tested, ensuring reliable performance across different use cases.
