# Function Removal Analysis for SensIAT Package

**Date**: December 16, 2025  
**Branch**: generalized-fit-optimization  
**Purpose**: Identify functions that could be removed to reduce testing coverage requirements

**Status Update**: ✅ **4 utility functions internalized** - `add_class`, `add_terminal_observations`, `extrapolate_from_last_observation`, `extrapolate_from_last_observation_multiple`

## Executive Summary

This report analyzes the SensIAT package's exported functions to identify candidates for removal or conversion to internal-only status. Functions are categorized by their usage in tests, dependencies, and role in the package API.

**Completed Actions:**
- ✅ Internalized 4 utility functions (removed from NAMESPACE)
- ✅ Added comprehensive test suite for `add_class()` (16 test cases)
- ✅ Verified all existing tests still pass (jackknife: 23, extrapolate: 47, add_class: 16)
- ✅ Updated documentation to mark functions as `@keywords internal`

**Remaining Candidates:**
- 3 high-priority removal candidates
- 2 medium-priority internalization candidates

## Methodology

1. **Identified all exported functions** from NAMESPACE
2. **Searched for function calls** in test files and R source files  
3. **Analyzed dependency chains** to understand which functions call others
4. **Categorized** functions by testing status and usage patterns

## Exported Functions Analysis

### Top-Level User-Facing Functions

These are the main entry points users would call:

#### 1. `fit_SensIAT_fulldata_model()` - **REMOVAL CANDIDATE**
- **Status**: NOT TESTED
- **Defined in**: `R/fulldata_model.R`
- **Purpose**: Fit sensitivity analysis for both treatment and control groups
- **Called by**: Only used in examples in `autoplot.R` documentation
- **Dependencies**: Calls `fit_SensIAT_within_group_model()`
- **Test usage**: 0 test files
- **Recommendation**: **REMOVE or mark as internal**
  - Thin wrapper around `fit_SensIAT_within_group_model()`
  - No test coverage
  - Users can easily call `fit_SensIAT_within_group_model()` twice
  - Only used in documentation examples

#### 2. `fit_SensIAT_within_group_model()` - **KEEP**
- **Status**: TESTED
- **Defined in**: `R/within_group_model.R`
- **Purpose**: Main fitting function for within-group sensitivity analysis
- **Test files**: 8 different test files
- **Recommendation**: **KEEP** - Core API function

#### 3. `fit_SensIAT_marginal_mean_model()` - **POSSIBLE REMOVAL**
- **Status**: NOT TESTED  
- **Defined in**: `R/fit_marginal_model.R`
- **Called by**: 
  - `fit_SensIAT_marginal_mean_model_generalized.R`
  - `within_group_model.R`
- **Test usage**: 0 test files
- **Recommendation**: **Consider internalizing**
  - No direct test coverage
  - Called internally by `fit_SensIAT_within_group_model()`
  - May have been superseded by generalized version

### Outcome Modelers

#### 4. `fit_SensIAT_single_index_fixed_coef_model()` - **KEEP**
- **Status**: TESTED (7 test files)
- **Purpose**: Outcome modeler with fixed first coefficient
- **Recommendation**: **KEEP** - Actively used and tested

#### 5. `fit_SensIAT_single_index_fixed_bandwidth_model()` - **REMOVAL CANDIDATE**
- **Status**: NOT TESTED
- **Defined in**: `R/sim_outcome_modeler.R`
- **Test usage**: 0 test files
- **Used in**: Only in `inst/examples/basic.R`
- **Recommendation**: **REMOVE or mark as internal**
  - No test coverage
  - Only used in one example file
  - Package has other outcome modelers that are tested

#### 6. `fit_SensIAT_single_index_norm1coef_model()` - **KEEP**
- **Status**: TESTED (test-SensIAT_sim_outcome_modeler_mave.R)
- **Purpose**: MAVE-based outcome modeler
- **Recommendation**: **KEEP** - Recently fixed and actively tested

### Support Functions

#### 7. `jackknife()` - **KEEP**
- **Status**: TESTED (3 test files, 23 passing assertions)
- **Purpose**: Jackknife variance estimation
- **Recommendation**: **KEEP** - Core statistical method

#### 8. `prepare_SensIAT_data()` - **KEEP**
- **Status**: TESTED (test-NW.R, test-splinebasis.R)
- **Purpose**: Data preparation with lag variables
- **Called by**: `within_group_model.R`
- **Recommendation**: **KEEP** - Essential preprocessing

#### 9. `compute_SensIAT_expected_values()` - **KEEP**
- **Status**: TESTED (test-vectorized-integration-diagnostic.R)
- **Called by**: Multiple influence calculation functions
- **Recommendation**: **KEEP** - Core computation

#### 10. `compute_influence_terms()` - **REMOVAL CANDIDATE**
- **Status**: NOT TESTED
- **Defined in**: `R/influence.R`
- **Test usage**: 0 test files
- **Recommendation**: **INVESTIGATE** - May be deprecated or internal-only

### Utility Functions - **INTERNALIZED** ✅

The following utility functions have been marked as `@keywords internal` and removed from NAMESPACE:

#### 11. `add_terminal_observations()` - **INTERNALIZED** ✅
- **Status**: TESTED (test-PCORI_within_group_model.R)
- **Purpose**: Add terminal rows to data for intensity modeling
- **Action taken**: Marked as `@keywords internal`
- **Rationale**: Called internally by `prepare_SensIAT_data()`, not needed in public API

#### 12. `add_class()` - **INTERNALIZED** ✅
- **Status**: NOW TESTED (test-add_class.R with 16 assertions)
- **Called by**: `jackknife.R`, `prune.R` (internal use only)
- **Purpose**: Add S3 class to objects (S3 system only, not S4/R6)
- **Action taken**: 
  - Marked as `@keywords internal`
  - Added comprehensive test suite (16 test cases)
  - Updated documentation to clarify S3-only usage
- **Rationale**: Pure utility function, only called internally

#### 13. `extrapolate_from_last_observation()` - **INTERNALIZED** ✅
- **Status**: TESTED (3 test files, 47 passing assertions)
- **Called by**: `fit_SensIAT_marginal_mean_model_generalized.R`
- **Action taken**: Marked as `@keywords internal`
- **Rationale**: Implementation detail for marginal model fitting

#### 14. `extrapolate_from_last_observation_multiple()` - **INTERNALIZED** ✅
- **Status**: TESTED (test-extrapolate_from_last_observation.R)
- **Action taken**: Marked as `@keywords internal`
- **Rationale**: Batch version of extrapolate function, internal utility

#### 15. `make_term2_integrand_fast()` - **REMOVAL CANDIDATE**
- **Status**: NOT TESTED
- **Defined in**: `R/term2_fast.R`
- **Called by**: `compute_term2_influence.R` (only)
- **Test usage**: 0 test files
- **Recommendation**: **Mark as @keywords internal**
  - Only called by one internal function
  - No direct testing
  - Implementation detail

## Dependency Tree Analysis

### `fit_SensIAT_fulldata_model` Dependency Chain

```
fit_SensIAT_fulldata_model (EXPORTED, NOT TESTED)
├── fit_SensIAT_within_group_model (EXPORTED, TESTED)
    ├── prepare_SensIAT_data (EXPORTED, TESTED)
    ├── fit_SensIAT_marginal_mean_model (EXPORTED, NOT TESTED)
    │   └── [internal computations]
    ├── outcome_modeler (user-provided)
    │   ├── fit_SensIAT_single_index_fixed_coef_model (EXPORTED, TESTED)
    │   ├── fit_SensIAT_single_index_fixed_bandwidth_model (EXPORTED, NOT TESTED)
    │   └── fit_SensIAT_single_index_norm1coef_model (EXPORTED, TESTED)
    └── compute_influence_terms (EXPORTED, NOT TESTED)
        ├── compute_SensIAT_expected_values (EXPORTED, TESTED)
        └── make_term2_integrand_fast (EXPORTED, NOT TESTED)
```

### `jackknife` Dependency Chain

```
jackknife (EXPORTED, TESTED)
├── add_class (INTERNALIZED ✅, TESTED with 16 assertions)
└── [S3 methods for different model types]
```

## Recommendations Summary

### COMPLETED - Internalized ✅

**Utility functions (4 functions internalized):**
1. ✅ **`add_class()`** - Internalized, tests added (16 assertions)
2. ✅ **`add_terminal_observations()`** - Internalized, already tested
3. ✅ **`extrapolate_from_last_observation()`** - Internalized, already tested (47 assertions)
4. ✅ **`extrapolate_from_last_observation_multiple()`** - Internalized, already tested

**Impact**: These functions remain available via `SensIAT:::function_name()` for internal use but are no longer part of the public API. All tests pass (jackknife: 23 assertions, extrapolate: 47 assertions, add_class: 16 assertions).

### HIGH PRIORITY - Consider Removing

5. ✂️ **`fit_SensIAT_fulldata_model()`** - Thin wrapper, not tested, low value
6. ✂️ **`fit_SensIAT_single_index_fixed_bandwidth_model()`** - Not tested, only in examples
7. ✂️ **`compute_influence_terms()`** - Not tested, may be internal-only

### MEDIUM PRIORITY - Consider Internalizing

8. 🔒 **`make_term2_integrand_fast()`** - Implementation detail, only called by one function
9. 🔒 **`fit_SensIAT_marginal_mean_model()`** - Not tested directly, called internally

### KEEP - Core API or Well-Tested

10. ✅ All other functions are either core API or well-tested

### HIGH PRIORITY - Consider Removing

1. **`fit_SensIAT_fulldata_model()`** - Thin wrapper, not tested, low value
2. **`fit_SensIAT_single_index_fixed_bandwidth_model()`** - Not tested, only in examples
3. **`compute_influence_terms()`** - Not tested, may be internal-only

### MEDIUM PRIORITY - Consider Internalizing

4. **`add_class()`** - Pure utility, only called internally
5. **`make_term2_integrand_fast()`** - Implementation detail, only called by one function
6. **`fit_SensIAT_marginal_mean_model()`** - Not tested directly, called internally

### KEEP - Core API or Well-Tested

7. **`fit_SensIAT_within_group_model()`** - Main entry point
8. **`fit_SensIAT_single_index_fixed_coef_model()`** - Well tested
9. **`fit_SensIAT_single_index_norm1coef_model()`** - Well tested
10. **`jackknife()`** - Core statistical method
11. **`prepare_SensIAT_data()`** - Essential preprocessing
12. **`compute_SensIAT_expected_values()`** - Core computation
13. **`add_terminal_observations()`** - Useful utility
14. **`extrapolate_from_last_observation()`** - Tested utility
15. **`extrapolate_from_last_observation_multiple()`** - Tested utility

## Impact Analysis

### If `fit_SensIAT_fulldata_model()` is Removed:

**Affected files:**
- `R/fulldata_model.R` - Could be removed entirely or marked internal
- `R/autoplot.R` - Update documentation examples
- Documentation examples would need updating

**User impact:**
- Minimal - users can call `fit_SensIAT_within_group_model()` separately for each group
- Migration: Split single call into two calls with filter

**Code example for users:**
```r
# Old way (removed):
model <- fit_SensIAT_fulldata_model(
    data = my_data,
    trt = treatment_group == 'active',
    ...
)

# New way (recommended):
control_model <- fit_SensIAT_within_group_model(
    group.data = filter(my_data, treatment_group != 'active'),
    ...
)
treatment_model <- fit_SensIAT_within_group_model(
    group.data = filter(my_data, treatment_group == 'active'),
    ...
)
```

### If `fit_SensIAT_single_index_fixed_bandwidth_model()` is Removed:

**Affected files:**
- `R/sim_outcome_modeler.R` - Remove function definition
- `inst/examples/basic.R` - Update example

**User impact:**
- Low - only appears in one example file
- Alternative: Use `fit_SensIAT_single_index_fixed_coef_model()` instead

## Files That Could Be Completely Removed

1. **`R/fulldata_model.R`** - If `fit_SensIAT_fulldata_model()` removed
   - Also contains `predict.SensIAT_fulldata_model()`
   - Consider keeping predict method but internalizing

## Next Steps

1. **Verify** that identified functions are truly unused in user code
2. **Check vignettes** for usage of candidate functions
3. **Review** git history to understand original intent
4. **Create deprecation plan** if removing functions:
   - Mark as deprecated in current version
   - Add `.Deprecated()` warnings
   - Remove in next major version
5. **Update** documentation and examples
6. **Add tests** for functions being kept (if currently missing)

## Notes

- Functions marked as `@keywords internal` still appear in coverage but signal they're not part of the public API
- Some functions may be exported for flexibility even if not directly tested
- Consider user feedback before removing any exported functions
- Check CRAN download stats if package is on CRAN to gauge actual usage
