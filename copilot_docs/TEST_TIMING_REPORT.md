# Test Timing Report

**Generated:** Branch `20-cut-down-time-on-testing`  
**Total Test Time:** 443.4 seconds (~7.4 minutes)

## Summary

- **Total Test Files:** 32
- **Tests > 30 seconds:** 3 (flagged for refactoring)
- **Tests with `skip_on_cran`:** 8 files (partially or fully skipped on CI)

## ⚠️ Tests Flagged for Refactoring (> 30 seconds)

| File | Time (sec) | Tests | Passed | Skipped | Runs on CI? |
|------|-----------|-------|--------|---------|-------------|
| `test-term2-integration-methods.R` | **153.8** | 4 | 26 | 0 | ❌ No (skip_on_cran) |
| `test-term2-debug-isolation.R` | **58.6** | 2 | 2 | 1 | ❌ No (skip_on_cran) |
| `test-fit_sensiat_marginal_mean_model_generalized-lp_mse.R` | **30.4** | 5 | 14 | 3 | ✅ Yes |

### Key Observations:
1. **test-term2-integration-methods.R** is the slowest at 2.5+ minutes, but already skipped on CRAN
2. **test-term2-debug-isolation.R** is debug/diagnostic in nature and already skipped on CRAN
3. **test-fit_sensiat_marginal_mean_model_generalized-lp_mse.R** runs on CI and takes 30+ seconds - **priority for optimization**

---

## Complete Test Timing (Sorted by Time)

| Rank | File | Time (s) | Tests | Passed | Skipped | Runs on CI? |
|------|------|----------|-------|--------|---------|-------------|
| 1 | test-term2-integration-methods.R | 153.80 | 4 | 26 | 0 | ❌ |
| 2 | test-term2-debug-isolation.R | 58.56 | 2 | 2 | 1 | ❌ |
| 3 | test-fit_sensiat_marginal_mean_model_generalized-lp_mse.R | 30.35 | 5 | 14 | 3 | ✅ |
| 4 | test-PCORI_within_group_model.R | 27.73 | 2 | 5 | 0 | ✅ |
| 5 | test-sim_outcome_modeler.R | 13.11 | 1 | 21 | 0 | ✅ |
| 6 | test-numerical-integration-methods.R | 12.24 | 4 | 26 | 0 | ✅ |
| 7 | test-jackknife-fast.R | 11.46 | 4 | 8 | 0 | ✅ |
| 8 | test-jackknife.R | 9.35 | 3 | 3 | 0 | ❌ |
| 9 | test-SensIAT_sim_outcome_modeler_fixed_bandwidth.R | 8.73 | 10 | 148 | 0 | ✅ |
| 10 | test-vectorized-integration-simple.R | 7.81 | 2 | 6 | 0 | Partial ❌ |
| 11 | test-influence_term1.R | 7.43 | 1 | 1 | 0 | ✅ |
| 12 | test-fit_sensiat_marginal_mean_model_generalized-identity.R | 7.16 | 3 | 7 | 0 | ✅ |
| 13 | test-influence_term2.R | 6.95 | 1 | 1 | 0 | ✅ |
| 14 | test-glm_outcome_models.R | 5.60 | 9 | 41 | 2 | ✅ |
| 15 | test-SensIAT_sim_outcome_modeler_mave.R | 5.35 | 3 | 31 | 0 | ✅ |
| 16 | test-extrapolate_from_last_observation.R | 5.08 | 10 | 47 | 0 | ✅ |
| 17 | test-add_class.R | 5.04 | 7 | 16 | 0 | ✅ |
| 18 | test-estimate_baseline_intensity.R | 5.03 | 2 | 3 | 0 | ✅ |
| 19 | test-simulate_SensIAT_data.R | 5.02 | 15 | 39 | 0 | ✅ |
| 20 | test-jackknife-comprehensive.R | 5.01 | 3 | 0 | 3 | ❌ |
| 21 | test-fit_sensiat_marginal_mean_model_generalized-quasi_likelihood.R | 4.76 | 3 | 5 | 2 | ✅ |
| 22 | test-integration-detailed.R | 4.75 | 1 | 0 | 1 | ✅ |
| 23 | test-fit_sensiat_marginal_mean_model_generalized-performance.R | 4.66 | 2 | 0 | 2 | ❌ |
| 24 | test-fit_sensiat_marginal_mean_model.R | 4.57 | 0 | 0 | 0 | ✅ |
| 25 | test-common.R | 4.57 | 3 | 0 | 0 | ✅ |
| 26 | test-vectorized-integration.R | 4.57 | 4 | 3 | 4 | ❌ |
| 27 | test-influence_term2-fixed.R | 4.27 | 1 | 0 | 4 | ✅ |
| 28 | test-term2_performance.R | 4.20 | 1 | 5 | 0 | ✅ |
| 29 | test-vectorized-integration-diagnostic.R | 4.17 | 1 | 0 | 1 | ❌ |
| 30 | test-splinebasis.R | 4.12 | 1 | 14 | 0 | ✅ |
| 31 | test-util-match.names.R | 4.05 | 1 | 0 | 1 | ✅ |
| 32 | test-NW.R | 3.86 | 2 | 12 | 0 | ✅ |

---

## Files with `skip_on_cran` Directives

These files have tests that are skipped when running on CRAN/CI:

| File | Skip Coverage |
|------|---------------|
| test-term2-integration-methods.R | All tests skipped |
| test-term2-debug-isolation.R | All tests skipped |
| test-vectorized-integration.R | All tests skipped |
| test-vectorized-integration-diagnostic.R | All tests skipped |
| test-vectorized-integration-simple.R | Partial (1 test skipped) |
| test-jackknife-comprehensive.R | All tests skipped |
| test-jackknife.R | All tests skipped |
| test-fit_sensiat_marginal_mean_model_generalized-performance.R | All tests skipped |

Note: `test-jackknife-fast.R` has `skip_on_cran()` commented out - tests run on CI.

---

## Recommendations for Optimization

### Priority 1: `test-fit_sensiat_marginal_mean_model_generalized-lp_mse.R` (30.4s)
- **Runs on CI** - directly impacts build time
- Consider reducing sample size in test data
- Reduce iterations in optimization routines for tests
- Could add `skip_on_cran()` for expensive tests

### Priority 2: `test-PCORI_within_group_model.R` (27.7s)
- **Runs on CI** - second largest CI impact
- Review test complexity and data size
- Consider extracting expensive tests with `skip_on_cran()`

### Priority 3: `test-sim_outcome_modeler.R` (13.1s)
- Consider reducing MAVE iterations for tests
- Use smaller test datasets

### General Strategies:
1. Use `testthat::skip_on_cran()` for comprehensive/slow tests
2. Create "fast" variants of expensive tests with smaller data
3. Mock expensive operations where feasible
4. Reduce tolerances/iterations for CI tests

---

## CI Impact Summary

| Category | Count | Total Time |
|----------|-------|------------|
| Tests that run on CI | ~24 files | ~220s |
| Tests skipped on CI | ~8 files | ~220s |
| **Effective CI test time** | - | **~220 seconds** |
