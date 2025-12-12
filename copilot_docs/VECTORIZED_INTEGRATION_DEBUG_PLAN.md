# Plan to Fix Vectorized Integration Implementation

## Problem Summary
The vectorized integration produces results that differ by **orders of magnitude** from the reference implementation:
- Reference (quadv): ~100-260
- Vectorized: ~0-1
- Relative error: ~100%

## Root Cause Analysis

### Key Differences Between Implementations

#### 1. **Original Method (`compute_influence_term_2_quadv_sim`)**
```r
# Integrand at time t:
integrand = B(t) * conditional_mean(t)

where:
- B(t) = spline basis vector at time t
- conditional_mean(t) = E[Y*exp(αY)|X(t)] / E[exp(αY)|X(t)]
- X(t) is linearly interpolated between observation boundaries
- Integration is done piecewise between observed times
```

#### 2. **Vectorized Method (Current)**
```r
# Integrand at time t:
integrand = W(t) * (conditional_mean(t) - marginal_mean(t))

where:
- W(t) = V_inv %*% B(t) * exp(-μ(t))
- μ(t) = B(t)' %*% marginal_beta
- marginal_mean(t) = inv_link(μ(t))
```

### **Critical Issue**: The vectorized method uses a WEIGHTED integrand with marginal mean subtraction, but the original method does NOT!

## Step-by-Step Fix Plan

### Phase 1: Understand the Mathematical Formulation (1-2 hours)
1. ✅ Review the mathematical derivation in the paper/documentation
2. ✅ Identify what term2 actually represents:
   - It's the influence function component for the marginal mean model
   - Should integrate: `B(t) * E[Y(t)|X(t), α]` over time
3. ✅ Determine if weighting and marginal mean subtraction are correct
   - **HYPOTHESIS**: The test is using `marginal_beta = 0` which makes the weight incorrect
   - With beta=0, the marginal mean should be 0, but the weight calculation might be wrong

### Phase 2: Create Diagnostic Test (30 min)
**File**: `tests/testthat/test-vectorized-integration-diagnostic.R`

```r
test_that("diagnostic: compare integrand values", {
  # Setup same as main test
  # ... (copy setup code)
  
  # Test 1: Evaluate integrand at specific time points
  test_times <- c(100, 200, 300)
  
  for (t in test_times) {
    # Original method: compute integrand manually
    imputed <- impute_patient_df(t, df_i, variables, centering.statistics, TRUE)
    expected <- compute_SensIAT_expected_values(outcome.model, 0, imputed)
    B <- evaluate(base, t)
    original_integrand <- as.vector(B * (expected$E_Yexp_alphaY / expected$E_exp_alphaY))
    
    # Vectorized method: what it computes
    weight <- weight_fn(t)
    marginal <- marginal_mean_fn(t)
    vectorized_integrand <- weight * (conditional_mean - marginal)
    
    cat(sprintf("Time %d:\n", t))
    cat("  Original integrand:", original_integrand, "\n")
    cat("  Vectorized integrand:", vectorized_integrand, "\n")
    cat("  Weight:", weight, "\n")
    cat("  Marginal mean:", marginal, "\n")
  }
})
```

### Phase 3: Fix the Integrand Calculation (2-3 hours)

#### Option A: Simplify to Match Original
**Most likely fix**: Remove the weighting and marginal mean subtraction

**File**: `src/vectorized_integration.cpp`
```cpp
// BEFORE (lines 75-83):
for (int i = 0; i < n_alphas; ++i) {
    double conditional_mean = E_Yexp_alphaY[i] / E_exp_alphaY[i];
    double integrand_scalar = conditional_mean - marginal_mean;
    
    for (int j = 0; j < weight_length; ++j) {
        result(j, i) = weight[j] * integrand_scalar;
    }
}

// AFTER:
for (int i = 0; i < n_alphas; ++i) {
    double conditional_mean = E_Yexp_alphaY[i] / E_exp_alphaY[i];
    
    for (int j = 0; j < weight_length; ++j) {
        result(j, i) = weight[j] * conditional_mean;  // Just B(t) * E[Y|X]
    }
}
```

**AND update R wrapper** (`R/vectorized_integration.R`):
```r
# Weight function should just be the basis
weight_fn <- function(t) {
  as.vector(pcoriaccel_evaluate_basis(base, t))
}
```

#### Option B: Fix the Mathematical Formulation
If weighting IS correct (needs verification from paper):
- Debug why `marginal_beta = 0` produces wrong results
- Check if `V_inv` calculation is correct
- Verify `exp(-mu)` weighting is appropriate

### Phase 4: Fix Interpolation Method (1-2 hours)

The original uses **piecewise linear interpolation** of X(t):
```r
a <- (time - lower)/(upper-lower)
xb_time <- (1-a)*xb_lower + a*xb_upper
```

The vectorized method calls `impute_fn(t, patient_data)` which uses a different approach.

**Fix**: Make vectorized method use same interpolation
- Modify C++ to track integration intervals
- Interpolate between boundaries like original

### Phase 5: Fix Integration Algorithm (1 hour)

Original uses `pracma::quadv` (adaptive Lobatto quadrature).
Vectorized uses custom adaptive Simpson.

**Potential issues**:
1. Simpson's rule might need different tolerance
2. Initial point selection might differ
3. Convergence criteria might be too strict/loose

**Fix**:
- Match convergence tolerance: `tol=.Machine$double.eps^(1/4)` ≈ 1.22e-4
- Verify Simpson's rule implementation
- Test with simple known integrals

### Phase 6: Add Integration Tests (1 hour)

**File**: `tests/testthat/test-vectorized-integration-simple.R`

Add tests for:
1. **Constant function**: `∫ c dt` should equal `c * (b-a)`
2. **Linear function**: `∫ t dt` from a to b should equal `(b²-a²)/2`
3. **Basis functions**: `∫ B(t) dt` should match known values
4. **Known outcome model**: Test with simple linear model where answer is analytical

### Phase 7: Validate Piecewise Integration (30 min)

Original integrates **piecewise** between observation times.
Vectorized might integrate whole interval at once.

**Check**: Does the test properly break integration into periods?
- Original: integrates between each pair of observation times
- Vectorized: should do the same

## Implementation Order

1. **Create diagnostic test** (Phase 2) - Highest priority
2. **Run diagnostics** - Identify exact mismatch
3. **Fix integrand** (Phase 3 Option A) - Most likely issue
4. **Add simple integration tests** (Phase 6) - Verify fix works
5. **Fix interpolation if needed** (Phase 4)
6. **Fine-tune integration algorithm** (Phase 5)
7. **Validate piecewise** (Phase 7)

## Success Criteria

- ✅ Diagnostic test shows matching integrand values at sample points
- ✅ Simple integration tests pass (constant, linear functions)
- ✅ Vectorized result matches original within tolerance (1e-4 relative error)
- ✅ Test passes for multiple alpha values (not just alpha=0)
- ✅ Performance is better than original (vectorized should be faster)

## Estimated Time

- **Minimum**: 3-4 hours (if Option A fix works immediately)
- **Expected**: 6-8 hours (with debugging and validation)
- **Maximum**: 12+ hours (if mathematical formulation needs rework)

## Next Steps

1. Create the diagnostic test to pinpoint the exact issue
2. Based on diagnostic output, choose Option A or B
3. Implement fix
4. Validate with simple tests
5. Re-enable main test

---

**Note**: The most likely issue is that the vectorized method incorrectly applies weighting and marginal mean subtraction that shouldn't be there for term2 calculation. Start with Option A.
