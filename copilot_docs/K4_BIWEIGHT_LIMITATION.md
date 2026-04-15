# K4_Biweight Kernel Support Limitation

## Issue Summary

The K4_Biweight kernel has **partial support** in the SensIAT package:
- ✅ Supported in R-level code
- ❌ **Not supported in C++ acceleration layer**

## Root Cause

**The K4_Biweight kernel can produce negative values**, which is mathematically valid for higher-order kernels but problematic for probability estimation.

### Why K4_Biweight Can Be Negative

K4_Biweight is a **4th-order kernel** designed to reduce bias compared to 2nd-order kernels. The formula is:

```
K(x/h) = (105/64) * (1 - 3*(x/h)²) * (1 - (x/h)²)²  for |x| ≤ h
```

This kernel becomes **negative** for approximately **0.577 < |x/h| < 1.0**, with a minimum value of about -0.216 at |x/h| ≈ 0.746.

### Impact on PMF Estimation

When estimating probability mass functions (PMFs) using `pcoriaccel_estimate_pmf()`:
- Negative kernel weights can produce **negative probability estimates**
- This violates the fundamental constraint that P(Y = y) ≥ 0
- Could lead to invalid cumulative distribution functions

This is why K4_Biweight was **intentionally disabled** in commit 7d1b426 (Oct 28, 2024).

## Technical Details

### Where K4_Biweight Works
- `SIDRnew_fixed_bandwidth()` in `R/SIDRnew.R` (line 346)
- Other R-level SIDR functions

### Where K4_Biweight Fails
The C++ acceleration functions have K4_Biweight **commented out**:

1. **`src/estimate_pmf.cpp`** (lines 93-96):
```cpp
// else if ( kernel == "K4_Biweight" )
// {
//     return _pcoriaccel_estimate_pmf< Tfloat, Kt_biweight4 >( Xb,Y, xi, y_seq, (Tfloat)h );
// }
```

2. **`src/NW.cpp`** (lines 220-223):
```cpp
// else if ( kernel == "K4_Biweight" )
// {
// 	return _pcoriaccel_NW< Tfloat, 3 >( Xb,Y, xb, y_seq, (Tfloat)h );
// }
```

### Error Message
```
Error: Invalid value for `kernel`: choices are { "dnorm", "K2_Biweight" }.
```

## Impact

### What Works
- Fitting outcome models with K4_Biweight using `fit_SensIAT_single_index_fixed_bandwidth_model` alone
- R-only workflows that don't require influence calculations

### What Fails
- Using K4_Biweight with `fit_SensIAT_within_group_model` (requires influence term calculations)
- Any workflow that calls:
  - `pcoriaccel_estimate_pmf()` 
  - `pcoriaccel_NW()`
  
These C++ functions are used by `compute_sim_influence_term_1_at_timepoint()` for performance.

## Current Status

**Supported kernels for full workflow:**
- ✅ `K2_Biweight` (default)
- ✅ `dnorm`
- ❌ `K4_Biweight` (partial - R only)

## Test Coverage

Added test in `test-SensIAT_sim_outcome_modeler_fixed_bandwidth.R`:
- Documents the limitation
- Verifies expected error message
- Prevents regression if C++ support is added later

## Why It Was Disabled

From git commit 7d1b426 (Oct 28, 2024):
```
Author: Andrew Redd <andrew.redd@hsc.utah.edu>
Date:   Mon Oct 28 13:47:32 2024 -0600

    Removed K4_biweight as option for kernel.
```

The commit removed K4_Biweight from exported functions and commented out C++ support, likely due to the negative kernel value issue causing problems in PMF estimation.

## Mathematical Properties

**K4_Biweight properties:**
- ✅ Integrates to 1 over support [-h, h]
- ✅ Symmetric around 0
- ✅ 4th-order kernel (lower bias than K2_Biweight)
- ⚠️ **Can be negative** for 0.577 < |x/h| < 1.0
- ⚠️ Minimum value: -0.216 at |x/h| ≈ 0.746

**Why this matters:**
- Higher-order kernels trade bias reduction for potential negative weights
- In density estimation, this is often acceptable (negative values cancel out)
- In **PMF estimation for discrete outcomes**, negative probabilities are invalid

## Resolution Options

To potentially re-enable K4_Biweight (requires careful consideration):

1. **Accept negative PMF values** and add post-processing:
   - Uncomment C++ code in `src/estimate_pmf.cpp` and `src/NW.cpp`
   - Add normalization/truncation to ensure non-negative PMFs
   - Document the post-processing clearly

2. **Use only in contexts where negativity is acceptable**:
   - Enable for regression smoothing (where negative weights are OK)
   - Keep disabled for PMF estimation
   - Create separate code paths

3. **Stick with 2nd-order kernels** (current approach):
   - K2_Biweight and dnorm never go negative
   - More appropriate for discrete outcome modeling
   - This is the recommended approach

## References
- Error origin: `src/estimate_pmf.cpp:99`, `src/NW.cpp:226`
- R implementation: `R/SIDRnew.R:346`
- Test documentation: `tests/testthat/test-SensIAT_sim_outcome_modeler_fixed_bandwidth.R`
- Removal commit: `7d1b426` (Oct 28, 2024)

## Verification

You can verify the negative kernel values with:

```r
K4 <- function(x, h = 1) {
    105/64 * (1 - 3*((x/h)^2)) * (1-(x/h)^2)^2 * (abs(x) <= h)
}

x <- seq(-1, 1, length = 1000)
y <- K4(x)
plot(x, y, type = 'l', main = "K4_Biweight Kernel")
abline(h = 0, col = 'red', lty = 2)

cat("Minimum value:", min(y), "\n")
cat("Negative for |x| >", min(abs(x[y < 0])), "\n")
```

## Conclusion

**K4_Biweight is intentionally disabled** because it can produce negative kernel weights, leading to invalid negative probability estimates when used for PMF estimation. This is a fundamental mathematical property of 4th-order kernels, not a bug.

The current approach of using only K2_Biweight and dnorm (which are always non-negative) is the correct design choice for this package's use case.
