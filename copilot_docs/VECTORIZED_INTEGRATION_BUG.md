# Vectorized Integration Bug

## Summary

The experimental vectorized C++ integration in `src/vectorized_integration.cpp` has a critical bug that causes it to return incorrect results (typically 10x too small).

## Root Cause

The adaptive Simpson implementation uses a **global convergence flag** that marks an alpha value as "converged" as soon as ANY segment converges. This causes subsequent segments to be skipped entirely.

### Code Flow

1. The integration interval [a, b] is initially split into 3 pieces
2. `adaptive_helper` is called for the first piece [x0, x2]
3. If that segment converges, it:
   - Adds the segment's integral to `final_integrals`
   - Sets `integration_state.alpha_states[alpha_idx].converged = true`
4. When `adaptive_helper` is called for the second piece [x2, x4], it immediately skips:
   ```cpp
   if (integration_state.alpha_states[alpha_idx].converged) {
       alpha_converged[alpha_idx] = true;
       continue;  // BUG: Skips this entire segment!
   }
   ```
5. Third piece also skipped
6. Result: Only ~1/3 of the interval is integrated

### Evidence

Test case `test-integration-detailed.R` shows:
- **Integrand at t=137**: Perfect match (ratio = 1.0)
- **Integral over [60, 214]**:
  - pracma::quadv: 109.47
  - Trapezoidal (1000 points): 109.47
  - C++ vectorized: **10.57** (10x too small!)
- **Function evaluations**: Only 25 (should be 100s for high accuracy)
- **Convergence**: Reports TRUE despite wrong answer

## Required Fix

The code needs to be rewritten to use **per-segment convergence** instead of global convergence:

1. Remove the global `integration_state.alpha_states[alpha_idx].converged` flag
2. Make `adaptive_helper` return the integral for its segment
3. Let the recursion naturally sum up leaf segment contributions
4. Only check global convergence after processing all initial segments

### Pseudocode for Correct Approach

```cpp
NumericVector adaptive_simpson(a, c, e, depth) {
    NumericVector Q1 = simpson_3point(a, c, e);
    NumericVector Q2 = simpson_5point(a, b, c, d, e);
    
    if (max_abs(Q2 - Q1) < tolerance * max_abs(Q2)) {
        // Segment converged - return its integral
        return Q2 + (Q2 - Q1) / 15.0;  // Romberg extrapolation
    } else {
        // Segment not converged - split and recurse
        return adaptive_simpson(a, b, c, depth+1) + 
               adaptive_simpson(c, d, e, depth+1);
    }
}
```

## Workaround

Until this is fixed, use the original `compute_influence_term_2_quadv_sim()` method which uses `pracma::quadv`. It is slower but correct.

## Testing

The test suite now includes:
- `test-vectorized-integration-simple.R`: End-to-end comparison (currently FAILS)
- `test-vectorized-integration-diagnostic.R`: Integrand component verification (PASSES)
- `test-integration-detailed.R`: Single-piece detailed analysis (shows the bug clearly)

All tests confirm:
- ✅ Piecewise interval setup is correct
- ✅ Linear interpolation of Xβ is correct  
- ✅ PMF calculation is correct
- ✅ Integrand formula is correct
- ❌ Integration algorithm has fundamental bug

## Impact

This bug affects any code using `compute_term2_influence_vectorized()`. The experimental vectorized path should be disabled until this is fixed.

## Next Steps

1. Rewrite C++ adaptive Simpson to use return-based accumulation
2. Add unit tests for the C++ integrator with known test functions
3. Verify against pracma::quadv with multiple test cases
4. Performance benchmark after correctness is established
