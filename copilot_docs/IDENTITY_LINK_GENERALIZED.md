# Identity Link Implementation in Generalized Version

## Summary of Changes

Modified `fit_SensIAT_marginal_mean_model_generalized()` to implement the identity link case using the same BB::sane optimization infrastructure as the other link functions (log and logit), rather than delegating to the linear version.

## What Changed

### Removed (Lines 115-128 in original)
```r
if (link == "identity") {
    return(
        fit_SensIAT_marginal_mean_model(
            data = data,
            id = !!id,
            alpha = alpha,
            knots = knots,
            outcome.model = outcome.model,
            intensity.model = intensity.model,
            spline.degree = spline.degree,
            ...
        )
    )
}
```

### Added (in `loss == "lp_mse"` branch for identity link)
```r
if (link == "identity") {
    # link.fun <- function(mu) mu
    inv.link <- function(eta) eta
    # d1.inv.link <- function(eta) rep(1, length(eta))
    V <- GramMatrix(base)
    V.inv <- solve(V)

    W <- function(t, beta) {
        B <- pcoriaccel_evaluate_basis(base, t)
        # For identity link: ds/dz = 1, so weight function is constant
        as.vector(V.inv %*% B)
    }
}
```

## Architecture Impact

### Before
- **identity link**: Delegated to linear version via early return
- **log/logit links**: Full BB::sane optimization in generalized version
- **Result**: Different code paths, different optimization strategies

### After
- **All links** (identity, log, logit): Use the same BB::sane optimization infrastructure
- **Unified codebase**: Single optimization loop for all link functions
- **Easier comparison**: Can now directly compare performance and results across links

## Key Implementation Details

### Identity Link Properties
- **Inverse link**: `η → η` (identity function)
- **Derivative**: `ds/dη = 1` (constant)
- **Weight function W(t, β)**: Constant, independent of β
  - For lp_mse loss: `W(t, β) = V⁻¹ B(t)`
  - For quasi-likelihood loss: Same (already implemented)

### Comparison with Other Links

| Aspect | Identity | Log | Logit |
|--------|----------|-----|-------|
| **Inverse Link** | `η` | `exp(η)` | `exp(η)/(1+exp(η))` |
| **Derivative** | 1 | `exp(η)` | `exp(η)/(1+exp(η))²` |
| **lp_mse W(t,β)** | `V⁻¹B(t)` | `(V⁻¹B(t))·exp(-μ)` | `(V⁻¹B(t))·(exp(η)+2+exp(-η))` |
| **quasi-likelihood W(t,β)** | `V⁻¹B(t)` | Requires integration | Requires integration |

## Testing Plan

After this change, we should:

1. **Verify equivalence**: Compare results of identity link in generalized version with linear version
2. **Performance comparison**: Benchmark the generalized version identity implementation vs. linear version
3. **Code validation**: Ensure term1 and term2 computations produce identical results
4. **Infrastructure consistency**: Verify that all three link functions work correctly within unified framework

## Files Modified
- `R/fit_SensIAT_marginal_mean_model_generalized.R`: Removed delegation, added identity link implementation

## Files NOT Modified (as requested)
- `R/fit_sensiat_marginal_mean_model.R`: Left completely untouched
