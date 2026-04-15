# Diagnosis: Anomalous Iteration Error in fit_SensIAT_marginal_mean_model_generalized

## Problem

Tests are failing with:

    Error: Optimization did not converge for alpha=0. Message: Anomalous iteration

## Root Cause

The issue occurs when using **identity link** with BB::sane
optimization. The diagnostic output shows:

    [Call 1] ||F(beta)|| = 1141.871 range: [ -22.47206 , 847.8221 ]
    [Call 2] ||F(beta)|| = 1141.869 range: [ -22.47204 , 847.8212 ]
    ...

The influence function has very large magnitude (~1142) and barely
decreases during optimization.

### Why This Happens

1.  **For identity link, the estimating equation is LINEAR in beta**

    - Weight function: `W(t, beta) = V^{-1} B(t)` (doesn’t depend on
      beta)
    - The influence equation reduces to:
      `sum_i sum_j W(t_ij) * (Y_ij - B(t_ij)'beta) / lambda(t_ij) = 0`
    - This is a linear system that can be solved directly!

2.  **The original `fit_SensIAT_marginal_mean_model` solves this
    analytically**:

    ``` r
    uncorrected.beta_hat <- (colSums(IT$term1) + colSums(IT$term2)) / length(IT$id)
    estimate <- as.vector(V_inverse %*% uncorrected.beta_hat)
    ```

3.  **BB::sane is designed for nonlinear systems**

    - Using it on a linear system is inefficient and can fail
    - The large influence function values suggest the system is poorly
      conditioned

## Proposed Solutions

### Solution 1: Use Direct Solution for Identity Link (RECOMMENDED)

For identity link (regardless of loss function), use the analytic
solution:

``` r
if (link == "identity") {
    # Identity link has linear influence equation - solve directly
    # Compute influence for all patients at beta = 0 (or any value - it's linear!)
    influence_by_patient <- map(unique_ids, compute_influence_by_patient, beta = rep(0, ncol(base)))
    
    # Sum influence terms
    total_influence <- reduce(map(influence_by_patient, getElement, 'total'), `+`)
    
    # Solve: V^{-1} * total_influence = 0
    # => beta = solve(V, total_influence)
    V <- GramMatrix(base)
    solution_par <- solve(V, total_influence)
    
    solution <- list(
        par = solution_par,
        convergence = 0,
        message = "Analytic solution (identity link)"
    )
} else {
    # Non-identity links require iterative optimization
    solution <- BB::sane(
        par = rep(0, ncol(base)),
        fn = influence,
        control = BBsolve.control
    )
}
```

### Solution 2: Better Initial Values for BB::sane

If we want to keep using BB::sane for all links, use better starting
values:

``` r
# Use the "naive" estimate as starting point
naive_influence_by_patient <- map(unique_ids, compute_influence_by_patient, beta = rep(0, ncol(base)))
naive_total <- reduce(map(naive_influence_by_patient, getElement, 'total'), `+`)
V <- GramMatrix(base)
initial_beta <- solve(V, naive_total)

solution <- BB::sane(
    par = initial_beta,  # Better starting point
    fn = influence,
    control = BBsolve.control
)
```

### Solution 3: Increase Tolerance / Iterations

Relax convergence criteria (NOT RECOMMENDED - doesn’t fix root cause):

``` r
BBsolve.control = list(maxit = 200, tol = 1e-2)
```

## Recommended Action

Implement **Solution 1** - it’s mathematically correct, faster, and more
robust. The identity link case should never need iterative optimization
since the problem is linear.

This is exactly what the original `fit_SensIAT_marginal_mean_model`
does, and we should maintain that behavior in the generalized version.
