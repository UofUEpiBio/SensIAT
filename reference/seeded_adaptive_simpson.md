# Seeded adaptive Simpson's rule

Adaptive Simpson's quadrature that starts with a pre-defined set of seed
points. The algorithm evaluates the function at seed points, then
recursively subdivides intervals where the error estimate exceeds the
tolerance.

## Usage

``` r
seeded_adaptive_simpson(f, seeds, tol = 1e-06, max_recursion = 8)
```

## Arguments

- f:

  Function to integrate (returns vector)

- seeds:

  Numeric vector of seed points (initial partition)

- tol:

  Error tolerance

- max_recursion:

  Maximum recursion depth

## Value

Integral approximation (vector)
