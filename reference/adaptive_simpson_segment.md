# Recursive adaptive Simpson for a single interval

Recursive adaptive Simpson for a single interval

## Usage

``` r
adaptive_simpson_segment(f, a, b, fa, fb, tol, max_recursion, recursion_level)
```

## Arguments

- f:

  Function to integrate

- a:

  Left endpoint

- b:

  Right endpoint

- fa:

  Function value at a

- fb:

  Function value at b

- tol:

  Tolerance for this interval

- max_recursion:

  Maximum recursion depth

- recursion_level:

  Current recursion level

## Value

Integral approximation
