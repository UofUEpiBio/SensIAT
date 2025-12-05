# Integrate function using adaptive Simpson quadrature.

Integrate function using adaptive Simpson quadrature.

## Usage

``` r
pcoriaccel_integrate_simp(integrand, lo, hi, tol = 1.490116e-08)
```

## Arguments

- integrand:

  The integrand, must take scalar argument, may return scalar, vector,
  or matrix.

- lo:

  Lower integration bound

- hi:

  Upper integration bound

- tol:

  Tolerance for integration, default `.Machine$double.eps^(1/2)`

## Value

integration result, list with elements `$Q` (the integral estimate),
`$fcnt` (the number of function evaluations), and `$estim.prec` (a
(pessimistic) estimate of the precision).
