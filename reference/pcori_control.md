# Control Parameters for Fitting the within Group Model

Control Parameters for Fitting the within Group Model

## Usage

``` r
pcori_control(
  integration.method = c("quadvcpp", "quadv", "linear", "numerical", "piecewise"),
  intensity.bandwidth = NULL,
  resolution = 1000,
  resolution.within.period = 50,
  tol = .Machine$double.eps^(1/4),
  ...
)
```

## Arguments

- integration.method:

  Method for integration when computing the second influence term.

- intensity.bandwidth:

  The bandwidth for the intensity model.

- resolution:

  The number of points to use for numerical integration.

- resolution.within.period:

  The number of points to use for numerical integration within a period.

- tol:

  The tolerance for numerical integration.

- ...:

  Currently ignored.

## Value

a list of control parameters.
