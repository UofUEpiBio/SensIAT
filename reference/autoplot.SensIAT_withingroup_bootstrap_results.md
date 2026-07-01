# Plot Bootstrap Marginal Mean Curves with Confidence Bands

Plot Bootstrap Marginal Mean Curves with Confidence Bands

## Usage

``` r
# S3 method for class 'SensIAT_withingroup_bootstrap_results'
autoplot(object, time = NULL, level = 0.95, ...)
```

## Arguments

- object:

  A `SensIAT_withingroup_bootstrap_results` object.

- time:

  Optional vector of times. If `NULL`, uses an internal regular grid.

- level:

  Confidence level used for bootstrap intervals.

- ...:

  Ignored.

## Value

A `ggplot2` object.
