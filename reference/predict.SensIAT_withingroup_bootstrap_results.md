# Predict Marginal Mean from Bootstrap Coefficient Replicates

Predict Marginal Mean from Bootstrap Coefficient Replicates

## Usage

``` r
# S3 method for class 'SensIAT_withingroup_bootstrap_results'
predict(object, time, include.var = TRUE, level = 0.95, ...)
```

## Arguments

- object:

  A `SensIAT_withingroup_bootstrap_results` object.

- time:

  Time points of interest.

- include.var:

  Logical. If `TRUE` and link is identity, include the original-model
  asymptotic variance column `var`.

- level:

  Confidence level used to produce `lower` and `upper` bounds.

- ...:

  Currently ignored.

## Value

A `tibble` with bootstrap summaries by `alpha` and `time`.
