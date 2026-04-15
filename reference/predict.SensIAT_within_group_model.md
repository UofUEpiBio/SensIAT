# Give the Marginal Mean Estimate and its Estimated Asymptotic Variance

Give the marginal mean model estimate

## Usage

``` r
# S3 method for class 'SensIAT_fulldata_model'
predict(object, time, ...)

# S3 method for class 'SensIAT_within_group_model'
predict(object, time, include.var = TRUE, ..., base = object$base)
```

## Arguments

- object:

  SensIAT_within_group_model object

- time:

  Time points of interest

- ...:

  Currently ignored.

- include.var:

  Logical. If TRUE, the variance of the outcome is also returned

- base:

  A `SplineBasis` object used to evaluate the basis functions.

## Value

If include.var is TRUE, a `tibble` with columns time, mean, and var is
returned. otherwise if include.var is FALSE, only the mean vector is
returned.

## Functions

- `predict(SensIAT_fulldata_model)`: For each combination of `time` and
  `alpha` estimate the mean response and variance for each group as well
  as estimate the mean treatment effect and variance.

## Examples

``` r
# \donttest{
model <-
    fit_SensIAT_within_group_model(
        group.data = SensIAT_example_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        End = 830,
        knots = c(60, 260, 460),
    )
predict(model, time = c(90, 180))
#> # A tibble: 10 × 4
#>    alpha  time  mean     var
#>    <dbl> <dbl> <dbl>   <dbl>
#>  1  -0.6    90  1.47 0.00983
#>  2  -0.6   180  1.58 0.00654
#>  3  -0.3    90  1.68 0.00996
#>  4  -0.3   180  1.84 0.00699
#>  5   0      90  1.95 0.0105 
#>  6   0     180  2.13 0.00858
#>  7   0.3    90  2.28 0.0120 
#>  8   0.3   180  2.47 0.0123 
#>  9   0.6    90  2.65 0.0130 
#> 10   0.6   180  2.83 0.0197 
# }
```
