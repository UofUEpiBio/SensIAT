# Plot for Estimated Treatment Effect for `SensIAT_fulldata_jackknife_results` Objects

The horizontal and vertical axes represent the sensitivity parameter
`alpha` for the control and treatment groups, respectively. The plot
shows at each combination of `alpha` values zero if the 95% confidence
interval contains zero, otherwise the bound of the confidence interval
that is closest to zero.

## Usage

``` r
# S3 method for class 'SensIAT_fulldata_jackknife_results'
autoplot(object, ..., include.rugs = NA)
```

## Arguments

- object:

  A `SensIAT_fulldata_jackknife_results` object.

- ...:

  Additional arguments passed to `predict`.

- include.rugs:

  If `TRUE`, adds rugs to the plot. If `FALSE`, no rugs are added. When
  `NA`, rugs are added only if the number of distinct values of
  `alpha_control` and `alpha_treatment` is less than or equal to 10.

## Examples

``` r
# Note: fitting the jackknife is computationally expensive,
#       so this example is here for reference.
if (FALSE) { # \dontrun{
full.object <-
    fit_SensIAT_fulldata_model(
        data = SensIAT_example_fulldata,
        trt = Treatment_group == 'treatment',
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        knots = c(60, 260, 460),
        alpha = c(-0.6, -0.3, 0, 0.3, 0.6)
    )
jk.full.model <- jackknife(full.object, time = 180)
ggplot2::autoplot(jk.full.model)
} # }
```
