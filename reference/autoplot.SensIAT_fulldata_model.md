# Plot for Estimated Treatment Effect for `SensIAT_fulldata_model` Objects

The horizontal and vertical axes represent the sensitivity parameter
`alpha` for the control and treatment groups, respectively. The contour
plot shows the estimated treatment effect at each combination of `alpha`
values.

## Usage

``` r
# S3 method for class 'SensIAT_fulldata_model'
autoplot(object, time, include.rugs = NA, ...)
```

## Arguments

- object:

  A `SensIAT_fulldata_model` object.

- time:

  Time at which to plot the estimates.

- include.rugs:

  If `TRUE`, adds rugs indicating the locations where the sensitivity
  was evaluated to the plot. If `FALSE`, no rugs are added. If `NA`,
  rugs are added only if the number of distinct values of
  `alpha_control` and `alpha_treatment` is less than or equal to 10.

- ...:

  Additional arguments passed to `predict`.

## Value

A `ggplot2` object.

## Examples

``` r
# \donttest{
full.object <-
    fit_SensIAT_fulldata_model(
        data = SensIAT_example_fulldata,
        trt = Treatment_group == "treatment",
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        knots = c(60, 260, 460),
        alpha = c(-0.6, -0.3, 0, 0.3, 0.6)
    )
ggplot2::autoplot(full.object, time = 180)

# }
```
