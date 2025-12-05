# Plot Estimates at Given Times for `SensIAT_withingroup_jackknife_results` Objects

Horizontal axis represents time, and the vertical axis represents the
outcome from the model. Point plotted is the mean estimate, and the
error bars show the 95% confidence interval using the variance estimated
from the jackknife.

## Usage

``` r
# S3 method for class 'SensIAT_withingroup_jackknife_results'
autoplot(object, width = NULL, ...)
```

## Arguments

- object:

  A `SensIAT_withingroup_jackknife_results` object produced from
  `SensIAT_jackknife`.

- width:

  Width of the dodge for position, default is half the minimum distance
  between time evaluation points.

- ...:

  Ignored.

## Value

A `ggplot2` object.

## Examples

``` r
# Note: fitting the jackknife is computationally expensive,
#       so this example is here for reference.
if (FALSE) { # \dontrun{
fitted <-
fit_SensIAT_within_group_model(
    group.data = SensIAT_example_data,
    outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
    alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
    id = Subject_ID,
    outcome = Outcome,
    time = Time,
    intensity.args=list(bandwidth = 30),
    knots = c(60,260,460),
    End = 830
)
jackknife.estimates <- SensIAT_jackknife(fitted, time = c(90, 180, 270, 360, 450))
ggplot2::autoplot(jackknife.estimates)
} # }
```
