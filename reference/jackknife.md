# Perform Jackknife Resampling on an Object

Perform Jackknife Resampling on an Object

## Usage

``` r
jackknife(object, ...)

# S3 method for class 'SensIAT_within_group_model'
jackknife(object, time, ...)

# S3 method for class 'SensIAT_fulldata_model'
jackknife(object, time, ...)
```

## Arguments

- object:

  An object to cross validate on.

- ...:

  Additional arguments passed to the method.

- time:

  Time points for which to estimate the response.

## Value

A data frame of the jackknife resampling results. For `SensIAT` objects,
a `tibble` with columns alpha, time, jackknife_mean, and jackknife_var,
where jackknife_mean is the mean of the jackknife estimates and
jackknife_var is the estimated variances of the response at the given
time points for the specified alpha values.

## Methods (by class)

- `jackknife(SensIAT_within_group_model)`: Perform jackknife resampling
  on a `SensIAT_within_group_model` object.

- `jackknife(SensIAT_fulldata_model)`: Perform jackknife resampling on a
  `SensIAT_fulldata_model` object.

## Examples

``` r
if (FALSE) { # \dontrun{
object <-
    fit_SensIAT_within_group_model(
        group.data = SensIAT_example_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        intensity.args = list(bandwidth = 30),
        knots = c(60, 260, 460),
        End = 830
    )
jackknife.estimates <- jackknife(object, time = c(90, 180, 270, 360, 450))
} # }
```
