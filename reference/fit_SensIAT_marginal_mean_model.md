# Fit the Marginal Means Model

Fit the Marginal Means Model

## Usage

``` r
fit_SensIAT_marginal_mean_model(
  data,
  id,
  alpha,
  knots,
  outcome.model,
  intensity.model,
  spline.degree = 3L,
  ...
)
```

## Arguments

- data:

  Data for evaluation of the model. Should match the data used to fit
  the intensity and outcome models.

- id:

  The subject identifier variable in the data. Lazy evaluation is used,
  so it can be a symbol or a string.

- alpha:

  Sensitivity parameter, a vector of values.

- knots:

  Location of spline knots. If a `SplineBasis` object is provided, it is
  used directly.

- outcome.model:

  The observed effects model.

- intensity.model:

  The assessment time intensity model.

- spline.degree:

  The degree of the spline basis, default is 3 (cubic splines).

- ...:

  Additional arguments passed to `compute_influence_terms`.

## Value

a list with the fitted model, including the coefficients and their
variances for each alpha value.

## Examples

``` r
# Note: example takes approximately 30 seconds to run.
# \donttest{
library(survival)
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
library(splines)
# Create followup data with lags
# added variables `..prev_time..`, `..delta_time..` and `..prev_outcome..`
# have special interpretations when computing the influence.
data_with_lags <- SensIAT_example_data |>
    dplyr::group_by(Subject_ID) |>
    dplyr::mutate(
        ..prev_outcome.. = dplyr::lag(Outcome, default = NA_real_, order_by = Time),
        ..prev_time.. = dplyr::lag(Time, default = 0, order_by = Time),
        ..delta_time.. = Time - dplyr::lag(.data$Time, default = NA_real_, order_by = Time)
    )

# Create the observation time intensity model
intensity.model <-
    coxph(Surv(..prev_time.., Time, !is.na(Outcome)) ~ ..prev_outcome.. + strata(Visit),
        data = data_with_lags |> dplyr::filter(.data$Time > 0)
    )

# Create the observed outcome model
outcome.model <-
    fit_SensIAT_single_index_fixed_coef_model(
        Outcome ~ ns(..prev_outcome.., df = 3) + ..delta_time.. - 1,
        id = Subject_ID,
        data = data_with_lags |> filter(Time > 0)
    )

# Fit the marginal outcome model
mm <- fit_SensIAT_marginal_mean_model(
    data = data_with_lags,
    id = Subject_ID,
    alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
    knots = c(60, 260, 460),
    intensity.model = intensity.model,
    time.vars = c("..delta_time.."),
    outcome.model = outcome.model
)
# }
```
