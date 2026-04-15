# Outcome Modeler for `SensIAT` Single Index Model.

Outcome Modeler for `SensIAT` Single Index Model.

## Usage

``` r
fit_SensIAT_single_index_fixed_coef_model(
  formula,
  data,
  kernel = "K2_Biweight",
  method = "nmk",
  id = ..id..,
  initial = NULL,
  ...
)

fit_SensIAT_single_index_fixed_bandwidth_model(
  formula,
  data,
  kernel = "K2_Biweight",
  method = "nmk",
  id = ..id..,
  initial = NULL,
  ...
)
```

## Arguments

- formula:

  The outcome model formula

- data:

  The data to fit the outcome model to. Should only include follow-up
  data, i.e. time \> 0.

- kernel:

  The kernel to use for the outcome model.

- method:

  The optimization method to use for the outcome model, either
  `"optim"`, `"nlminb"`, or `"nmk"`.

- id:

  The patient identifier variable for the data.

- initial:

  Either a vector of initial values or a function to estimate initial
  values. If NULL (default), the initial values are estimated using the
  [`MAVE::mave.compute`](https://rdrr.io/pkg/MAVE/man/mave.html)
  function.

- ...:

  Currently ignored, included for future compatibility.

## Value

Object of class `SensIAT::Single-index-outcome-model` which contains the
outcome model portion.

## Functions

- `fit_SensIAT_single_index_fixed_bandwidth_model()`: for fitting with a
  fixed bandwidth

## Examples

``` r
# \donttest{
# A basic example using fixed intensity bandwidth.
object <-
    fit_SensIAT_within_group_model(
        group.data = SensIAT_example_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_bandwidth_model,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        knots = c(60, 260, 460),
        End = 830,
        intensity.args = list(bandwidth = 30)
    )

# A basic example using variable bandwidth but with fixed first coefficient.
object.bw <-
    fit_SensIAT_within_group_model(
        group.data = SensIAT_example_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        knots = c(60, 260, 460),
        End = 830,
        intensity.args = list(bandwidth = 30)
    )
# }
```
