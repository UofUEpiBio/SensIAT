# Extrapolate Multiple Time Points from Last Observation

Vectorized version of
[`extrapolate_from_last_observation()`](https://uofuepibio.github.io/SensIAT/reference/extrapolate_from_last_observation.md)
that extrapolates to multiple time points at once.

## Usage

``` r
extrapolate_from_last_observation_multiple(
  target_times,
  data,
  time_var,
  slopes = NULL,
  strict = TRUE
)
```

## Arguments

- target_times:

  Numeric vector. Multiple time points at which to extrapolate values.

- data:

  A data.frame containing the observations, with one row per time point.
  Must be sorted by the time variable.

- time_var:

  Character string. The name of the time variable in `data`.

- slopes:

  A named numeric vector. Names should match column names in `data` for
  variables to extrapolate, and values are the slopes (rate of change
  per unit time). Variables not in `slopes` will be carried forward
  unchanged from the last observation.

- strict:

  Logical. If `TRUE` (default), requires that `target_time` is greater
  than or equal to the last observed time. If `FALSE`, allows
  extrapolation to earlier times.

## Value

A data.frame with one row per target time, containing extrapolated
values.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- data.frame(
    time = c(0, 1, 2, 3),
    x1 = c(10, 12, 14, 16),
    x2 = c(5, 5.5, 6, 6.5)
)

extrapolate_from_last_observation_multiple(
    target_times = c(3.5, 4, 5),
    data = df,
    time_var = "time",
    slopes = c(x1 = 2, x2 = 0.5)
)
} # }
```
