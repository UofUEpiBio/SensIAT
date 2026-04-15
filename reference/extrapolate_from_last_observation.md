# Extrapolate Variables from Last Observation

Extrapolates variables from the last observation before a given time
point using linear extrapolation with specified slopes.

## Usage

``` r
extrapolate_from_last_observation(
  target_time,
  data,
  time_var,
  slopes = NULL,
  strict = TRUE
)
```

## Arguments

- target_time:

  Numeric scalar. The time point at which to extrapolate values.

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

A single-row data.frame with extrapolated values at `target_time`.

## Details

The function finds the last observation in `data` where the time is less
than or equal to `target_time`, then extrapolates variables using the
formula:

\$\$x\_{extrapolated} = x\_{last} + (target\\time - time\_{last}) \times
slope\$\$

For variables not specified in `slopes`, the last observed value is used
(slope = 0).

## Examples

``` r
# \donttest{
# Create example data
df <- data.frame(
    time = c(0, 1, 2, 3),
    x1 = c(10, 12, 14, 16),
    x2 = c(5, 5.5, 6, 6.5),
    id = 1
)

# Extrapolate to time = 5 with known slopes
extrapolate_from_last_observation(
    target_time = 5,
    data = df,
    time_var = "time",
    slopes = c(x1 = 2, x2 = 0.5)
)
#>   time x1  x2 id
#> 4    5 20 7.5  1
# Expected: x1 = 16 + (5-3)*2 = 20, x2 = 6.5 + (5-3)*0.5 = 7.5

# Extrapolate with only some variables having slopes
extrapolate_from_last_observation(
    target_time = 4,
    data = df,
    time_var = "time",
    slopes = c(x1 = 2)
)
#>   time x1  x2 id
#> 4    4 18 6.5  1
# Expected: x1 = 16 + (4-3)*2 = 18, x2 = 6.5 (carried forward), id = 1
# }
```
