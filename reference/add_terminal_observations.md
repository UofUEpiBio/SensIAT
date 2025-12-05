# Add Terminal Observations to a Dataset

This function adds terminal observations to a dataset. For each subject
given by `id`, if that subject has less than the maximum number of
observations, A row is added with the `end` time value, leaving all
other variables as `NA`.

## Usage

``` r
add_terminal_observations(
  data,
  id,
  time,
  end = max(pull(data, {
     {
         time
     }
 }))
)
```

## Arguments

- data:

  A data frame containing the dataset.

- id:

  A variable in `data` that identifies the subject.

- time:

  A variable in `data` that identifies the time of the observation.

- end:

  The value to use for the `time` variable in the terminal observation.
  If end is less that the maximum in the dataset resulting data will be
  filtered such that `time` is less than or equal to `end`.

## Value

A data frame with terminal observations added.

## See also

[`tidyr::complete()`](https://tidyr.tidyverse.org/reference/complete.html)

## Examples

``` r
exdata <- tibble::tibble(
  patient = rep(1:3, 3:5),
  day = c(0, 30, 60,
          0, 30, 60, 90,
          0, 30, 60, 90, 120),
  value = TRUE
)
add_terminal_observations(exdata, patient, day)
#> # A tibble: 14 × 3
#> # Groups:   patient [3]
#>    patient   day value
#>      <int> <dbl> <lgl>
#>  1       1     0 TRUE 
#>  2       1    30 TRUE 
#>  3       1    60 TRUE 
#>  4       1   120 NA   
#>  5       2     0 TRUE 
#>  6       2    30 TRUE 
#>  7       2    60 TRUE 
#>  8       2    90 TRUE 
#>  9       2   120 NA   
#> 10       3     0 TRUE 
#> 11       3    30 TRUE 
#> 12       3    60 TRUE 
#> 13       3    90 TRUE 
#> 14       3   120 TRUE 
```
