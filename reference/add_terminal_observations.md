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
if (FALSE) { # \dontrun{
exdata <- tibble::tibble(
    patient = rep(1:3, 3:5),
    day = c(
        0, 30, 60,
        0, 30, 60, 90,
        0, 30, 60, 90, 120
    ),
    value = TRUE
)
add_terminal_observations(exdata, patient, day)
} # }
```
