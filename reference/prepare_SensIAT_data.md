# Prepare Data for Sensitivity Analysis with Irregular Assessment Times

This function prepares the data for SensIAT analysis by transforming it
into a format suitable for the SensIAT models.

## Usage

``` r
prepare_SensIAT_data(
  data,
  id.var,
  time.var,
  outcome.var,
  End,
  add.terminal.observations = TRUE
)
```

## Arguments

- data:

  A data frame containing the data to be prepared.

- id.var:

  The variable in `data` that identifies the subject.

- time.var:

  The variable in `data` that identifies the time of the observation.

- outcome.var:

  The variable in `data` that contains the outcome of interest.

- End:

  The end time for the analysis. Observations with time greater than
  `End` will be filtered out.

- add.terminal.observations:

  Logical indicating whether to add terminal observations to the data
  (`TRUE`), or terminal observations have already been added (`FALSE`).

## Value

A data frame with the following transformations:

- Data filtered to time less than or equal to `End`.

- Observations are arranged by `id.var` and `time.var`.

- Terminal observations added if `add.terminal.observations` is `TRUE`,
  with `..time..` set to `End` and `..outcome..` set to `NA`, if the
  subject has less observations than the maximum number of observations.

- New variables created:

  - `..id..` aliases `id.var`,

  - `..time..` aliases `time.var`,

  - `..outcome..` aliases `outcome.var`,

  - `..visit_number..` is the visit number within each subject derived
    from `time.var`,

  - `..prev_outcome..`, i.e. lag-outcome, the outcome from the previous
    visit,

  - `..prev_time..`, i.e. lag-time, the time from the previous visit,

  - `..delta_time..`, the difference in time between the current and
    previous visit.

## Examples

``` r
prepare_SensIAT_data( SensIAT_example_data, Subject_ID, Time, Outcome, 830)
#> # A tibble: 1,000 × 11
#>    ..id.. ..visit_number.. Subject_ID Visit  Time Outcome ..time.. ..outcome..
#>     <int>            <int>      <int> <dbl> <dbl>   <dbl>    <dbl>       <dbl>
#>  1      1                0          1     0     0   3            0       3    
#>  2      1                1          1     1   214   4.5        214       4.5  
#>  3      1                2          1     2   292   4.17       292       4.17 
#>  4      1                3          1     3   370   1.33       370       1.33 
#>  5      1                4          1     4   441   0.833      441       0.833
#>  6      2                0          2     0     0   3            0       3    
#>  7      2                1          2     1    72   0.5         72       0.5  
#>  8      2                2          2     2   181   2          181       2    
#>  9      2                3          2     3   297   1.5        297       1.5  
#> 10      2                4          2     4   366   1.83       366       1.83 
#> # ℹ 990 more rows
#> # ℹ 3 more variables: ..prev_outcome.. <dbl>, ..prev_time.. <dbl>,
#> #   ..delta_time.. <dbl>

exdata <- tibble::tibble(ID=rep(1:2, c(3,5)),
                         Time=c(0, 30, 60,
                                0, 30, 60, 90, 120),
                         Outcome=floor(runif(8, 1, 100)))

prepare_SensIAT_data(exdata, ID, Time, Outcome, 120)
#> # A tibble: 10 × 10
#>    ..id.. ..visit_number..    ID  Time Outcome ..time.. ..outcome..
#>     <int>            <int> <int> <dbl>   <dbl>    <dbl>       <dbl>
#>  1      1                0     1     0       6        0           6
#>  2      1                1     1    30      53       30          53
#>  3      1                2     1    60      69       60          69
#>  4      1                3     1   120      NA      120          NA
#>  5      1                4     1   120      NA      120          NA
#>  6      2                0     2     0      69        0          69
#>  7      2                1     2    30       4       30           4
#>  8      2                2     2    60      23       60          23
#>  9      2                3     2    90      30       90          30
#> 10      2                4     2   120      64      120          64
#> # ℹ 3 more variables: ..prev_outcome.. <dbl>, ..prev_time.. <dbl>,
#> #   ..delta_time.. <dbl>
```
