# Generate Next Observation Time

Generate the time of the next observation using a Cox proportional
hazards model. The hazard depends on the previous outcome value and is
stratified by visit number.

## Usage

``` r
generate_next_observation_time(
  current_time,
  current_outcome,
  visit_num,
  intensity_coef,
  baseline_hazard,
  End
)
```

## Arguments

- current_time:

  Current time point.

- current_outcome:

  Current outcome value.

- visit_num:

  Current visit number (for stratification).

- intensity_coef:

  Coefficient for the effect of previous outcome on observation
  intensity. Can be a scalar or vector (one per visit number stratum).

- baseline_hazard:

  Baseline hazard function. Either a function of time and visit number,
  or a numeric value for constant baseline hazard.

- End:

  Maximum follow-up time.

## Value

Time of next observation.
