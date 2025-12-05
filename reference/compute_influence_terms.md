# Compute Influence Terms

This function computes the influence terms for the marginal outcome
model sensitivity analysis. It is a generic function that can handle
different types of outcome models.

## Usage

``` r
compute_influence_terms(outcome.model, intensity.model, alpha, data, ...)

# Default S3 method
compute_influence_terms(
  outcome.model,
  intensity.model,
  alpha,
  data,
  id,
  base,
  ...
)

# S3 method for class '`SensIAT::Single-index-outcome-model`'
compute_influence_terms(
  outcome.model,
  intensity.model,
  alpha,
  data,
  base,
  tolerance = .Machine$double.eps^(1/3),
  na.action = na.fail,
  id = NULL,
  time = NULL,
  ...
)
```

## Arguments

- outcome.model:

  The outcome model fitted to the data.

- intensity.model:

  The intensity model fitted to the data.

- alpha:

  A numeric vector representing the sensitivity parameter.

- data:

  A data frame containing the observations.

- ...:

  Additional arguments passed to the method.

- id:

  A variable representing the patient identifier.

- base:

  A spline basis object.

- tolerance:

  Numeric value indicating the tolerance for integration, default is
  `.Machine$double.eps^(1/3)`.

- na.action:

  Function to handle missing values, default is `na.fail`.

- time:

  Variable indicating the time variable in the data, by Default will be
  extracted from the intensity model response.

## Methods (by class)

- `compute_influence_terms(default)`: Generic method, which throws a not
  implemented error.

- `` compute_influence_terms(`SensIAT::Single-index-outcome-model`) ``:
  Optimized method for the single index model.
