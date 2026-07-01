# Produce fitted model for group (treatment or control)

Produces a fitted model that may be used to produce estimates of mean
and variance for the given group.

## Usage

``` r
fit_SensIAT_fulldata_model(data, trt, ...)

fit_SensIAT_within_group_model(
  group.data,
  outcome_modeler,
  id,
  outcome,
  time,
  knots,
  alpha = 0,
  End = NULL,
  intensity.args = list(),
  outcome.args = list(),
  influence.args = list(),
  spline.degree = 3,
  add.terminal.observations = TRUE,
  loss = c("lp_mse", "quasi-likelihood"),
  link = c("identity", "log", "logit"),
  term2_method = c("fast", "original", "fixed_grid", "seeded_adaptive", "gauss_legendre"),
  impute_data = NULL,
  progress = FALSE
)
```

## Arguments

- data:

  the full data set.

- trt:

  an expression that determine what is treated as the treatment.
  Everything not treatment is considered control.

- ...:

  common arguments passed to `fit_SensIAT_within_group_model`.

- group.data:

  The data for the group that is being analyzed. Preferably passed in as
  a single `tibble` that internally is subsetted/filtered as needed.

- outcome_modeler:

  function for fitting the outcome model. Called with a formula, data
  argument and `outcome.args` list.

- id:

  The variable that identifies the patient.

- outcome:

  The variable that contains the outcome.

- time:

  The variable that contains the time.

- knots:

  knot locations for defining the spline basis.

- alpha:

  The sensitivity parameter.

- End:

  The end time for this data analysis, we need to set the default value
  as the max value of the time.

- intensity.args:

  A list of optional arguments for intensity model. See the Intensity
  Arguments section.

- outcome.args:

  parameters as needed passed into the `outcome_modeler`. One special
  element may be `'model.modifications'` which, if present, should be a
  formula that will be used to modify the outcome model per,
  [update.formula](https://rdrr.io/r/stats/update.formula.html).

- influence.args:

  A list of optional arguments used when computing the influence. See
  the Influence Arguments section.

- spline.degree:

  The degree of the spline basis.

- add.terminal.observations:

  Logical indicating whether to add terminal observations to the data.
  If TRUE, data may not contain any `NA`s. if FALSE, data will be
  assumed to already include the terminal observations

- loss:

  The loss function to use. Options are `"lp_mse"` (default) and
  `"quasi-likelihood"`. Only used when `link` is not `"identity"`.

- link:

  The link function to use. Options are `"identity"` (default), `"log"`,
  and `"logit"`. When `"identity"`, the original analytic formula is
  used; otherwise, the generalized iterative solver is used.

- term2_method:

  Method for computing term2 influence components. Options are: `"fast"`
  (default), `"original"`, `"fixed_grid"`, `"seeded_adaptive"`,
  `"gauss_legendre"`. Only used when `link` is not `"identity"`.

- impute_data:

  A function that takes `(t, df)` and returns the imputed data at time
  `t`. If `NULL` (default), a standard imputation function is generated
  automatically.

- progress:

  Logical indicating whether to print stage progress messages. Default
  is `FALSE`.

## Value

a list with class `SensIAT-fulldata-fitted-model` with two components,
`control` and `treatment`, each of which is an independently fitted
`SensIAT-within-group-fitted-model` fit with the fit_within_group_model
function.

A `SensIAT_within_group_model` object containing:

- `models`: List with `intensity` and `outcome` sub-models

- `data`: The transformed data used for fitting

- `variables`: List of variable names (id, outcome, time)

- `End`: The end time for the analysis

- `influence`: Influence function values by patient

- `alpha`: Sensitivity parameter values

- `coefficients`: Estimated spline coefficients for each alpha

- `coefficient.variance`: Variance of coefficients for each alpha

- `base`: The
  [`SplineBasis`](https://rdrr.io/pkg/orthogonalsplinebasis/man/SplineBasis.html)
  object

- `V_inverse`: Inverse of the Gram matrix

- `link`: The link function used

- `loss`: The loss function used

- `term2_method`: The term2 integration method used

## Details

This function should be agnostic to whether it is being provided a
treatment or control group.

## Functions

- `fit_SensIAT_fulldata_model()`: Fit the Marginal Mean Sensitivity
  Analysis with Tilting Assumption for Both Treatment and Control
  Groups.

## Intensity Arguments

The `intensity.args` list may contain the following elements:

- **`model.modifications`** A formula that will be used to modify the
  intensity model from it's default, per
  [update.formula](https://rdrr.io/r/stats/update.formula.html).

- **`kernel`** The kernel function for the intensity model. Default is
  the Epanechnikov kernel.

- **`bandwidth`** The bandwidth for the intensity model kernel.

## Influence Arguments

The `influence.args` list may contain the following elements:

- **`method`** The method for integrating, adaptive or fixed quadrature.
  Default is `'adaptive'`.

- **`tolerance`** The tolerance when using adaptive quadrature.

- **`delta`** The bin width for fixed quadrature.

- **`resolution`** alternative to `delta` by specifying the number of
  bins.

- **`fix_discontinuity`** Whether to account for the discontinuity in
  the influence at observation times.

## Examples

``` r
# \donttest{
# Example 1: Default identity link (original analytic solution)
model <-
    fit_SensIAT_within_group_model(
        group.data = SensIAT_example_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        End = 830,
        knots = c(60, 260, 460),
    )

# Example 2: Log link with generalized iterative solver
model_log <-
    fit_SensIAT_within_group_model(
        group.data = SensIAT_example_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        alpha = 0,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        End = 830,
        knots = c(60, 260, 460),
        link = "log",
        loss = "lp_mse"
    )
# }
```
