# Single Index Model using `MAVE` and Optimizing Bandwidth.

Single index model estimation using minimum average variance estimation
(`MAVE`). A direction is estimated using `MAVE`, and then the bandwidth
is selected by minimization of the cross-validated pseudo-integrated
squared error. Optionally, the initial coefficients of the outcome model
can be re-estimated by optimization on a spherical manifold. This option
requires the ManifoldOptim package (see
[`ManifoldOptim::manifold.optim()`](https://rdrr.io/pkg/ManifoldOptim/man/manifold.optim.html)).

## Usage

``` r
fit_SensIAT_single_index_norm1coef_model(
  formula,
  data,
  kernel = "K2_Biweight",
  mave.method = "meanMAVE",
  id = ..id..,
  bw.selection = c("ise", "mse"),
  bw.method = c("optim", "grid", "optimize"),
  bw.range = c(0.01, 1.5),
  reestimate.coef = 0,
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

- mave.method:

  The method to use for the `MAVE` estimation.

- id:

  The patient identifier variable for the data.

- bw.selection:

  The criteria for bandwidth selection, either `'ise'` for Integrated
  Squared Error or `'mse'` for Mean Squared Error.

- bw.method:

  The method for bandwidth selection, either `'optim'` for using
  optimization or `'grid'` for grid search.

- bw.range:

  A numeric vector of length 2 indicating the range of bandwidths to
  consider for selection as a multiple of the standard deviation of the
  single index predictor.

- reestimate.coef:

  number of iterations to go through.

- ...:

  Additional arguments to be passed to
  [optim](https://rdrr.io/r/stats/optim.html).

## Value

Object of class `SensIAT::Single-index-outcome-model` which contains the
outcome model portion.
