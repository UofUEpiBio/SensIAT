# Fit Single Index Dimension Reduction with Fixed Bandwidth

Fit Single Index Dimension Reduction with Fixed Bandwidth

## Usage

``` r
SIDRnew_fixed_bandwidth(
  X,
  Y,
  ids,
  Y.CP = NULL,
  initial = NULL,
  kernel = "K2_Biweight",
  method = c("optim", "nlminb", "nmk"),
  optim_method = "BFGS",
  abs.tol = 1e-04
)
```
