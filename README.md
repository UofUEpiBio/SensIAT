
[![R-CMD-check](https://github.com/UofUEpi/pcoriRpackage/actions/workflows/r.yml/badge.svg)](https://github.com/UofUEpi/pcoriRpackage/actions/workflows/r.yml)

# SensIAT

R Package for sensitivity analysis with irregular assessment times.

## Setup

Installation of required dependencies:

``` r
install.packages(c( "assertthat", "dplyr", "orthogonalsplinebasis", "pracma", "purrr", "Rcpp", "rlang", "roxygen2", "tibble", "tidyr" ))
```

On Windows, you will also need
[`RTools`](https://cran.r-project.org/bin/windows/Rtools/).

To build and install, from the project root directory, do:

``` sh
R CMD INSTALL .
```
