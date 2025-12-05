# Directly estimate the probability mass function of Y.

Directly estimate the probability mass function of Y.

## Usage

``` r
pcoriaccel_estimate_pmf(Xb, Y, xi, y_seq, h, kernel = "K2_Biweight")
```

## Arguments

- Xb:

  Numeric vector of individual linear predictors from the data

- Y:

  Numeric vector of individual responses from the data

- xi:

  value of the individuals linear predictor at the point of estimation

- y_seq:

  Numeric vector of unique values of `Y`.

- h:

  bandwidth of the kernel

- kernel:

  character string specifying the kernel to use, either `"dnorm"`,
  `"K2_Biweight"`, or `"K4_Biweight"`
