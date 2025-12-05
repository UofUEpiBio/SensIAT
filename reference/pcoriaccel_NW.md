# Runs an optimized implementation of the "NW" function.

Runs an optimized implementation of the "NW" function.

## Usage

``` r
pcoriaccel_NW(Xb, Y, xb, y_seq, h, kernel = "K2_Biweight")
```

## Arguments

- Xb:

  a vector (expected to be about 500 elements)

- Y:

  a vector (same size as `Xb`)

- xb:

  a vector

- y_seq:

  a vector

- h:

  a scalar, the bandwidth of kernel

- kernel:

  a string, denoting the kernel function to use, either `"dnorm"`,
  `"K2_Biweight"`, or `"K4_Biweight"`

## Value

A matrix of the same size as `xb` by `y_seq`.
