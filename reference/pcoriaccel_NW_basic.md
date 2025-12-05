# Runs a *basic* implementation of the "NW" function with the "K2_Biweight" kernel, just as a proof-of-concept.

Runs a *basic* implementation of the "NW" function with the
"K2_Biweight" kernel, just as a proof-of-concept.

## Usage

``` r
pcoriaccel_NW_basic(Xb, Y, xb, y_seq, h)
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

## Value

A matrix of the same size as `xb` by `y_seq`.
