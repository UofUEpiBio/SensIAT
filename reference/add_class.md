# Adds an S3 class to an object

Adds an S3 class to an object

## Usage

``` r
add_class(x, class)
```

## Arguments

- x:

  An object to which the class should be added.

- class:

  A character vector of class names to be added.

## Examples

``` r
add_class(TRUE, 'flag')
#> [1] TRUE
#> attr(,"class")
#> [1] "flag"    "logical"
```
