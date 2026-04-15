# Adds an S3 class to an object

This is a utility function for working with the S3 class system only. It
does not work with S4 or R6 class systems.

## Usage

``` r
add_class(x, class)
```

## Arguments

- x:

  An object to which the class should be added.

- class:

  A character vector of class names to be added.

## Value

The object with the new class(es) prepended to existing classes.

## Examples

``` r
if (FALSE) { # \dontrun{
# Internal use only
obj <- add_class(list(a = 1), "my_class")
class(obj)  # c("my_class", "list")
} # }
```
