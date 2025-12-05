# Compute Conditional Expected Values based on Outcome Model

Compute Conditional Expected Values based on Outcome Model

## Usage

``` r
compute_SensIAT_expected_values(
  model,
  alpha = 0,
  new.data = model.frame(model),
  ...
)

# S3 method for class 'lm'
compute_SensIAT_expected_values(model, alpha, new.data, ...)

# S3 method for class 'glm'
compute_SensIAT_expected_values(
  model,
  alpha,
  new.data,
  ...,
  y.max = NULL,
  eps = .Machine$double.eps
)

# S3 method for class 'negbin'
compute_SensIAT_expected_values(
  model,
  alpha,
  new.data,
  ...,
  y.max = NULL,
  eps = .Machine$double.eps^(1/4)
)
```

## Arguments

- model:

  An object representing the output of the outcome model.

- alpha:

  The sensitivity parameter

- new.data:

  Data to compute conditional means for, defaults to the model frame for
  the fitted model.

- ...:

  passed onto methods.

- y.max:

  The maximum value of the outcome variable for the Poisson and Negative
  Binomial models. If omitted it is chosen from the quantile function
  for the distribution at `1-eps`.

- eps:

  The tolerance for the quantile function used to estimate `y.max`,
  default is `.Machine$double.eps`.

## Value

The `new.data` frame with additional columns `alpha`, `E_Yexp_alphaY`,
and `E_exp_alphaY` appended.

## Details

Compute the conditional expectations needed for predictions in the
models. Two additional values/expectations are computed:

- `$E \big[ Y(t) \exp \{ \alpha Y(t) \} | A(t)=1, \bar{O}(t) \big]$`,
  returned as `E_Yexp_alphaY`, and

- `$E \big[ \exp \{ \alpha Y(t) \} \ | A(t)=1, \bar{O}(t) \big]$`,
  returned as `E_exp_alphaY`.

For the methods shown here

## Methods (by class)

- `compute_SensIAT_expected_values(lm)`: (Gaussian) Linear Model method
  The [stats::integrate](https://rdrr.io/r/stats/integrate.html) method
  is used to compute the conditional expectations.

- `compute_SensIAT_expected_values(glm)`: Generalized Linear Model
  method

- `compute_SensIAT_expected_values(negbin)`: Negative Binomial Model
  method

## Examples

``` r
model <- lm(mpg ~ as.factor(cyl)+disp+wt, data=mtcars)
compute_SensIAT_expected_values(model, alpha= c(-0.3, 0, 0.3), new.data = mtcars[1:5, ])
#>                         mpg cyl disp  hp drat    wt  qsec vs am gear carb
#> Mazda RX4...1          21.0   6  160 110 3.90 2.620 16.46  0  1    4    4
#> Mazda RX4 Wag...2      21.0   6  160 110 3.90 2.875 17.02  0  1    4    4
#> Datsun 710...3         22.8   4  108  93 3.85 2.320 18.61  1  1    4    1
#> Hornet 4 Drive...4     21.4   6  258 110 3.08 3.215 19.44  1  0    3    1
#> Hornet Sportabout...5  18.7   8  360 175 3.15 3.440 17.02  0  0    3    2
#> Mazda RX4...6          21.0   6  160 110 3.90 2.620 16.46  0  1    4    4
#> Mazda RX4 Wag...7      21.0   6  160 110 3.90 2.875 17.02  0  1    4    4
#> Datsun 710...8         22.8   4  108  93 3.85 2.320 18.61  1  1    4    1
#> Hornet 4 Drive...9     21.4   6  258 110 3.08 3.215 19.44  1  0    3    1
#> Hornet Sportabout...10 18.7   8  360 175 3.15 3.440 17.02  0  0    3    2
#> Mazda RX4...11         21.0   6  160 110 3.90 2.620 16.46  0  1    4    4
#> Mazda RX4 Wag...12     21.0   6  160 110 3.90 2.875 17.02  0  1    4    4
#> Datsun 710...13        22.8   4  108  93 3.85 2.320 18.61  1  1    4    1
#> Hornet 4 Drive...14    21.4   6  258 110 3.08 3.215 19.44  1  0    3    1
#> Hornet Sportabout...15 18.7   8  360 175 3.15 3.440 17.02  0  0    3    2
#>                        E_Yexp_alphaY E_exp_alphaY
#> Mazda RX4...1           4.335670e-02 2.244828e-03
#> Mazda RX4 Wag...2       5.339876e-02 2.890980e-03
#> Datsun 710...3          1.153856e-02 4.705305e-04
#> Hornet 4 Drive...4      6.745769e-02 3.851517e-03
#> Hornet Sportabout...5   1.249134e-01 8.367589e-03
#> Mazda RX4...6           2.134680e+01 1.000000e+00
#> Mazda RX4 Wag...7       2.050358e+01 1.000000e+00
#> Datsun 710...8          2.655522e+01 1.000000e+00
#> Hornet 4 Drive...9      1.954734e+01 1.000000e+00
#> Hornet Sportabout...10  1.696102e+01 1.000000e+00
#> Mazda RX4...11          1.916457e+04 8.197146e+02
#> Mazda RX4 Wag...12      1.434446e+04 6.365034e+02
#> Datsun 710...13         1.117999e+05 3.910732e+03
#> Hornet 4 Drive...14     1.031021e+04 4.777646e+02
#> Hornet Sportabout...15  4.176927e+03 2.199103e+02
model <- glm(cyl ~ mpg+disp+wt, data=mtcars, family=poisson())
compute_SensIAT_expected_values(model, alpha= c(-0.3, 0, 0.3), new.data = mtcars[1:5, ]) |>
    dplyr::mutate('E(y|alpha)' = .data$E_Yexp_alphaY/.data$E_exp_alphaY)
#> # A tibble: 15 × 15
#>      mpg   cyl  disp    hp  drat    wt  qsec    vs    am  gear  carb alpha
#>    <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#>  1  21       6   160   110  3.9   2.62  16.5     0     1     4     4  -0.3
#>  2  21       6   160   110  3.9   2.88  17.0     0     1     4     4  -0.3
#>  3  22.8     4   108    93  3.85  2.32  18.6     1     1     4     1  -0.3
#>  4  21.4     6   258   110  3.08  3.22  19.4     1     0     3     1  -0.3
#>  5  18.7     8   360   175  3.15  3.44  17.0     0     0     3     2  -0.3
#>  6  21       6   160   110  3.9   2.62  16.5     0     1     4     4   0  
#>  7  21       6   160   110  3.9   2.88  17.0     0     1     4     4   0  
#>  8  22.8     4   108    93  3.85  2.32  18.6     1     1     4     1   0  
#>  9  21.4     6   258   110  3.08  3.22  19.4     1     0     3     1   0  
#> 10  18.7     8   360   175  3.15  3.44  17.0     0     0     3     2   0  
#> 11  21       6   160   110  3.9   2.62  16.5     0     1     4     4   0.3
#> 12  21       6   160   110  3.9   2.88  17.0     0     1     4     4   0.3
#> 13  22.8     4   108    93  3.85  2.32  18.6     1     1     4     1   0.3
#> 14  21.4     6   258   110  3.08  3.22  19.4     1     0     3     1   0.3
#> 15  18.7     8   360   175  3.15  3.44  17.0     0     0     3     2   0.3
#> # ℹ 3 more variables: E_Yexp_alphaY <dbl>, E_exp_alphaY <dbl>,
#> #   `E(y|alpha)` <dbl>
```
