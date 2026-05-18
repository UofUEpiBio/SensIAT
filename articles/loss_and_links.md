# Fitting with Different Loss and Link Functions

``` r

library(SensIAT)
```

## Introduction

In this framework let $`A=a`$ denote the treatment arm, $`Y`$ the
outcome variable, and $`t`$ the time variable. The SensIAT package
allows users to specify different loss functions and link functions when
fitting models to their data. This flexibility enables users to tailor
the modeling approach to their specific research questions and data
characteristics. Let $`\mu_a(t) = E\{Y(t)\mid A=a\}`$ denote the
population mean outcome at time $`t`$ were all participants assigned to
treatment arm $`A=a`$. The mean function is modeled as \$\_a(t) =
g^({-1}((t))\_a) = s((t)^\_a) \$ where $`g`$ is a link function,
$`\mathbf{B}`$ is a vector values basis function and $`\gamma_a`$ is a
vector of coefficients, distinct for each treatment arm. The function
$`s = g^{-1}`$ is the inverse link function. We consider three loss
functions to use in order to fit the coefficient $`\gamma_a`$:

- $`\mathcal{L}_1(\boldsymbol{\gamma}) = \int_{t=t_1}^{t_2} \left[ g\left\{\mu(t)\right\} - \boldsymbol{B}(t)^\prime \boldsymbol{\gamma}\right]^2 dt`$,
  squared error loss in the transformed space;
- $`\mathcal{L}_2(\boldsymbol{\gamma}) = \int_{t=t_1}^{t_2} \left[ \mu(t) - s\left\{\boldsymbol{B}(t)^\prime \boldsymbol{\gamma}\right\}\right]^2 dt`$,
  squared error loss in the original space;
- $`\mathcal{L}_3(\boldsymbol{\gamma}) = \int_{t=t_1}^{t_2} \left[ b\{\boldsymbol{B}(t)^\prime \boldsymbol{\gamma}\} - \mu(t)\boldsymbol{B}(t)^\prime \boldsymbol{\gamma}\right] dt`$,
  quasi-likelihood loss, where
  $`s(z) = \frac{\partial b(z)}{\partial z }`$. In the SensIAT package,
  we have implemented the following link functions:
- Identity link: $`g(\mu) = \mu`$, $`s(z) = z`$,
- Log link: $`g(\mu) = \log(\mu)`$, $`s(z) = \exp(z)`$,
- Logit link: $`g(\mu) = \log\left(\frac{\mu}{1-\mu}\right)`$,
  $`s(z) = \frac{\exp(z)}{1+\exp(z)}`$.

## Details

The estimate for the treatment group marginal mean function,
$`\widehat\mu_{\mathcal{L}, g}(t) = g^{-1}\big\{\boldsymbol{B}(t)^\prime \widehat{\boldsymbol{\beta}}\big\}`$
is found by solving
$`\frac{1}{n}\sum_{i=1}^n \widehat{\boldsymbol{\Psi}}_{\mathcal{L},g}(\boldsymbol{O}_i;\boldsymbol{\beta}) = 0`$,
for $`\beta`$, where

``` math
   \widehat{\boldsymbol{\Psi}}_{\mathcal{L},g}(\boldsymbol{O}_i;\boldsymbol{\beta})
    =
    \sum_{k=1}^{K_i}  \Bigg\{ W_{\mathcal{L}, g}(T_{ik};\boldsymbol{\beta})\frac{\big[Y_i(T_{ik}) - \widehat{\mathbb{E}} \big\{Y(T_{ik}) \mid  \overline{\boldsymbol{O}}_i(T_{ik}) \big\} \big] }{\widehat{\rho}
    \big \{T_{ik} \mid \overline{\boldsymbol{O}}_i(T_{ik}), Y_i(T_{ik}) \big\}} \Bigg\}  
    +  
    \int_{t=t_1}^{t_2} W_{\mathcal{L}, g}(t\mid\boldsymbol{\beta})
    \left[
        \widehat{\mathbb{E}} \left\{
            Y(t)\mid \overline{\boldsymbol{O}}_i(t)
        \right\} 
      -  s\left\{ 
                \boldsymbol{B}(t)^\prime \boldsymbol{\beta}
          \right\}
    \right]dt
```
The $`W_{\mathcal{L}, g}(t;\boldsymbol{\beta})`$ term depends on the
choice of loss function and link function. Then we solve
$`\frac{1}{n}\sum_{i=1}^n \widehat{\boldsymbol{\Psi}}_{\mathcal{L},g}(\boldsymbol{O}_i;\boldsymbol{\beta}) = 0`$
for $`\beta`$ using a vectorized root-finding algorithm.

### Squared Error Loss in the Transformed Space

For squared error loss in the transformed space, $`\mathcal{L}_1`$, in
general we have
``` math
    W_{\mathcal{L},g}(t;\boldsymbol{\beta}) = W_1(t;\boldsymbol{\beta}) = \boldsymbol{V}_1^{-1}\boldsymbol{B}(t)\left.\frac{\partial g(z)}{\partial z} \right|_{z = s\left\{\boldsymbol{B}(t)^\prime \boldsymbol{\beta}\right\}}
    , \quad 
    \boldsymbol{V}_1 = \int_{t=t_1}^{t_2} \boldsymbol{B}(t) \boldsymbol{B}(t)^\prime dt
```

#### Identity Link

For the identity link, $`g(\mu) = \mu`$,
$`\frac{\partial g(z)}{\partial z}\equiv 1`$, so
``` math
W_1(t;\boldsymbol{\beta}) = \boldsymbol{V}_1^{-1}\boldsymbol{B}(t).
```

#### Log Link

For the log link,
``` math
g(\mu) = \log(\mu),
```

$`\frac{\partial g(z)}{\partial z} = \frac{1}{z}`$,

``` math
 
W_1(t;\boldsymbol{\beta}) = 
    \frac{\boldsymbol{V}_1^{-1}\boldsymbol{B}(t)} {\boldsymbol{B}(t)^\prime \boldsymbol{\beta}}.
```

#### Logit Link

For the logit link,
``` math
\begin{align}
g(\mu) &= \log\left(\frac{\mu}{1-\mu}\right), \\
\frac{\partial g(z)}{\partial z} &= \frac{1}{z(1-z)}, \\
W_1(t;\boldsymbol{\beta}) &= \boldsymbol{V}_1^{-1}\boldsymbol{B}(t) \left.\frac{1}{z(1-z)}\right|_{z=s\{\boldsymbol{B}(t)^\prime\boldsymbol{\beta}\}} = \boldsymbol{V}_1^{-1}\boldsymbol{B}(t) \frac{\left[1+\exp\{\boldsymbol{B}(t)^\prime\boldsymbol{\beta}\}\right]^2}{\exp\{\boldsymbol{B}(t)^\prime\boldsymbol{\beta}\}}.
\end{align}
```

### Quasi-likelihood Loss

For quasi-likelihood loss, $`\mathcal{L}_3`$, in general we have
``` math
    W_{\mathcal{L}_3,g}(t;\boldsymbol{\beta}) = W_3(t;\boldsymbol{\beta}) = 
        \boldsymbol{V}_3(\boldsymbol{\beta})^{-1} 
        \boldsymbol{B}(t) 
    , \quad 
    \boldsymbol{V}_3(\boldsymbol{\beta}) = \int_{t=t_1}^{t_2} 
        \boldsymbol{B}(t) 
        \boldsymbol{B}(t) ^ \prime
        \left.
            \frac{\partial s(z)}{\partial z} 
        \right|_{z=\boldsymbol{B}(t)^\prime \boldsymbol{\beta}}
        dt
```

#### Identity Link

For the identity link, $`s(z) = z`$, so
``` math
\frac{\partial s(z)}{\partial z} \equiv 1,
```
and
``` math
W_3(t;\boldsymbol{\beta}) = \boldsymbol{V}_3^{-1}\boldsymbol{B}(t),
```
where
``` math
\boldsymbol{V}_3 = \int_{t=t_1}^{t_2} \boldsymbol{B}(t) \boldsymbol{B}(t)^\prime dt.
```

#### Log Link

For the log link, $`s(z) = \exp(z)`$, so
``` math
\frac{\partial s(z)}{\partial z} = \exp(z),
```
and
``` math
W_3(t;\boldsymbol{\beta}) = \boldsymbol{V}_3(\boldsymbol{\beta})^{-1} \boldsymbol{B}(t),
```
where
``` math
\boldsymbol{V}_3(\boldsymbol{\beta}) = \int_{t=t_1}^{t_2} \boldsymbol{B}(t) \boldsymbol{B}(t)^\prime \exp\left(\boldsymbol{B}(t)^\prime \boldsymbol{\beta}\right) dt.
```

#### Logit Link

For the logit link, $`s(z) = \frac{\exp(z)}{1+\exp(z)}`$, so
``` math
\frac{\partial s(z)}{\partial z} = \frac{\exp(z)}{\left\{1+\exp(z)\right\}^2},
```
and
``` math
W_3(t;\boldsymbol{\beta}) = \boldsymbol{V}_3(\boldsymbol{\beta})^{-1} \boldsymbol{B}(t),
```
where
``` math
\boldsymbol{V}_3(\boldsymbol{\beta}) = \int_{t=t_1}^{t_2} \boldsymbol{B}(t) \boldsymbol{B}(t)^\prime \frac{\exp\left(\boldsymbol{B}(t)^\prime \boldsymbol{\beta}\right)}{\left\{1+\exp\left(\boldsymbol{B}(t)^\prime \boldsymbol{\beta}\right)\right\}^2} dt.
```
