# Fitting with Different Loss and Link Functions

``` r
library(SensIAT)
```

## Introduction

In this framework let $A = a$ denote the treatment arm, $Y$ the outcome
variable, and $t$ the time variable. The SensIAT package allows users to
specify different loss functions and link functions when fitting models
to their data. This flexibility enables users to tailor the modeling
approach to their specific research questions and data characteristics.
Let $\mu_{a}(t) = E\{ Y(t) \mid A = a\}$ denote the population mean
outcome at time $t$ were all participants assigned to treatment arm
$A = a$. The mean function is modeled as \$\_a(t) = g^({-1}((t))\_a) =
s((t)^\_a) \$ where $g$ is a link function, $\mathbf{B}$ is a vector
values basis function and $\gamma_{a}$ is a vector of coefficients,
distinct for each treatment arm. The function $s = g^{- 1}$ is the
inverse link function. We consider three loss functions to use in order
to fit the coefficient $\gamma_{a}$:

- $\mathcal{L}_{1}({\mathbf{γ}}) = \int_{t = t_{1}}^{t_{2}}\left\lbrack g\left\{ \mu(t) \right\} - \mathbf{B}(t)^{\prime}{\mathbf{γ}} \right\rbrack^{2}dt$,
  squared error loss in the transformed space;
- $\mathcal{L}_{2}({\mathbf{γ}}) = \int_{t = t_{1}}^{t_{2}}\left\lbrack \mu(t) - s\left\{ \mathbf{B}(t)^{\prime}{\mathbf{γ}} \right\} \right\rbrack^{2}dt$,
  squared error loss in the original space;
- $\mathcal{L}_{3}({\mathbf{γ}}) = \int_{t = t_{1}}^{t_{2}}\left\lbrack b\{\mathbf{B}(t)^{\prime}{\mathbf{γ}}\} - \mu(t)\mathbf{B}(t)^{\prime}{\mathbf{γ}} \right\rbrack dt$,
  quasi-likelihood loss, where
  $s(z) = \frac{\partial b(z)}{\partial z}$. In the SensIAT package, we
  have implemented the following link functions:
- Identity link: $g(\mu) = \mu$, $s(z) = z$,
- Log link: $g(\mu) = \log(\mu)$, $s(z) = \exp(z)$,
- Logit link: $g(\mu) = \log\left( \frac{\mu}{1 - \mu} \right)$,
  $s(z) = \frac{\exp(z)}{1 + \exp(z)}$.

## Details

The estimate for the treatment group marginal mean function,
${\widehat{\mu}}_{\mathcal{L},g}(t) = g^{- 1}\{\mathbf{B}(t)^{\prime}\widehat{\mathbf{β}}\}$
is found by solving
$\frac{1}{n}\sum_{i = 1}^{n}{\widehat{\mathbf{\Psi}}}_{\mathcal{L},g}\left( \mathbf{O}_{i};{\mathbf{β}} \right) = 0$,
for $\beta$, where

$${\widehat{\mathbf{\Psi}}}_{\mathcal{L},g}\left( \mathbf{O}_{i};{\mathbf{β}} \right) = \sum\limits_{k = 1}^{K_{i}}\{ W_{\mathcal{L},g}\left( T_{ik};{\mathbf{β}} \right)\frac{\lbrack Y_{i}\left( T_{ik} \right) - \widehat{\mathbb{E}}\{ Y\left( T_{ik} \right) \mid {\overline{\mathbf{O}}}_{i}\left( T_{ik} \right)\}\rbrack}{\widehat{\rho}\{ T_{ik} \mid {\overline{\mathbf{O}}}_{i}\left( T_{ik} \right),Y_{i}\left( T_{ik} \right)\}}\} + \int_{t = t_{1}}^{t_{2}}W_{\mathcal{L},g}(t \mid {\mathbf{β}})\left\lbrack \widehat{\mathbb{E}}\left\{ Y(t) \mid {\overline{\mathbf{O}}}_{i}(t) \right\} - s\left\{ \mathbf{B}(t)^{\prime}{\mathbf{β}} \right\} \right\rbrack dt$$
The $W_{\mathcal{L},g}(t;{\mathbf{β}})$ term depends on the choice of
loss function and link function. Then we solve
$\frac{1}{n}\sum_{i = 1}^{n}{\widehat{\mathbf{\Psi}}}_{\mathcal{L},g}\left( \mathbf{O}_{i};{\mathbf{β}} \right) = 0$
for $\beta$ using a vectorized root-finding algorithm.

### Squared Error Loss in the Transformed Space

For squared error loss in the transformed space, $\mathcal{L}_{1}$, in
general we have
$$W_{\mathcal{L},g}(t;{\mathbf{β}}) = W_{1}(t;{\mathbf{β}}) = \mathbf{V}_{1}^{- 1}\mathbf{B}(t)\left. \frac{\partial g(z)}{\partial z} \right|_{z = s{\{\mathbf{B}{(t)}^{\prime}{\mathbf{β}}\}}},\quad\mathbf{V}_{1} = \int_{t = t_{1}}^{t_{2}}\mathbf{B}(t)\mathbf{B}(t)^{\prime}dt$$

#### Identity Link

For the identity link, $g(\mu) = \mu$,
$\frac{\partial g(z)}{\partial z} \equiv 1$, so
$$W_{1}(t;{\mathbf{β}}) = \mathbf{V}_{1}^{- 1}\mathbf{B}(t).$$

#### Log Link

For the log link, $$g(\mu) = \log(\mu),$$

$\frac{\partial g(z)}{\partial z} = \frac{1}{z}$,

$$W_{1}(t;{\mathbf{β}}) = \frac{\mathbf{V}_{1}^{- 1}\mathbf{B}(t)}{\mathbf{B}(t)^{\prime}{\mathbf{β}}}.$$

#### Logit Link

For the logit link, $$\begin{aligned}
{g(\mu)} & {= \log\left( \frac{\mu}{1 - \mu} \right),} \\
\frac{\partial g(z)}{\partial z} & {= \frac{1}{z(1 - z)},} \\
{W_{1}(t;{\mathbf{β}})} & {= \mathbf{V}_{1}^{- 1}\mathbf{B}(t)\left. \frac{1}{z(1 - z)} \right|_{z = s\{\mathbf{B}{(t)}^{\prime}{\mathbf{β}}\}} = \mathbf{V}_{1}^{- 1}\mathbf{B}(t)\frac{\left\lbrack 1 + \exp\{\mathbf{B}(t)^{\prime}{\mathbf{β}}\} \right\rbrack^{2}}{\exp\{\mathbf{B}(t)^{\prime}{\mathbf{β}}\}}.}
\end{aligned}$$

### Quasi-likelihood Loss

For quasi-likelihood loss, $\mathcal{L}_{3}$, in general we have
$$W_{\mathcal{L}_{3},g}(t;{\mathbf{β}}) = W_{3}(t;{\mathbf{β}}) = \mathbf{V}_{3}({\mathbf{β}})^{- 1}\mathbf{B}(t),\quad\mathbf{V}_{3}({\mathbf{β}}) = \int_{t = t_{1}}^{t_{2}}\mathbf{B}(t)\mathbf{B}(t)^{\prime}\left. \frac{\partial s(z)}{\partial z} \right|_{z = \mathbf{B}{(t)}^{\prime}{\mathbf{β}}}dt$$

#### Identity Link

For the identity link, $s(z) = z$, so
$$\frac{\partial s(z)}{\partial z} \equiv 1,$$ and
$$W_{3}(t;{\mathbf{β}}) = \mathbf{V}_{3}^{- 1}\mathbf{B}(t),$$ where
$$\mathbf{V}_{3} = \int_{t = t_{1}}^{t_{2}}\mathbf{B}(t)\mathbf{B}(t)^{\prime}dt.$$

#### Log Link

For the log link, $s(z) = \exp(z)$, so
$$\frac{\partial s(z)}{\partial z} = \exp(z),$$ and
$$W_{3}(t;{\mathbf{β}}) = \mathbf{V}_{3}({\mathbf{β}})^{- 1}\mathbf{B}(t),$$
where
$$\mathbf{V}_{3}({\mathbf{β}}) = \int_{t = t_{1}}^{t_{2}}\mathbf{B}(t)\mathbf{B}(t)^{\prime}\exp\left( \mathbf{B}(t)^{\prime}{\mathbf{β}} \right)dt.$$

#### Logit Link

For the logit link, $s(z) = \frac{\exp(z)}{1 + \exp(z)}$, so
$$\frac{\partial s(z)}{\partial z} = \frac{\exp(z)}{\left\{ 1 + \exp(z) \right\}^{2}},$$
and
$$W_{3}(t;{\mathbf{β}}) = \mathbf{V}_{3}({\mathbf{β}})^{- 1}\mathbf{B}(t),$$
where
$$\mathbf{V}_{3}({\mathbf{β}}) = \int_{t = t_{1}}^{t_{2}}\mathbf{B}(t)\mathbf{B}(t)^{\prime}\frac{\exp\left( \mathbf{B}(t)^{\prime}{\mathbf{β}} \right)}{\left\{ 1 + \exp\left( \mathbf{B}(t)^{\prime}{\mathbf{β}} \right) \right\}^{2}}dt.$$
