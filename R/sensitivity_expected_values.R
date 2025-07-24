

#' Compute Conditional Expected Values based on Outcome Model
#'
#' @param model An object representing the output of the outcome model.
#' @param alpha The sensitivity parameter
#' @param new.data Data to compute conditional means for, defaults to the model frame for the fitted model.
#' @param ... passed onto methods.
#' @param y.max The maximum value of the outcome variable for the Poisson and Negative Binomial models.
#'          If omitted it is chosen from the quantile function for the distribution at `1-eps`.
#' @param eps The tolerance for the quantile function used to estimate `y.max`, default is `.Machine$double.eps`.
#'
#' @details
#' Compute the conditional expectations needed for predictions in the models.
#' Two additional values/expectations are computed:
#'
#' * `$E \big[ Y(t) \exp \{  \alpha Y(t) \}   | A(t)=1, \bar{O}(t) \big]$`, returned as `E_Yexp_alphaY`, and
#' * `$E \big[ \exp \{ \alpha Y(t) \} \  | A(t)=1, \bar{O}(t) \big]$`, returned as `E_exp_alphaY`.
#'
#'  For the methods shown here
#'
#' @return The `new.data` frame with additional columns `alpha`, `E_Yexp_alphaY`, and `E_exp_alphaY` appended.
#' @export
sensitivity_expected_values <-
function(
    model,
    alpha = 0, #< sensitivity parameter
    new.data = model.frame(model),
    ...){
    UseMethod('sensitivity_expected_values')
}

#' @describeIn sensitivity_expected_values (Gaussian) Linear Model method
#' The [stats::integrate] method is used to compute the conditional expectations.
#' @examples
#' model <- lm(mpg ~ as.factor(cyl)+disp+wt, data=mtcars)
#' sensitivity_expected_values(model, alpha= c(-0.3, 0, 0.3), new.data = mtcars[1:5, ])
#' @export
sensitivity_expected_values.lm <- function(model, alpha, new.data, ...){
    if(length(alpha) > 1L){
        return(
            purrr::list_rbind(
                purrr::map(alpha, \(alpha1)sensitivity_expected_values.lm(model=model, alpha=alpha1, new.data=new.data, ...))
            )
        )
    }
    if(nrow(new.data) > 1L){
        return(
            purrr::list_rbind(
                purrr::map(1:nrow(new.data), \(i)sensitivity_expected_values.lm(model=model, alpha=alpha, new.data=new.data[i, , drop=FALSE], ...))
            )
        )
    }

    beta <- coef(model)
    sigma <- summary(model)$sigma

    Xi <- model.matrix(model, new.data)
    mu <- predict(model, newdata = new.data, type = "response")

    # compute the conditional expectations
    pmf_estimator <- function(y){
        dnorm(y, mu, sd = sigma)
    }

    E_exp_alphaY <- stats::integrate(
        \(y)if_else(pmf_estimator(y) == 0, 0, exp(alpha*y)*pmf_estimator(y)),
        lower = -Inf,
        upper = Inf
    )$value
    E_Yexp_alphaY <- stats::integrate(
        \(y)if_else(pmf_estimator(y) == 0, 0, y*exp(alpha*y)*pmf_estimator(y)),
        lower = -Inf,
        upper = Inf
    )$value
    data.frame(
        new.data,
        E_Yexp_alphaY = E_Yexp_alphaY,
        E_exp_alphaY = E_exp_alphaY
    )
}


#' @describeIn sensitivity_expected_values Generalized Linear Model method
#' @examples
#' model <- glm(cyl ~ mpg+disp+wt, data=mtcars, family=poisson())
#' sensitivity_expected_values(model, alpha= c(-0.3, 0, 0.3), new.data = mtcars[1:5, ]) |>
#'     dplyr::mutate('E(y|alpha)' = .data$E_Yexp_alphaY/.data$E_exp_alphaY)
#' @export
sensitivity_expected_values.glm <- function(model, alpha, new.data, ..., y.max = NULL, eps=.Machine$double.eps){
    # Recursion
    if(length(alpha) > 1L){
        return(
            purrr::list_rbind(
                purrr::map(alpha, \(alpha1)sensitivity_expected_values.glm(model=model, alpha=alpha1, new.data=new.data, ...))
            )
        )
    }
    if(nrow(new.data) > 1L){
        return(
            purrr::list_rbind(
                purrr::map(1:nrow(new.data), \(i)sensitivity_expected_values.glm(model=model, alpha=alpha, new.data=new.data[i, , drop=FALSE], ...))
            )
        )
    }

    if(family(model)$family == "gaussian"){
        return(sensitivity_expected_values.lm(model, alpha, new.data, ...))
    } else
    if (family(model)$family == "binomial"){
        mu <- predict(model, newdata = new.data, type = "response")
        y <- 0:1
        pmf <- dbinom(y, size = 1, prob = mu)
        E_exp_alphaY <- sum(exp(alpha*y)*pmf)
        E_Yexp_alphaY <- sum(y*exp(alpha*y)*pmf)
    } else
    if (family(model)$family == "poisson"){
        mu <- predict(model, newdata = new.data, type = "response")
        if(is.null(y.max)){
            y.max <- qpois(1-eps,lambda = mu)
        }
        y <- seq(0, y.max)
        pmf <- dpois(0:y.max, lambda = mu)
        E_exp_alphaY <- sum(exp(alpha*y)*pmf)
        E_Yexp_alphaY <- sum(y*exp(alpha*y)*pmf)
    } else {
        stop("Model family not supported")
    }
    return(
        tibble(
            new.data,
            alpha = alpha,
            E_Yexp_alphaY = E_Yexp_alphaY,
            E_exp_alphaY = E_exp_alphaY
        )
    )
}

#' @describeIn sensitivity_expected_values Negative Binomial Model method
#' @export
sensitivity_expected_values.negbin <- function(model, alpha, new.data, ..., y.max = NULL, eps=.Machine$double.eps^(1/4)){
    if(length(alpha) > 1L){
        return(
            purrr::list_rbind(
                purrr::map(alpha, \(alpha1)sensitivity_expected_values.negbin(model=model, alpha=alpha1, new.data=new.data, ...))
            )
        )
    }
    if(nrow(new.data) > 1L){
        return(
            purrr::list_rbind(
                purrr::map(1:nrow(new.data), \(i)sensitivity_expected_values.negbin(model=model, alpha=alpha, new.data=new.data[i, , drop=FALSE], ...))
            )
        )
    }

    mu <- predict(model, newdata = new.data, type = "response")
    theta <- model$theta

    pmf_estimator <- function(y){
        dnbinom(y, mu = mu, size = theta)
    }


    if(is.null(y.max)){
        y.max <- qnbinom(1-.Machine$double.eps, mu=mu, size=theta)
    }

    y <- seq(0, y.max)
    pmf <- pmf_estimator(y)
    E_exp_alphaY=sum( exp(alpha*y)*pmf )
    E_Yexp_alphaY=sum( y*exp(alpha*y)*pmf )

    return(
        tibble(
            new.data,
            alpha = alpha,
            E_Yexp_alphaY = E_Yexp_alphaY,
            E_exp_alphaY = E_exp_alphaY
        )
    )
}


