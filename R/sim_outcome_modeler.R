#' Outcome Modeler for PCORI Single Index Model.
#'
#' @param formula The outcome model formula
#' @param data The data to fit the outcome model to.
#'             Should only include follow-up data, i.e. time > 0.
#' @param kernel The kernel to use for the outcome model.
#' @param method The optimization method to use for the outcome model, either "optim", "nlminb", or "nmk".
#' @param ... Currently ignored, included for future compatibility.
#'
#' @return Object of class `PCORI::Single-index-outcome-model` which contains the outcome model portion.
#' @export
PCORI_sim_outcome_modeler <-
function(formula, data, kernel = "K2_Biweight", method = "nmk", id = ..id.., ...){
  id <- ensym(id)
  mf <- rlang::inject(model.frame(formula, data = data, id = !!id))
  Xi <- model.matrix(formula, data = mf)

  Yi <- model.response(mf)

  initial = estimate_starting_coefficients(Xi, Yi)

  val <- SIDR_Ravinew(X = Xi, Y = Yi, index_ID= mf[['(id)']],
                      initial=initial,
                      kernel = kernel,
                      method = method,
                      ...)
  structure(
      append(
        val,
          list(
              frame = mf,
              data = data
          )
      ),
      class = c('PCORI::outcome-model', 'PCORI::Single-index-outcome-model'),
      kernel = kernel,
      terms = terms(mf))
}

#' @export
`model.frame.PCORI::Single-index-outcome-model` <-
    function(formula, data=NULL, ...){
        if(is.null(data))
            data <- formula$data
        NextMethod('model.frame', data=data, ...)
    }
#' @export
`model.matrix.PCORI::Single-index-outcome-model` <-
    function(object, data = model.frame(object), ...){
        model.matrix(terms(object), data = data, ...)
    }
#' @export
`formula.PCORI::Single-index-outcome-model` <-
    function(x, ...){
        as.formula(terms(x))
    }
#' @export
`coef.PCORI::Single-index-outcome-model` <-
    function(x, ...)x$coef

#' @export
`predict.PCORI::Single-index-outcome-model` <-
    function( object
            , newdata = NULL
            , type = c('response', 'terms')
            , ...){
        if(is.null(newdata)) newdata = model.frame(object)
        type = match.arg(type)

        frame <-


        predict(object$formula, data = newdata, ...)

        if(type == 'terms'){}
    }

estimate_starting_coefficients <- function(X,Y, eps = 1e-7){
    X <- as.matrix(X)
    Y <- as.vector(Y)

    number_n <- dim(X)[1]
    number_p <- dim(X)[2]

    Y.CP <- outer(Y, Y, "<=")

    # centralizing covariates
    X.cs <- t(t(X)-colMeans(X))

    # calculating m(y)=\E[X_i 1(Y_i\leq y)]
    # m.y <- (t(X)-colMeans(X)) %*% Y.CP/number_n
    # calculating K=\E[m(Y_i)m(Y_i)^T]
    m.y <- crossprod(X.cs, Y.CP)/number_n
    Km <- tcrossprod(m.y)/number_n

    eigen(solve(var(X) + eps*diag(number_p), Km))$vectors[,1]
}


K2_Biweight_kernel <- function(x, h){15/16*(1-(x/h)^2)^2 * (abs(x) <= h)}
K4_Biweight_kernel <- function(x, h){105/64*(1-3*((x/h)^2))*(1-(x/h)^2)^2 * (abs(x) <= h) }

NW_new <-
function(Xb, Y, xb, y, h, kernel = "K2_Biweight"){

    if(kernel == "dnorm"){
        K <- function(x, h){dnorm(x/h, 0, 1)} # Gaussian
    } else if(kernel == "K2_Biweight"){
        K <- function(x, h){15/16*(1-(x/h)^2)^2 * (abs(x) <= h)} # K2_biweight
    } else if(kernel=="K4_Biweight"){
        K <- function(x, h){105/64*(1-3*((x/h)^2))*(1-(x/h)^2)^2 * (abs(x) <= h) }# K4_biweight
    }

    Kxb <- sapply(xb, function(x, Xb) K(Xb-x, h), Xb=Xb)

    Ylty <- sapply(y, function(x, Y) 1*(Y <= x), Y=Y)

    denom <- colSums(Kxb)

    fyxb <- (denom!=0)*crossprod(Kxb, Ylty)/(denom + (denom==0))

    return(fyxb)

}

Cond_mean_fn_single2 <-
    function( alpha #< sensitivity parameter
            , X     #< Matrix of covariates for all observations, including the spline basis as well as other covariates such as lag(time) and lag(outcome)
            , Y     #< Outcome vector for all observations
            , x     #< vector of covariates for the observation of interest
            , beta
            , bandwidth
            , ...  #< for passing kernel forward
            ){


        y <- sort(unique(Y))

        # conditional distribution
        #start <- Sys.time()
        Fhat <- pcoriaccel_NW(
            Xb = X %*% beta, Y = Y,
            xb = x %*% beta, y_seq = y,
            h = bandwidth,
            ...)
        #end <- Sys.time()
        #end - start

        # density function
        Fhat1 <- c(0, Fhat[1:(length(y) - 1)])
        pmf <- Fhat - Fhat1

        # Question: Are we assuming Y is finite with support range_y or are we approximating an integral here?
        E_exp_alphaY <- sum( exp(alpha*y)*pmf )

        E_Yexp_alphaY <- sum( y*exp(alpha*y)*pmf )

        E_Y_past <- E_Yexp_alphaY/E_exp_alphaY

        return(list(
            E_Y_past = E_Y_past,
            E_exp_alphaY = E_exp_alphaY,
            E_Yexp_alphaY = E_Yexp_alphaY
        ))

    }



#' @export
`pcori_conditional_means.PCORI::Single-index-outcome-model` <-
function(
    model,
    alpha,
    # gamma,
    new.data = model.frame(model),
    ...
    )
{
    assert_that(
        is(model, 'PCORI::Single-index-outcome-model'),
        is.numeric(alpha)
    )
    if(length(alpha) > 1){
        return(
            purrr::map_dfr(
                alpha,
                `pcori_conditional_means.PCORI::Single-index-outcome-model`,
                model = model, new.data = new.data,
                ...
            )
        )
    }
    if (nrow(new.data)==0) return(mutate(
        new.data,
        alpha = alpha,
        E_Y_past = numeric(0),
        E_exp_alphaY = numeric(0),
        E_Yexp_alphaY = numeric(0)
    ))

    Xi <- model.matrix(terms(model), model$data)
    Yi <- model.response(model.frame(model))
    for(var in setdiff(all.vars(terms(model)), tbl_vars(new.data)))
        new.data[[var]] <- NA
    Xi_new <- model.matrix(terms(model), data=new.data)

    if(nrow(Xi_new)==0) return(mutate(
        new.data,
        alpha = alpha,
        E_Y_past = NA_real_,
        E_exp_alphaY = NA_real_,
        E_Yexp_alphaY = NA_real_
    ))

    E_Y_past <- numeric(nrow(Xi_new))
    E_exp_alphaY <- numeric(nrow(Xi_new))
    E_Yexp_alphaY <- numeric(nrow(Xi_new))

    for(k in 1:nrow(Xi_new)){
        # df_k <- new.data[k, ]
        # x = model.matrix(terms(model), data = df_k)
        temp <- Cond_mean_fn_single2(alpha,
                                     X = Xi,
                                     Y = Yi,
                                     x = Xi_new[k,,drop=FALSE],
                                     beta = model$coef,
                                     bandwidth = model$bandwidth,
                                     kernel = attr(model, 'kernel')
        )

        E_Y_past[k] <- temp$E_Y_past
        E_exp_alphaY[k] <- temp$E_exp_alphaY
        E_Yexp_alphaY[k] <- temp$E_Yexp_alphaY

    }

    tibble(new.data, alpha, E_Y_past, E_exp_alphaY, E_Yexp_alphaY)
}
