#' Outcome Modeler for PCORI Negative Binomial Model.
#'
#' @param formula The outcome model formula
#' @param data The data to fit the outcome model to.
#' @param ... Currently ignored, included for future compatibility.
#'
#' @return Object of class `PCORI::Single-index-outcome-model` which contains the outcome model portion.
#' @export
PCORI_nb_outcome_modeler <- function(formula, 
                                     data, 
                                     ...){
  mf <- model.frame(formula, data = data)
  # Xi <- model.matrix(formula, data = mf)
  # Yi <- model.response(mf)
  model <- glm.nb(formula, data=data)
  structure(
      model,
      class = c('PCORI::outcome-model', 'PCORI::Negative-binomial-model', class(model))
      )
}

#' @export
`model.frame.PCORI::Negative-binomial-model` <-
    function(formula, data=NULL, ...){
        if(is.null(data))
            data <- formula$data
        NextMethod('model.frame', data=data, ...)
    }
#' @export
`model.matrix.PCORI::Negative-binomial-model` <-
    function(object, data = model.frame(object), ...){
        model.matrix(terms(object), data = data, ...)
    }#' @export
`formula.PCORI::Negative-binomial-model` <-
    function(x, ...){
        as.formula(terms(x))
    }

#' @export
`predict.PCORI::Negative-binomial-model` <-
    function( object
            , newdata = NULL
            , type = c('response', 'terms')
            , ...){
        if(is.null(newdata)) newdata = model.frame(object)
        type = match.arg(type)

        frame <-


        predict(object$formula, data = data, ...)

        if(type == 'terms'){}
    }


###### need to change ######

Cond_mean_fn_nb <-
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
            xb = x %*% beta, y = y,
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

# nuisance parameter: Upper bound of Y = default setting is observed maximum value of Y (2 * maximum(Y))

#' @export
`pcori_conditional_means.PCORI::Negative-binomial-model` <-
function(
    model,
    alpha,
    # gamma,
    new.data = model.frame(model),
    ...
    )
{
    assert_that(
        is(model, 'PCORI::Negative-binomial-model'),
        is.numeric(alpha)
    )
    if(length(alpha) > 1){
        return(
            purrr::map_dfr(
                alpha,
                `pcori_conditional_means.PCORI::Negative-binomial-model`,
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
        temp <- Cond_mean_fn_nb(alpha,
                                     # X = Xi,
                                     # Y = Yi,
                                     # x = Xi_new[k,,drop=FALSE],
                                     # beta = model$coef,
                                     # bandwidth = model$bandwidth,
                                     # kernel = attr(model, 'kernel'
                                                   )
        )

        E_Y_past[k] <- temp$E_Y_past
        E_exp_alphaY[k] <- temp$E_exp_alphaY
        E_Yexp_alphaY[k] <- temp$E_Yexp_alphaY

    }

    tibble(new.data, alpha, E_Y_past, E_exp_alphaY, E_Yexp_alphaY)
}
