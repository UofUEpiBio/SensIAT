#' Title
#'
#' @param data Data for evaluation of the model. Should match the data used to fit the intensity and outcome models.
#' @param alpha Sensitivity parameter, a vector of values.
#' @param knots Location of spline knots. If a `SplineBasis` object is provided, it is used directly.
#' @param intensity.model The assessment time intensity model.
#' @param outcome.model The observed effects model.
#' @param spline.degree The degree of the spline basis, default is 3 (cubic splines).
#' @param ... Additional arguments passed to `compute_influence_terms`.
#'
#' @returns a list with the fitted model, including the coefficients and their variances for each alpha value.
#' @export
#'
#' @examples
#' data_with_lags <- SensIAT_example_data |>
#'     group_by(Subject_ID) |>
#'     mutate(
#'         Prev_Outcome = lag(Outcome, default = NA_real_),
#'         Prev_time = lag(Time, default = NA_real_),
#'         Delta_Time = Time - Prev_time
#'
#'     )
#'
#' intensity.model <- coxph(Surv(Outcome))
#'
SensIAT_fit_marginal_model <-
function(data, alpha,
         knots,
         intensity.model,
         outcome.model,
         spline.degree = 3L,
         ...){

    if(is(knots, 'SplineBasis')){
        base <- knots
    } else {
        knots <- c(
            rep(head(knots,1), spline.degree),
            knots,
            rep(tail(knots, 1), spline.degree)
        )
        base <- SplineBasis(knots, order=spline.degree+1L)
    }
    V_inverse <- solve(GramMatrix(base))


    needed.vars <- c( all.vars(terms(intensity.model))
                    , all.vars(terms(outcome.model))
                    )
    stopifnot(all(needed.vars %in% names(data)))

    for_each_alpha <- \(a, ...){
        compute_influence_terms(
            data,
            base = base,
            alpha = a,
            outcome.model = outcome.model,
            intensity.model = intensity.model,
            ...
        )
    }
    influence.terms <- purrr::map(alpha, for_each_alpha)

    # Results
    Beta = map(influence.terms, \(IT){
        uncorrected.beta_hat <- (colSums(IT$term1) + colSums(IT$term2))/length(IT$id)
        estimate <- as.vector(V_inverse %*% uncorrected.beta_hat)
        variance <- tcrossprod(V_inverse %*% (t(IT$term1 + IT$term2) - uncorrected.beta_hat))/c(length(IT$id)^2)
        list(estimate = estimate, variance = variance)
    })


    structure(list(
        models = list(
            intensity = intensity.model,
            outcome = outcome.model
        ),
        data = data,
        influence = influence.terms,
        alpha = alpha,
        coefficients = map(Beta, getElement, 'estimate'),
        coefficient.variance = map(Beta, getElement, 'variance'),
        influence.args = list(...),
        base=base,
        V_inverse = V_inverse
    ), class = "SensIAT_marginal_outcome_model",
    call = match.call(expand.dots = TRUE)
    )
}
