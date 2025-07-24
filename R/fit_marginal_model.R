#' Title
#'
#' @param data Data for evaluation of the model. Should match the data used to fit the intensity and outcome models.
#' @param id The subject identifier variable in the data. Lazy evaluation is used, so it can be a symbol or a string.
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
#' # Note: example takes approximately 30 seconds to run.
#' \donttest{
#' library(survival)
#' library(dplyr)
#' library(splines)
#' # Create followup data with lags
#' # added variables `..prev_time..`, `..delta_time..` and `..prev_outcome..`
#' # have special interpretations when computing the influence.
#' data_with_lags <- SensIAT_example_data |>
#'     dplyr::group_by(Subject_ID) |>
#'     dplyr::mutate(
#'         ..prev_outcome.. = dplyr::lag(Outcome, default = NA_real_, order_by = Time),
#'         ..prev_time.. = dplyr::lag(Time, default = 0, order_by = Time),
#'         ..delta_time.. = Time - dplyr::lag(.data$Time, default = NA_real_, order_by = Time)
#'     )
#'
#' # Create the observation time intensity model
#' intensity.model <-
#'     coxph(Surv(..prev_time.., Time, !is.na(Outcome)) ~ ..prev_outcome.. + strata(Visit),
#'     data = data_with_lags |> dplyr::filter(.data$Time > 0))
#'
#' # Create the observed outcome model
#' outcome.model <-
#'     SensIAT_sim_outcome_modeler(
#'         Outcome ~ ns(..prev_outcome.., df=3) + ..delta_time.. - 1,
#'         id = Subject_ID,
#'         data = data_with_lags |> filter(Time > 0))
#'
#' # Fit the marginal outcome model
#' mm <- SensIAT_fit_marginal_model(
#'     data = data_with_lags,
#'     id = Subject_ID,
#'     alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'     knots = c(60, 260, 460),
#'     intensity.model = intensity.model,
#'     time.vars = c('..delta_time..'),
#'     outcome.model = outcome.model)
#' }
SensIAT_fit_marginal_model <-
function(data,
         id,
         alpha,
         knots,
         outcome.model,
         intensity.model,
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

    for_each_alpha <- \(a){
        compute_influence_terms(
            data=data, id={{id}},
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
