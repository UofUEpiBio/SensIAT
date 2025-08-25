#' Give the Marginal Mean Estimate and its Estimated Asymptotic Variance
#'
#' Give the marginal mean model estimate
#'
#'
#' @param object SensIAT_within_group_model object
#' @param time Time points of interest
#' @param include.var Logical. If TRUE, the variance of the outcome is also returned
#' @param ... Currently ignored.
#' @param base A `SplineBasis` object used to evaluate the basis functions.
#'
#' @return
#' If include.var is TRUE, a `tibble` with columns time, mean, and var is returned.
#' otherwise if include.var is FALSE, only the mean vector is returned.
#' @export
#'
#' @examples
#' \donttest{
#' model <-
#'     fit_SensIAT_within_group_model(
#'         group.data = SensIAT_example_data,
#'         outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
#'         alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'         id = Subject_ID,
#'         outcome = Outcome,
#'         time = Time,
#'         End = 830,
#'         knots = c(60,260,460),
#'     )
#' predict(model, time = c(90, 180))
#' }
predict.SensIAT_within_group_model <-
function(object, time, include.var= TRUE, ..., base = object$base){
    B <- do.call(rbind, map(time, pcoriaccel_evaluate_basis, spline_basis = base))
    tmp <- purrr::map2(object$coefficients, object$coefficient.variance,
         function(beta, var_beta){
             mean <- as.vector(B %*% beta)
             if(!include.var) return(tibble(time, mean))
             var  <- apply(B, 1, function(b) t(b) %*% var_beta %*% b)
             tibble(time, mean, var)
         }
    )

    tibble(alpha = object$alpha, tmp) |>
        tidyr::unnest(tmp)
}
