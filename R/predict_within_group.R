#' Predict mean and variance of the outcome for a PCORI within-group model
#'
#' @param object PCORI_within_group_model object
#' @param time Time points of interest
#' @param include.var Logical. If TRUE, the variance of the outcome is also returned
#'
#' @return
#' If include.var is TRUE, a tibble with columns time, mean, and var is returned.
#' otherwise if include.var is FALSE, only the mean vector is returned.
#' @export
#'
#' @examples
predict.PCORI_within_group_model <-
function(object, time, include.var= TRUE, ..., base = object$base){
    B <- evaluate(base, time)
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
