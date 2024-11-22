
#' @describeIn fit_SensIAT_within_group_model
#' Fit the sensitivity analysis for both treatment and control groups.
#'
#' @param data the full data set.
#' @param trt  an expression that determine what is treated as the treatment.
#'              Everything not treatment is considered control.
#' @param ... common arguments passed to `fit_SensIAT_within_group_model`.
#'
#' @return a list with class `SensIAT-fulldata-fitted-model` with two components,
#'      `control` and `treatment`, each of which is an independently fitted
#'      `SensIAT-within-group-fitted-model` fit with the fit_within_group_model
#'      function.
#' @export
fit_SensIAT_fulldata_model <- function(data, trt, ...){
    trt <- rlang::enquo(trt) #< diffuses evaluation for tidy expressions

    structure(list(
        control   = fit_SensIAT_within_group_model(filter(data, !as.logical(!!trt)), ...),
        treatment = fit_SensIAT_within_group_model(filter(data,  as.logical(!!trt)), ...)
    ), class = 'SensIAT_fulldata_model')
}

#' @describeIn predict.SensIAT_within_group_model
#' For each combination of `time` and `alpha` estimate the mean response and
#' variance for each group as well as estimate the mean treatment effect and
#' variance.
#' @export
`predict.SensIAT_fulldata_model` <-
    function(object, time, ...){

        control.predicted    <- predict(object$control  , time, ...)
        treatment.predicted  <- predict(object$treatment, time, ...)

        full_join(
            control.predicted,
            treatment.predicted,
            by= c('time', 'alpha'),
            suffix = c('_control', '_treatment')
        ) |>
            mutate(
                mean_effect = mean_treatment - mean_control,
                var_effect = var_treatment + var_control
            )
    }
globalVariables(c('mean_effect', 'mean_treatment', 'mean_control', 'var_effect', 'var_treatment', 'var_control'))
