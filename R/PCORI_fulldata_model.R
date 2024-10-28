
#' Fit the sensitivity analysis model for a single data set.
#'
#' Fit the sensitivity analysis for both treatment and control groups.
#'
#' This function is agnostic to whether it is being given the original or a
#' bootstrap replication.
#'
#' @param data the full data set.
#' @param trt  an expression that determine what is treated as the treatment.
#'              Everything not treatment is considered control.
#' @param ... Passed on to fit_within_group_model.
#'
#' @return a list with class `SensIAT-fulldata-fitted-model` with two components,
#'      `control` and `treatment`, each of which is an independently fitted
#'      `SensIAT-within-group-fitted-model` fit with the fit_within_group_model
#'      function.
#' @export
#' @examples
#' \dontrun{
#' model <-
#'     fit_SensIAT_within_group_model(
#'         group.data = SensIAT_example_data,
#'         outcome_modeler = SensIAT_sim_outcome_modeler,
#'         alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'         id.var = Subject_ID,
#'         outcome.var = Outcome,
#'         time.var = Time,
#'         End = 830,
#'         knots = c(60,60,60,60,260,460,460,460,460),
#'     )
#' }
fit_SensIAT_fulldata_model <- function(data, trt, ...){
    trt <- rlang::enquo(trt) #< diffuses evaluation for tidy expressions

    structure(list(
        control   = fit_SensIAT_within_group_model(filter(data, !as.logical(!!trt)), ...),
        treatment = fit_SensIAT_within_group_model(filter(data,  as.logical(!!trt)), ...)
    ), class = 'SensIAT_fulldata_model')
}

#' Prediction method for `SensIAT_fulldata_model` objects.
#'
#' @param object an object of class `SensIAT_fulldata_model`
#' @param time  time points to evaluate at
#' @param ... additional arguments passed to `predict.SensIAT_within_group_model`
#'
#' @description
#' For each combination of `time` and `alpha` estimate the mean response and
#' variance for each group as well as estimate the mean treatment effect and
#' variance.
#'
#'
#' @return a `tibble`/`data.frame` with the following components.
#'  * `time`, the time point
#' @export
#' @examples
#' \dontrun{
#' model <-
#'     fit_SensIAT_within_group_model(
#'         group.data = SensIAT_example_data,
#'         outcome_modeler = SensIAT_sim_outcome_modeler,
#'         alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'         id.var = Subject_ID,
#'         outcome.var = Outcome,
#'         time.var = Time,
#'         End = 830,
#'         knots = c(60,60,60,60,260,460,460,460,460),
#'     )
#' predict(model, time = c(90, 180))
#' }
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
