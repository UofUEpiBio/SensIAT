
#' Fit the PCORI model for a single data set.
#'
#' Fit the PCORI model for both treatment and control groups.
#'
#' This function is agnostic to whether it is being given the original or a
#' bootstrap replication.
#'
#' @param data the full data set.
#' @param trt  an expression that determine what is treated as the treatment.
#'              Everything not treatment is considered control.
#' @param End
#' @param ... Passed on to fit_within_group_model.
#'
#' @return a list with class `PCORI-fulldata-fitted-model` with two components,
#'      `control` and `treatment`, each of which is an indepently fitted
#'      `PCORI-within-group-fitted-model` fit with the fit_within_group_model
#'      function.
#' @export
#'
#' @examples
#'
#' fitted.arc.nb.model <-
#'     fit_PCORI_fulldata_model(
#'         data = ARC_data,
#'         Trt = trt=='home_visits',
#'         # Passed to group model
#'         outcome_modeler = pcori_nb_outcome_modeler  #< not yet defined.
#'     )
#' fitted.arc.sim.model <-
#'     fit_PCORI_fulldata_model(
#'         data = ARC_data,
#'         Trt = trt=='home_visits',
#'         # Passed to group model
#'         outcome_modeler = PCORI_sim_outcome_modeler  #< not yet defined.
#'     )
#'
fit_PCORI_fulldata_model <- function(data, trt, ...){
    trt <- rlang::enquo(trt) #< diffuses evaluation for tidy expressions

    structure(list(
        control   = fit_PCORI_within_group_model(filter(data, !(!!trt)), ...),
        treatment = fit_PCORI_within_group_model(filter(data, {{trt}} ) , ...))
    ), class = 'PCORI_fulldata_model')
}

#' Prediction method for `PCORI_fulldata_model` objects.
#'
#' @param object an object of class `PCORI_fulldata_model`
#' @param time  time points to evaluate at
#' @param alpha alpha values to evaluate at.
#'
#' @description
#' For each combination of `time` and `alpha` estimate the mean response and
#' variance for each group as well as estimate the mean treatment effect and
#' variance.
#'
#'
#' @return a tibble/data.frame with the following components.
#'  * `time`, the time point
#' @export
#'
#' @examples
`predict.PCORI_fulldata_model` <-
    function(object, time, alpha){

        control.predicted    <- predict(object$control  , time, alpah)
        treatment.predicted  <- predict(object$treatment, time, alpha)

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

