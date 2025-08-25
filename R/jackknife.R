globalVariables(c('alpha', add=TRUE))

cross_validate <- function(original.object, progress = interactive(), prune = TRUE){
    data <- original.object$data
    ids <- na.omit(unique(pull(data, original.object$variables$id)))


    run_without <- function(id){
        if(progress)on.exit(try(pb$tick(), silent = TRUE))
        replication <-
        data |> filter(!!original.object$variables$id != id) |>
            fit_SensIAT_within_group_model(
                outcome_modeler = original.object$outcome_modeler,
                id = !!original.object$variables$id,
                outcome = !!original.object$variables$outcome,
                time = !!original.object$variables$time,
                knots = unique(original.object$base@knots),
                alpha = original.object$alpha,
                End = original.object$End,
                intensity.args = original.object$args$intensity,
                outcome.args = original.object$args$outcome,
                influence.args = original.object$args$influence,
                spline.degree = original.object$base@order-1L,
                add.terminal.observations = FALSE
            )
        if(prune)
            replication <- prune(replication)
        replication$jackknife_excluded_id <- id
        return(replication)
    }

    if(progress && rlang::is_installed("progress")){
        pb <- progress::progress_bar$new(
            format = "  cross-validation [:bar] :current/:total(:percent) eta: :eta",
            total = length(ids)
        )
        pb$tick(0)
        on.exit(pb$terminate(), add = TRUE)
    }
    structure(
        map(ids, run_without),
        ids = ids
    )
}



#' Perform Jackknife resampling on an object.
#'
#' @param object An object to cross validate on.
#' @param ... Additional arguments passed to the method.
#'
#' @return A data frame of the jackknife resampling results.
#' @export
jackknife <- function(object, ...){
    UseMethod('jackknife')
}

#' @describeIn jackknife Perform jackknife resampling on a `SensIAT_within_group_model` object.
#' @export
jackknife.SensIAT_within_group_model <- function(object, ...){
    SensIAT_jackknife(object, ...)
}

#' @describeIn jackknife Perform jackknife resampling on a `SensIAT_fulldata_model` object.
#' @export
jackknife.SensIAT_fulldata_model <- function(object, ...){
    SensIAT_jackknife_fulldata(object, ...)
}


#' Estimate response with jackknife resampling
#'
#' @param object A SensIAT_within_group_model object.
#' @param time Time points for which to estimate the response.
#' @param ... currently ignored.
#'
#' @return
#' A `tibble` with columns alpha, time, jackknife_mean, and jackknife_var,
#' where jackknife_mean is the mean of the jackknife estimates and jackknife_var
#' is the estimated variances of the response at the given time points for the
#' specified alpha values.
#' @export
#'
#' @examples
#' \dontrun{
#' object <-
#' fit_SensIAT_within_group_model(
#'     group.data = SensIAT_example_data,
#'     outcome_modeler = SensIAT_sim_outcome_modeler,
#'     alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'     id = Subject_ID,
#'     outcome = Outcome,
#'     time = Time,
#'     intensity.args=list(bandwidth = 30),
#'     knots = c(60,260,460),
#'     End = 830
#' )
#' jackknife.estimates <- SensIAT_jackknife(object, time = c(90, 180, 270, 360, 450))
#' }
SensIAT_jackknife <- function(object, time, ...){
    replications <- cross_validate(object)

    summarize_jackknife_replications(replications, object, time, ...)
}
summarize_jackknife_replications <- function(replications, original.object, time, ...){
    original.estimates <- predict.SensIAT_within_group_model(original.object, time=time, ...)
    estimates <- map(replications, predict.SensIAT_within_group_model, time=time,
                     include.var=FALSE, base = original.object$base)
    estimates |> bind_rows(.id='.rep') |>
        group_by(alpha, time) |>
        summarize(
            # n = n(),
            jackknife_mean = mean(mean),
            jackknife_var = (n()-1)/n() * sum((mean-mean(mean))^2),
            , .groups='drop') |>
        ungroup() |>
        full_join(original.estimates, by=c('alpha', 'time')) |>
        add_class('SensIAT_withingroup_jackknife_results') |>
        structure(original.object = original.object)
}

#' @describeIn SensIAT_jackknife Estimate variance of the treatment effect with jackknife resampling for full data models.
#' @export
SensIAT_jackknife_fulldata <- function(object, time, ...){
    assertthat::assert_that(is.numeric(time), length(time) > 0,
                            msg = "time must be a numeric vector of length > 0")

    replications_treatment <- cross_validate(object$treatment)
    replications_control <- cross_validate(object$control)

    summary_treatment <- summarize_jackknife_replications(replications_treatment, object$treatment, time)
    summary_control <- summarize_jackknife_replications(replications_control, object$control, time)

    full_join(
        summary_treatment,
        summary_control,
        by=c('time'),
        suffix=c('_treatment', '_control'),
        relationship = "many-to-many"
    ) |> as_tibble() |>
        add_class('SensIAT_fulldata_jackknife_results') |>
        structure(
            original.object = object,
            summary_treatment = summary_treatment,
            summary_control = summary_control
        ) |>
        mutate(
            mean_effect = .data$mean_treatment - .data$mean_control,
            mean_effect_jackknife_var = .data$jackknife_var_treatment + .data$jackknife_var_control
        )
}
