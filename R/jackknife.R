globalVariables(c('alpha', add=TRUE))

cross_validate <- function(original.object, progress = interactive()){
    data <- original.object$data
    ids <- unique(pull(data, original.object$variables$id))


    run_without <- function(id){
        if(progress)on.exit(try(pb$tick(), silent = TRUE))
        data |> filter(!!original.object$variables$id != id) |>
            fit_SensIAT_within_group_model(
                outcome_modeler = attr(original.object, 'call')$outcome_modeler,
                knots = original.object$base@knots,
                id.var = !!original.object$variables$id,
                outcome.var = !!original.object$variables$outcome,
                time.var = !!original.object$variables$time,
                alpha = original.object$alpha,
                intensity.covariates = attr(original.object$intensity.model, 'additional.covariates'),
                outcome.covariates = attr(original.object$outcome.model, 'additional.covariates'),
                End = original.object$End
            ) |>
            prune_bootstrap_replication()
    }

    if(progress && rlang::is_installed("progress")){
        pb <- progress::progress_bar$new(
            format = "  cross-validation [:bar] :current/:total(:percent) eta: :eta",
            total = length(ids)
        )
    }
    if(progress) pb$tick(0)

    map(ids, run_without)
}


#' Estimate response with jackknife resampling
#'
#' @param original.object A SensIAT_within_group_model object.
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
#' original.object <-
#' fit_SensIAT_within_group_model(
#'     group.data = SensIAT_example_data,
#'     outcome_modeler = SensIAT_sim_outcome_modeler,
#'     alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'     id.var = Subject_ID,
#'     outcome.var = Outcome,
#'     time.var = Time,
#'     intensity.bandwidth = 30,
#'     knots = c(60,60,60,60,260,460,460,460,460),
#'     End = 830
#' )
#' jackknife.estimates <- SensIAT_jackknife(original.object, time = c(90, 180, 270, 360, 450))
#' }
SensIAT_jackknife <- function(original.object, time, ...){
    replications <- cross_validate(original.object)

    estimates <- map(replications, predict.SensIAT_within_group_model, time=time,
                     include.var=FALSE, base = original.object$base)
    estimates |> bind_rows(.id='.rep') |>
        group_by(alpha, time) |>
        summarize(
            # n = n(),
            jackknife_mean = mean(mean),
            jackknife_var = (n()-1)/n() * sum((mean-mean(mean))^2),
        , .groups='drop')
}
