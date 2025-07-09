globalVariables(c('alpha', add=TRUE))

cross_validate <- function(original.object, progress = interactive(), prune = TRUE){
    data <- original.object$data
    ids <- unique(pull(data, original.object$variables$id))


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
                spline.degree = original.object$base@order-1L
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
    }
    if(progress){
        pb$tick(0)
        on.exit(pb$terminate(), add = TRUE)
    }
    structure(
        map(ids, run_without),
        ids = ids
    )
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
#'     intensity.args=list(bandwidth = 30),
#'     knots = c(60,260,460),
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


## TODO:
## Give class to jackknife results object
## Implement the autoplot function.
# 4 plots
# 1. within group lamp plot
# 2. Within group jackknife plot (Adds error bars to point estimates)
# 3. Full model effect contour plot.
# 4. Full model jackknife plot (Adds contour significance lines to full model contour plot.)
