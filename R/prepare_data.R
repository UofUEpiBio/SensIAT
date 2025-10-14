#' Prepare Data for Sensitivity Analysis with Irregular Assessment Times
#'
#' This function prepares the data for SensIAT analysis by transforming it into
#' a format suitable for the SensIAT models.
#'
#' @param data A data frame containing the data to be prepared.
#' @param id.var The variable in `data` that identifies the subject.
#' @param time.var The variable in `data` that identifies the time of the observation.
#' @param outcome.var The variable in `data` that contains the outcome of interest.
#' @param End The end time for the analysis. Observations with time greater than `End` will be filtered out.
#' @param add.terminal.observations Logical indicating whether to add terminal observations to the data (`TRUE`), or terminal observations have already been added (`FALSE`).
#' @return A data frame with the following transformations:
#' \itemize{
#' \item Data filtered to time less than or equal to `End`.
#' \item Observations are arranged by `id.var` and `time.var`.
#' \item Terminal observations added if `add.terminal.observations` is `TRUE`,
#'       with `..time..` set to `End` and `..outcome..` set to `NA`, if the
#'       subject has less observations than the maximum number of observations.
#' \item New variables created: \itemize{
#'     \item `..id..` aliases `id.var`,
#'     \item `..time..` aliases `time.var`,
#'     \item `..outcome..` aliases `outcome.var`,
#'     \item `..visit_number..` is the visit number within each subject derived from `time.var`,
#'     \item `..prev_outcome..`, i.e. lag-outcome,  the outcome from the previous visit,
#'     \item `..prev_time..`, i.e. lag-time, the time from the previous visit,
#'     \item `..delta_time..`, the difference in time between the current and previous visit.
#'     }
#' }
#' @export
#' @examples
#'
#' prepare_SensIAT_data( SensIAT_example_data, Subject_ID, Time, Outcome, 830)
#'
#' exdata <- tibble::tibble(ID=rep(1:2, c(3,5)),
#'                          Time=c(0, 30, 60,
#'                                 0, 30, 60, 90, 120),
#'                          Outcome=floor(runif(8, 1, 100)))
#'
#' prepare_SensIAT_data(exdata, ID, Time, Outcome, 120)
#'
prepare_SensIAT_data <-
function(data,
         id.var, time.var, outcome.var,
         End,
         add.terminal.observations = TRUE){
    id.var <- ensym(id.var)
    outcome.var <- ensym(outcome.var)
    time.var <- ensym(time.var)

    if(add.terminal.observations && anyNA(data)) # Ensure no NAs in the data
        rlang::abort("Data contains missing values, cannot add terminal observations.")
    data_all_with_transforms <- data |>
        filter((!!time.var) <= !!End) |>
        arrange(!!id.var, !!time.var) |>
        mutate(
            ..id.. = !!id.var,
            ..time.. = !!time.var,
            ..outcome.. = !!outcome.var
        ) |>
        group_by(..id.., !!id.var) |>
        mutate(
            ..visit_number.. = seq_along(..time..) - 1L
        ) |>
        ungroup()
    if(add.terminal.observations){
        data_all_with_transforms <- data_all_with_transforms |>
            complete(..id.., ..visit_number..,
                     fill = tibble::lst(
                                ..time.. = End,..outcome.. = NA_real_,
                                !!time.var := End, !!outcome.var := NA_real_
                                 )
            ) |>
            mutate(!!id.var := .data$..id..)
    }
    data_all_with_transforms <- data_all_with_transforms |>
        group_by(..id..) |>
        arrange(..id.., ..visit_number..) |>
        mutate(
            ..time..            := as.double(..time..),
            ..prev_outcome..    := lag(..outcome.., order_by = ..time..),
            ..prev_time..       := lag(..time.., order_by =  ..time.., default = 0),
            ..delta_time..      := ..time.. - lag(..time.., order_by =  ..time.., default = 0)
        ) |>
        ungroup()

    data_all_with_transforms
}
