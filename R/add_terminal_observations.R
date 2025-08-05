#' Add Terminal Observations to a Dataset
#'
#' This function adds terminal observations to a dataset.
#' For each subject given by `id`, if that subject has less than the maximum
#' number of observations, A row is added with the `end` time value, leaving all
#' other variables as `NA`.
#'
#' @param data A data frame containing the dataset.
#' @param id A variable in `data` that identifies the subject.
#' @param time A variable in `data` that identifies the time of the observation.
#' @param end The value to use for the `time` variable in the terminal observation.
#'            If end is less that the maximum in the dataset resulting data will
#'            be filtered such that `time` is less than or equal to `end`.
#' @return A data frame with terminal observations added.
#' @seealso [tidyr::complete()]
#' @examples
#' exdata <- tibble::tibble(
#'   patient = rep(1:3, 3:5),
#'   day = c(0, 30, 60,
#'           0, 30, 60, 90,
#'           0, 30, 60, 90, 120),
#'   value = TRUE
#' )
#' add_terminal_observations(exdata, patient, day)
#' @export
add_terminal_observations <-
function(data, id, time, end=max(pull(data, {{time}}))){
    id <- rlang::ensym(id)
    time <- rlang::ensym(time)

    force(end)

    data |>
        filter({{time}} <= end) |>
        group_by({{id}}) |>
        arrange({{time}}) |>
        mutate(
            ..visit_number.. = seq_along({{time}})
        ) |>
        ungroup() |>
        tidyr::complete(
            {{id}},
            ..visit_number..,
            fill = structure(list(end), names = rlang::as_name(time)),
            explicit = FALSE
            )|>
        group_by({{id}}) |>
        filter({{time}} != dplyr::lag({{time}}, default = -Inf)) |>
        select(-..visit_number..)
}
