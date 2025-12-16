#' Extrapolate Variables from Last Observation
#'
#' Extrapolates variables from the last observation before a given time point
#' using linear extrapolation with specified slopes.
#'
#' @param target_time Numeric scalar. The time point at which to extrapolate values.
#' @param data A data.frame containing the observations, with one row per time point.
#'   Must be sorted by the time variable.
#' @param time_var Character string. The name of the time variable in `data`.
#' @param slopes A named numeric vector. Names should match column names in `data`
#'   for variables to extrapolate, and values are the slopes (rate of change per unit time).
#'   Variables not in `slopes` will be carried forward unchanged from the last observation.
#' @param strict Logical. If `TRUE` (default), requires that `target_time` is greater than
#'   or equal to the last observed time. If `FALSE`, allows extrapolation to earlier times.
#'
#' @return A single-row data.frame with extrapolated values at `target_time`.
#'
#' @details
#' The function finds the last observation in `data` where the time is less than or
#' equal to `target_time`, then extrapolates variables using the formula:
#'
#' \deqn{x_{extrapolated} = x_{last} + (target\_time - time_{last}) \times slope}
#'
#' For variables not specified in `slopes`, the last observed value is used (slope = 0).
#'
#' @keywords internal
#'
#' @examples
#' # Create example data
#' df <- data.frame(
#'     time = c(0, 1, 2, 3),
#'     x1 = c(10, 12, 14, 16),
#'     x2 = c(5, 5.5, 6, 6.5),
#'     id = 1
#' )
#'
#' # Extrapolate to time = 5 with known slopes
#' extrapolate_from_last_observation(
#'     target_time = 5,
#'     data = df,
#'     time_var = "time",
#'     slopes = c(x1 = 2, x2 = 0.5)
#' )
#' # Expected: x1 = 16 + (5-3)*2 = 20, x2 = 6.5 + (5-3)*0.5 = 7.5
#'
#' # Extrapolate with only some variables having slopes
#' extrapolate_from_last_observation(
#'     target_time = 4,
#'     data = df,
#'     time_var = "time",
#'     slopes = c(x1 = 2)
#' )
#' # Expected: x1 = 16 + (4-3)*2 = 18, x2 = 6.5 (carried forward), id = 1
#'
#' @keywords internal
extrapolate_from_last_observation <- function(
  target_time,
  data,
  time_var,
  slopes = NULL,
  strict = TRUE
) {
    # Input validation
    if (!rlang::is_scalar_double(target_time)) {
        rlang::abort("`target_time` must be a single numeric value")
    }

    if (!is.data.frame(data)) {
        rlang::abort("`data` must be a data.frame")
    }

    if (!rlang::is_string(time_var) || !time_var %in% names(data)) {
        rlang::abort(sprintf("`time_var` must be a character string naming a column in `data`. Got: %s", time_var))
    }

    times <- data[[time_var]]
    if (!is.numeric(times)) {
        rlang::abort(sprintf("The time variable `%s` must be numeric", time_var))
    }

    if (length(times) == 0) {
        rlang::abort("`data` must have at least one row")
    }

    # Find the last observation before or at target_time
    valid_indices <- which(times <= target_time)

    if (length(valid_indices) == 0) {
        if (strict) {
            rlang::abort(sprintf(
                "No observations found at or before target_time = %g (earliest time is %g)",
                target_time, min(times)
            ))
        } else {
            # Use the first observation
            period <- 1
        }
    } else {
        period <- max(valid_indices)
    }

    # Get the last observation
    last_obs <- data[period, , drop = FALSE]
    last_time <- last_obs[[time_var]]

    # Calculate time difference
    delta_time <- target_time - last_time

    # Initialize result with last observation
    result <- last_obs

    # Apply slopes to specified variables
    if (!is.null(slopes)) {
        if (!is.numeric(slopes) || is.null(names(slopes))) {
            rlang::abort("`slopes` must be a named numeric vector")
        }

        # Check that all slope names are in data
        missing_vars <- setdiff(names(slopes), names(data))
        if (length(missing_vars) > 0) {
            rlang::abort(sprintf(
                "Variables in `slopes` not found in `data`: %s",
                paste(missing_vars, collapse = ", ")
            ))
        }

        # Apply linear extrapolation for each variable with a slope
        for (var_name in names(slopes)) {
            slope_value <- slopes[[var_name]]
            result[[var_name]] <- last_obs[[var_name]] + delta_time * slope_value
        }
    }

    # Update the time variable to target_time
    result[[time_var]] <- target_time

    # Add attributes with metadata
    attr(result, "source_period") <- period
    attr(result, "source_time") <- last_time
    attr(result, "delta_time") <- delta_time

    return(result)
}


#' Extrapolate Multiple Time Points from Last Observation
#'
#' Vectorized version of [extrapolate_from_last_observation()] that extrapolates
#' to multiple time points at once.
#'
#' @inheritParams extrapolate_from_last_observation
#' @param target_times Numeric vector. Multiple time points at which to extrapolate values.
#'
#' @return A data.frame with one row per target time, containing extrapolated values.
#'
#' @examples
#' df <- data.frame(
#'     time = c(0, 1, 2, 3),
#'     x1 = c(10, 12, 14, 16),
#'     x2 = c(5, 5.5, 6, 6.5)
#' )
#'
#' extrapolate_from_last_observation_multiple(
#'     target_times = c(3.5, 4, 5),
#'     data = df,
#'     time_var = "time",
#'     slopes = c(x1 = 2, x2 = 0.5)
#' )
#'
#' @keywords internal
extrapolate_from_last_observation_multiple <- function(
  target_times,
  data,
  time_var,
  slopes = NULL,
  strict = TRUE
) {
    if (!is.numeric(target_times)) {
        rlang::abort("`target_times` must be a numeric vector")
    }

    results <- purrr::map(
        target_times,
        \(t) extrapolate_from_last_observation(
            target_time = t,
            data = data,
            time_var = time_var,
            slopes = slopes,
            strict = strict
        )
    )

    purrr::list_rbind(results)
}
