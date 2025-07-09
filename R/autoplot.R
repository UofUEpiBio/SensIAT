


#' Plot a `SensIAT_within_group_model` object
#'
#' This creates a lamp plot for a `SensIAT_within_group_model` object.
#' The horizontal axis represents time, and the vertical axis represents the
#' expected marginal outcome given the sensitivity parameter `alpha`.
#'
#' @param object A `SensIAT_within_group_model` object.
#' @param ... currently ignored
#' @return A `ggplot2` object.
#' @export
#' @examples
#' object <-
#'     fit_SensIAT_within_group_model(
#'         group.data = SensIAT_example_data,
#'         outcome_modeler = SensIAT_sim_outcome_modeler_fbw,
#'         id = Subject_ID,
#'         outcome = Outcome,
#'         time = Time,
#'         knots = c(60,260,460),
#'         End = 830,
#'         alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'         intensity.args=list(bandwidth=30)
#'     )
#' autoplot(object) +
#'     # Title not included
#'     ggplot2::ggtitle("SensIAT within group model") +
#'     # Nor are bounds on reasonable values of alpha
#'     ggplot2::geom_hline(yintercept = c(1.2, 3), linetype = "dotted", linewidth = 1.5)
autoplot.SensIAT_within_group_model <- function(object, ...) {
  lower <- object$base@knots[object$base@order]
  upper <- object$base@knots[length(object$base@knots) - object$base@order + 1]

  x <- seq(lower, upper, length.out = 100)
  df <- predict(object, time=x, type = "response") |>
      mutate(alpha_factor = factor(alpha))

  ggplot2::ggplot(data=df, ggplot2::aes(x = time, y = mean, col = alpha_factor)) +
    ggplot2::geom_line() +
    ggplot2::labs(
        x = rlang::as_string(object$variables$time),
        y = rlang::as_string(object$variables$outcome),
        col = expression(alpha)
    )
}
