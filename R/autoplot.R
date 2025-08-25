


#' Plot a `SensIAT_within_group_model` Object
#'
#' This creates a line plot for a `SensIAT_within_group_model` object.
#' The horizontal axis represents time, and the vertical axis represents the
#' expected marginal outcome given the sensitivity parameter `alpha`.
#'
#' @param object A `SensIAT_within_group_model` object.
#' @param ... currently ignored
#' @return A `ggplot2` object.
#' @export
#' @examples
#' # Note: example takes a few seconds to run.
#' \donttest{
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
#' ggplot2::autoplot(object) +
#'     # Title not included
#'     ggplot2::ggtitle("SensIAT within group model") +
#'     # Nor are bounds on reasonable values of alpha
#'     ggplot2::geom_hline(yintercept = c(1.2, 3), linetype = "dotted", linewidth = 1.5)
#' }
autoplot.SensIAT_within_group_model <- function(object, ...) {
  lower <- object$base@knots[object$base@order]
  upper <- object$base@knots[length(object$base@knots) - object$base@order + 1]

  x <- seq(lower, upper, length.out = 100)
  df <- predict(object, time=x, type = "response") |>
      mutate(alpha_factor = factor(alpha))

  ggplot2::ggplot(data=df, ggplot2::aes(x = .data$time, y = .data$mean, col = .data$alpha_factor)) +
    ggplot2::geom_line() +
    ggplot2::labs(
        x = rlang::as_string(object$variables$time),
        y = rlang::as_string(object$variables$outcome),
        col = expression(alpha)
    )
}

#' Plot Estimates at Given Times for `SensIAT_withingroup_jackknife_results` Objects
#'
#' Horizontal axis represents time, and the vertical axis represents the outcome
#' from the model. Point plotted is the mean estimate, and the error bars
#' show the 95% confidence interval using the variance estimated from the jackknife.
#'
#' @param object A `SensIAT_withingroup_jackknife_results` object produced from
#'                      `SensIAT_jackknife`.
#' @param width Width of the dodge for position, default is half the minimum
#'              distance between time evaluation points.
#' @param ... Ignored.
#' @return A `ggplot2` object.
#' @export
#' @examples
#' # Note: fitting the jackknife is computationally expensive,
#' #       so this example is here for reference.
#' \dontrun{
#' fitted <-
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
#' jackknife.estimates <- SensIAT_jackknife(fitted, time = c(90, 180, 270, 360, 450))
#' ggplot2::autoplot(jackknife.estimates)
#' }
autoplot.SensIAT_withingroup_jackknife_results <- function(object, width = NULL, ...) {

    if(is.null(width)) {
        width <- min(diff(sort(unique(object$time))))/2
        # /n_distinct(object$alpha)
    }

    dodge <- ggplot2::position_dodge(width = width)

    ggplot2::ggplot(data=dplyr::mutate(object, alpha_factor = factor(alpha)),
                    ggplot2::aes(
                        x = .data$time,
                        y = .data$mean,
                        col = .data$alpha_factor,
                        ymin = mean + qnorm(0.025) * sqrt(.data$jackknife_var),
                        ymax = mean + qnorm(0.975) * sqrt(.data$jackknife_var),
                        group= .data$alpha_factor
                    )) +
        ggplot2::geom_point(size = 2, position = dodge) +
        ggplot2::geom_errorbar(position = dodge) +
        ggplot2::labs(
            x = rlang::as_string(attr(object, "original.object")$variables$time),
            y = substitute(expression(outcome %+-% sigma),
                           attr(object, "original.object")$variables),
            col = expression(alpha)
        )
}

#' Plot for Estimated Treatment Effect for `SensIAT_fulldata_model` Objects
#'
#' The horizontal and vertical axes represent the sensitivity parameter `alpha`
#' for the control and treatment groups, respectively. The contour plot shows
#' the estimated treatment effect at each combination of `alpha` values.
#'
#' @param object A `SensIAT_fulldata_model` object.
#' @param time Time at which to plot the estimates.
#' @param include.rugs If `TRUE`, adds rugs indicating the locations where the
#'    sensitivity was evaluated to the plot. If `FALSE`, no rugs are added.
#'    If `NA`, rugs are added only if the number of distinct values of
#'    `alpha_control` and `alpha_treatment` is less than or equal to 10.
#' @param ... Additional arguments passed to `predict`.
#'
#' @return A `ggplot2` object.
#' @export
#' @examples
#' \donttest{
#' full.object <-
#'     fit_SensIAT_fulldata_model(
#'         data = SensIAT_example_fulldata,
#'         trt = Treatment_group == 'treatment',
#'         outcome_modeler = SensIAT_sim_outcome_modeler,
#'         id = Subject_ID,
#'         outcome = Outcome,
#'         time = Time,
#'         knots = c(60, 260, 460),
#'         alpha = c(-0.6, -0.3, 0, 0.3, 0.6)
#'     )
#' ggplot2::autoplot(full.object, time = 180)
#' }
autoplot.SensIAT_fulldata_model <- function(object, time, include.rugs = NA, ...) {
    df <- predict(object, time, ...)

    rslt <- ggplot2::ggplot(data = df,
                            ggplot2::aes(x = .data$alpha_control,
                                         y = .data$alpha_treatment,
                                         z = .data$mean_effect)) +
        ggplot2::labs(
            x = expression(alpha[control]),
            y = expression(alpha[treatment]),
            fill = "Treatment Effect"
        )

    if(rlang::is_installed("metR")){
      rslt <- rslt +
        metR::geom_contour_fill() +
        metR::scale_fill_divergent(mid="white",low="forestgreen",high="orange")
    } else {
      rlang::inform("Package 'metR' is recomended for creating 'SensIAT' contour plots, please install it.",
                    .frequency = 'once')
      rslt <- rslt +
        ggplot2::geom_contour_filled()
    }

    if (length(time)> 1){
        rslt <- rslt + ggplot2::facet_wrap(~time, scales = 'free')
    } else {
        rslt <- rslt + ggplot2::ggtitle(paste("Treatment Effect at Time =", time))
    }

    if(isFALSE(include.rugs)) return(rslt)
    if(is.na(include.rugs)
       && ( n_distinct(df$alpha_control) > 10
          || n_distinct(df$alpha_treatment) > 10
          )
        ) return(rslt)
    rslt + ggplot2::geom_rug(sides = 'bl')
}


#' Plot for Estimated Treatment Effect for `SensIAT_fulldata_jackknife_results` Objects
#'
#' The horizontal and vertical axes represent the sensitivity parameter `alpha`
#' for the control and treatment groups, respectively. The plot shows
#' at each combination of `alpha` values zero if the 95% confidence interval
#' contains zero, otherwise the bound of the confidence interval that is closest
#' to zero.
#'
#' @param object A `SensIAT_fulldata_jackknife_results` object.
#' @param ... Additional arguments passed to `predict`.
#' @param include.rugs If `TRUE`, adds rugs to the plot. If `FALSE`, no rugs are added.
#' When `NA`, rugs are added only if the number of distinct values of `alpha_control`
#' and `alpha_treatment` is less than or equal to 10.
#'
#' @export
#' @examples
#' # Note: fitting the jackknife is computationally expensive,
#' #       so this example is here for reference.
#' \dontrun{
#' full.object <-
#'     fit_SensIAT_fulldata_model(
#'         data = SensIAT_example_fulldata,
#'         trt = Treatment_group == 'treatment',
#'         outcome_modeler = SensIAT_sim_outcome_modeler,
#'         id = Subject_ID,
#'         outcome = Outcome,
#'         time = Time,
#'         knots = c(60, 260, 460),
#'         alpha = c(-0.6, -0.3, 0, 0.3, 0.6)
#'     )
#' jk.full.model <- jackknife(full.object, time = 180)
#' ggplot2::autoplot(jk.full.model)
#' }
autoplot.SensIAT_fulldata_jackknife_results <-
  function(object, ..., include.rugs = NA) {
      rslt <- ggplot2::ggplot(data = object |>
                                  mutate(
                                      lower_95 = mean_effect + qnorm(0.025) * sqrt(.data$mean_effect_jackknife_var),
                                      upper_95 = mean_effect + qnorm(0.975) * sqrt(.data$mean_effect_jackknife_var),
                                      plot_value = pmax(.data$lower_95, pmin(.data$upper_95, 0))
                                  ),
                              ggplot2::aes(x = .data$alpha_control,
                                           y = .data$alpha_treatment,
                                           z = .data$plot_value)) +
          ggplot2::labs(
              x = expression(alpha[control]),
              y = expression(alpha[treatment]),
              fill = expression(delta),
              caption = expression(
                paste(delta == 0, italic(' if '), 0 %in% 'CI'['95%'],
                      plain(' otherwise '), delta ,
                      plain(' is the bound of the 95% CI closest to zero.')))
          )


    if(rlang::is_installed("metR")){
      rslt <- rslt +
        metR::geom_contour_fill() +
        metR::scale_fill_divergent(mid="white",low="forestgreen",high="orange")
    } else {
      rlang::inform("Package 'metR' is recomended for creating 'SensIAT' contour plots, please install it.",
                    .frequency = 'once')
      rslt <- rslt +
          ggplot2::geom_contour_filled()
    }
      time <- unique(object$time)
      if (length(time)> 1){
          rslt <- rslt + ggplot2::facet_wrap(~time, scales = 'free')
      } else {
          rslt <- rslt + ggplot2::ggtitle(paste("Treatment Effect at Time =", time))
      }

      if(isFALSE(include.rugs)) return(rslt)
      if(is.na(include.rugs)
         && ( n_distinct(object$alpha_control) > 10
              || n_distinct(object$alpha_treatment) > 10
         )
      ) return(rslt)


    rslt + ggplot2::geom_rug(sides = 'bl')
}
