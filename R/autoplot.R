


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

#' Plot estimates at given times for `SensIAT_withingroup_jackknife_results` objects
#'
#' Horizontal axis represents time, and the vertical axis represents the outcome
#' from the model. Point plotted is the mean estimate, and the error bars
#' show the standard error estimated from the jackknife.
#'
#' @param object A `SensIAT_withingroup_jackknife_results` object produced from
#'                      `SensIAT_jackknife`.
#' @param width Width of the dodge for position, default is half the minimum
#'              distance between time evaluation points.
#' @param ... Ignored.
#' @return A `ggplot2` object.
#' @export
#' @examples
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
#' autoplot(jackknife.estimates)
#' }
autoplot.SensIAT_withingroup_jackknife_results <- function(object, width = NULL, ...) {

    if(is.null(width)) {
        width <- min(diff(sort(unique(object$time))))/2
        # /n_distinct(object$alpha)
    }

    dodge <- ggplot2::position_dodge(width = width)

    ggplot2::ggplot(data=dplyr::mutate(object, alpha_factor = factor(alpha)),
                    ggplot2::aes(
                        x = time, y = mean, col = alpha_factor,
                        ymin = mean - sqrt(jackknife_var),
                        ymax = mean + sqrt(jackknife_var),
                        group= alpha_factor
                    )) +
        ggplot2::geom_point(size = 2, position = dodge) +
        ggplot2::geom_errorbar(position = dodge) +
        ggplot2::labs(
            x = rlang::as_string(attr(object, "original.object")$variables$time),
            y = expr(
                !!rlang::as_string(attr(object, "original.object")$variables$outcome)
                %+-% sigma),
            col = expression(alpha)
        )
}

#' Plot for estimated treatment effect for `SensIAT_fulldata_model` objects
#'
#' The horizontal and vertical axes represent the sensitivity parameter `alpha`
#' for the control and treatment groups, respectively. The contour plot shows
#' the estimated treatment effect at each combination of `alpha` values.
#'
#' @param object A `SensIAT_fulldata_model` object.
#' @param time Time at which to plot the estimates.
#' @param ... Additional arguments passed to `predict`.
#'
#' @return A `ggplot2` object.
#' @export
#' @examples
#' full.object <-
#'     fit_SensIAT_fulldata_model(
#'         data = ARC_data,
#'         trt = Trt == 'home_visits',
#'         outcome_modeler = SensIAT_sim_outcome_modeler,
#'         id = elig_pid,
#'         outcome = Asthma_control,
#'         time = time,
#'         knots = c(60, 260, 460),
#'         alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'         End = 830
#'     )
#' autoplot(full.object, time = 180) +
#'      ggplot2::scale_fill_brewer(palette = "Spectral", guide = 'coloursteps')
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


#' Plot for estimated treatment effect for `SensIAT_fulldata_jackknife_results` objects
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
autoplot.SensIAT_fulldata_jackknife_results <-
  function(object, ..., include.rugs = NA) {
      rslt <- ggplot2::ggplot(data = object |>
                                  mutate(
                                      lower_95 = mean_effect + qnorm(0.025) * sqrt(mean_effect_jackknife_var),
                                      upper_95 = mean_effect + qnorm(0.975) * sqrt(mean_effect_jackknife_var),
                                      plot_value = pmax(lower_95, pmin(upper_95, 0))
                                  ),
                              ggplot2::aes(x = .data$alpha_control,
                                           y = .data$alpha_treatment,
                                           z = plot_value)) +
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
