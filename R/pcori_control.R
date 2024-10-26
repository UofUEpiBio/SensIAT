#' Control Parameters for Fitting the within Group Model
#'
#' @param integration.method        Method for integration when computing the second influence term.
#' @param intensity.bandwidth       The bandwidth for the intensity model.
#' @param resolution                The number of points to use for numerical integration.
#' @param resolution.within.period  The number of points to use for numerical integration within a period.
#' @param tol                       The tolerance for numerical integration.
#' @param ...                       Currently ignored.
#'
#' @return a list of control parameters.
#' @export
#'
#' @examples
#' pcori_control("piecewise", intensity.bandwidth = 30, resolution.within.period = 50)
#' pcori_control("numerical", intensity.bandwidth = 30, resolution = 1000)
#' pcori_control("quadv", intensity.bandwidth = 30, tol = 1e-6)
#' pcori_control("quadvcpp", intensity.bandwidth = 30, tol = 1e-6)
pcori_control <-
    function(
        integration.method = c('quadvcpp', 'quadv', "linear", "numerical", "piecewise"),
        intensity.bandwidth = NULL,
        resolution = 1000,
        resolution.within.period = 50,
        tol=.Machine$double.eps^(1/4),
        ...
    ){
        integration.method = match.arg(integration.method)
        assert_that(
            is.null(intensity.bandwidth) || is.numeric(intensity.bandwidth),
            is.count(resolution),
            is.count(resolution.within.period),
            rlang::is_scalar_double(tol)
        )
        lst(
            integration.method,
            intensity.bandwidth,
            resolution,
            resolution.within.period,
            tol,
            ...
        )
    }
