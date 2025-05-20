globalVariables(c('..visit_number..', 'term1', 'term2', 'IF', 'IF_ortho',
                  'estimate', 'Beta_hat', 'Var_beta'))

#' Produce fitted model for group (treatment or control)
#'
#' Produces a fitted model that may be used to produce estimates of mean and
#' variance for the given group.
#'
#' This function should be agnostic to whether it is being provided a
#' treatment or control group.
#'
#' @param group.data The data for the group that is being analyzed.
#'          Preferably passed in as a single `tibble` that internally is
#'          subsetted/filtered as needed.
#' @param outcome_modeler A separate function that may be swapped out to switch
#'          between negative-binomial, single index model, or another we will
#'          dream up in the future.
#' @param knots knot locations for defining the spline basis.
#' @param id The variable that identifies the patient.
#' @param outcome The variable that contains the outcome.
#' @param time The variable that contains the time.
#' @param alpha The sensitivity parameter.
#' @param End The end time for this data analysis, we need to set the default value as the
#'          max value of the time.
#' @param intensity.args A list of optional arguments for intensity model.
#'          See the Intensity Arguments section.
#' @param outcome.args parameters as needed passed into the `outcome_modeler`.
#'          One special element may be `'model.modifications'` which, if present,
#'          should be a formula that will be used to modify the outcome model per, [update.formula].
#' @param influence.args A list of optional arguments used when computing the influence.
#'         See the Influence Arguments section.
#' @param spline.degree The degree of the spline basis.
#'
#' @return
#'  Should return everything needed to define the fit of the model.
#'  This can then be used for producing the estimates of mean, variance,
#'  and in turn treatment effect.  For the full data model a list with two
#'  models one each for the treatment and control groups.
#'
#' @export
#'
#' @examples
#' \donttest{
#' model <-
#'     fit_SensIAT_within_group_model(
#'         group.data = SensIAT_example_data,
#'         outcome_modeler = SensIAT_sim_outcome_modeler,
#'         alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'         id = Subject_ID,
#'         outcome = Outcome,
#'         time = Time,
#'         End = 830,
#'         knots = c(60,260,460),
#'     )
#' model <-
#'     fit_SensIAT_within_group_model(
#'         group.data = SensIAT_example_data |> mutate(Outcome=Outcome*6),
#'         outcome_modeler = MASS::glm.nb,
#'         alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'         id = Subject_ID,
#'         outcome = Outcome,
#'         time = Time,
#'         End = 830,
#'         knots = c(60,260,460),
#'     )
#' }
fit_SensIAT_within_group_model <- function(
        group.data,
        outcome_modeler,
        id,
        outcome,
        time,
        knots,
        alpha = 0,
        End = NULL,
        intensity.args = list(),
        outcome.args = list(),
        influence.args = list(),
        spline.degree = 3
){
    ###### Input clean and capture -------------------------------------------
    # Variables
    id.var <- ensym(id)
    outcome.var <- ensym(outcome)
    time.var <- ensym(time)
    vars <- list(
        outcome = outcome.var,
        id = id.var,
        time = time.var
    )

    assert_that(is.numeric(knots))
    if(anyDuplicated(knots)){
        rlang::warn("Duplicate knots may be ignored.")
        knots <- unique(knots)
    }
    if(is.unsorted(knots)){
        rlang::warn("Knots are not sorted.")
        knots <- sort(knots)
    }
    knots <- c(
        rep(head(knots,1), spline.degree),
        knots,
        rep(tail(knots, 1), spline.degree)
    )

    # Pass through Argument Lists
    intensity.args <- match.names(intensity.args, c("model.modifications", 'bandwidth'), FALSE)
    outcome.args   <- match.names(outcome.args  , c("model.modifications"))
    influence.args <- match.names(influence.args, c("tolerance"))

    # Function
    outcome_modeler <- match.fun(outcome_modeler)
    if(is.null(End)){
        End <- rlang::expr(max({{time}}, na.rm = TRUE) + 1) %>%
            rlang::eval_tidy(data = group.data, env = parent.frame())
    }

    ###### Create Model Frame -----------------------------------------------
    outcome.extra.vars <- all.vars(outcome.args$model.modifications) |>
        setdiff(c('..id..', '..time..', '..prev_outcome..', '..delta_time..', '..visit_number..', '..outcome..'))
    intensity.extra.vars <- all.vars(intensity.args$model.modifications) |>
        setdiff(c('..id..', '..time..', '..prev_outcome..', '..delta_time..', '..visit_number..', '..outcome..'))

    mf <- group.data |>
        filter((!!time.var) <= !!End) |>
        select(!!id.var, !!time.var, !!outcome.var,
               any_of(outcome.extra.vars),
               any_of(intensity.extra.vars)) |>
        arrange(!!id.var, !!time.var)

    data_all_with_transforms <- mf |>
        rename(
            ..id.. = !!id.var,
            ..time.. = !!time.var,
            ..outcome.. = !!outcome.var
        ) |>
        group_by(..id..) |>
        mutate(
            ..visit_number.. = seq_along(..time..)
        ) |>
        ungroup() |>
        complete(..id.., ..visit_number..,
                 fill = list(..time.. = End,..outcome.. = NA_real_)
        ) |>
        group_by(..id..) |>
        arrange(..id.., ..visit_number..) |>
        mutate(
            ..time..            := as.double(..time..),
            ..prev_outcome..    := lag(..outcome.., order_by = ..time..),
            ..prev_time..       := lag(..time.., order_by =  ..time.., default = 0),
            ..delta_time..      := ..time.. - lag(..time.., order_by =  ..time.., default = 0)
        ) |>
        ungroup()


    ######   Intensity model  ##################################################
    #' @section Intensity Arguments:
    #' The `intensity.args` list may contain the following elements:
    #'
    #' * **`model.modifications`** A formula that will be used to modify the intensity model from it's default, per [update.formula].
    #' * **`kernel`** The kernel function for the intensity model. Default is the Epanechnikov kernel.
    #' * **`bandwidth`** The bandwidth for the intensity model kernel.
    #'
    followup_data <- data_all_with_transforms |>
        filter(..time.. > 0, !is.na(..prev_outcome..))

    intensity.formula <- Surv(..prev_time.., ..time..,  !is.na(..outcome..)) ~ ..prev_outcome.. + strata(..visit_number..)
    if(!is.null(intensity.args$model.modifications))
        intensity.formula <- update.formula(intensity.formula, intensity.args$model.modifications)

    intensity.model <-  rlang::inject(
        coxph(intensity.formula, id = ..id.., data = followup_data,
                             !!!purrr::discard_at(intensity.args, c("bandwidth", "kernel", "model.modifications")))
        )
    # baseline_intensity_all <-
    #     estimate_baseline_intensity(
    #         intensity.model = intensity.model,
    #         data = followup_data,
    #         kernel = intensity.args$kernel %||% \(x) 0.75*(1 - (x)**2) * (abs(x) < 1),
    #         bandwidth = intensity.args$bandwidth
    #     )
    attr(intensity.model, 'bandwidth') <- intensity.args$bandwidth
    attr(intensity.model, 'kernel') <- intensity.args$kernel %||% \(x) 0.75*(1 - (x)**2) * (abs(x) < 1)
    # followup_data$baseline_intensity <- baseline_intensity_all$baseline_intensity

    ######   Outcome model #####################################################
    outcome.formula <-
        ..outcome..~ -1 +
            ns(..prev_outcome.., df=3) +
            scale(..time..) +
            scale(..delta_time..)
    if(!is.null(outcome.args$model.modifications))
        outcome.formula <- update.formula(outcome.formula, outcome.args$model.modifications)

    outcome.model <- rlang::inject(
        outcome_modeler(outcome.formula,
                        data = filter(followup_data, !is.na(..outcome..)),
                        !!!purrr::discard_at(outcome.args, "model.modifications"))
        )

    base <- SplineBasis(knots, order=spline.degree+1L)
    V_inverse <- solve(GramMatrix(base))

    # Compute value of the influence function: -----------------------------
    #' @section Influence Arguments:
    #'
    #' The `influence.args` list may contain the following elements:
    #'
    #' * **`method`** The method for integrating, adaptive or fixed quadrature. Default is `'adaptive'`.
    #' * **`tolerance`** The tolerance when using adaptive quadrature.
    #' * **`delta`** The bin width for fixed quadrature.
    #' * **`resolution`** alternative to `delta` by specifying the number of bins.
    #' * **`fix_discontinuity`** Whether to account for the discontinuity in the influence at observation times.
    influence.terms <- rlang::inject(purrr::map(alpha,\(a){
        compute_influence_terms(
            left_join(
                data_all_with_transforms,
                followup_data |> select('..id..', '..time..'),
                by=join_by('..id..', '..time..')
            ),
            base = base,
            alpha = a,
            outcome.model = outcome.model,
            intensity.model = intensity.model,
            !!!influence.args
        )
    }))

    # Results
    Beta = map(influence.terms, \(IT){
        uncorrected.beta_hat <- (colSums(IT$term1) + colSums(IT$term2))/length(IT$id)
        estimate <- as.vector(V_inverse %*% uncorrected.beta_hat)
        variance <- tcrossprod(V_inverse %*% (t(IT$term1 + IT$term2) - uncorrected.beta_hat))/c(length(IT$id)^2)
        list(estimate = estimate, variance = variance)
    })


    structure(list(
        models = list(
            intensity = intensity.model,
            outcome = outcome.model
        ),
        data = mf,
        variables = vars,
        End = End,
        influence = influence.terms,
        outcome_modeler = outcome_modeler,
        alpha = alpha,
        coefficients = map(Beta, getElement, 'estimate'),
        coefficient.variance = map(Beta, getElement, 'variance'),
        args = list(
            intensity = intensity.args,
            outcome = outcome.args,
            influence = influence.args
        ),
        base=base,
        V_inverse = V_inverse
    ), class = "SensIAT_within_group_model",
    call = match.call(expand.dots = TRUE))
}



