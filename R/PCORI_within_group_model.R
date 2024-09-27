globalVariables(c('visit.number', 'term1', 'term2', 'IF', 'IF_ortho',
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
#'          Preferably passed in as a single tibble that internally is
#'          subsetted/filtered as needed.
#' @param outcome_modeler A separate function that may be swapped out to switch
#'          between negative-binomial, single index model, or another we will
#'          dream up in the future.
#' @param knots knot locations for defining the spline basis.
#' @param id.var The variable that identifies the patient.
#' @param outcome.var The variable that contains the outcome.
#' @param time.var The variable that contains the time.
#' @param alpha The sensitivity parameter.
#' @param intensity.covariates A formula representing modifications to the intensity model.
#' @param outcome.covariates A formula representing modifications to the outcome model.  The default removes the intercept term.
#' @param End The end time for this data analysis, we need to set the default value as the
#'           max value of the time
#' @param integration.tolerance The tolerance for the integration.
#' @param intensity.bandwidth The bandwidth for the intensity model kernel.
#' @param ... add parameters as needed or use this to pass forward into the
#'          outcome_modeler.
#'
#'
#' @return
#'      Should return everything needed to define the fit of the model.
#'      This can then be used for producing the estimates of mean, variance,
#'      and in turn treatment effect.
#'
#' @export
#'
#' @examples
#'
#' model <-
#'     fit_PCORI_within_group_model(
#'         group.data = PCORI_example_data,
#'         outcome_modeler = PCORI_sim_outcome_modeler,
#'         alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'         id.var = Subject_ID,
#'         outcome.var = Outcome,
#'         time.var = Time,
#'         End = 830,
#'         knots = c(60,60,60,60,260,460,460,460,460),
#'     )
#'
fit_PCORI_within_group_model <- function(
        group.data,
        outcome_modeler,
        knots,
        id.var,
        outcome.var,
        time.var,
        alpha = 0,
        intensity.covariates = ~.,
        outcome.covariates = ~.-1,
        End = max({{time.var}}, na.rm = TRUE) + 1,
        integration.tolerance = .Machine$double.eps^(1/3),
        intensity.bandwidth = NULL,
        ...
){
    id.var <- ensym(id.var)
    outcome.var <- ensym(outcome.var)
    time.var <- ensym(time.var)

    mf <- rlang::inject(!!outcome.var ~ !!id.var + !!time.var + !!rlang::f_rhs(intensity.covariates)) |>
        model.frame(data=filter(group.data, (!!time.var) <= !!End), na.action = na.pass) |>
        arrange(!!id.var, !!time.var)

    outcome_modeler <- match.fun(outcome_modeler)
    End <- rlang::enexpr(End)
    End <- rlang::eval_tidy(End, data = group.data, env =parent.frame())

    vars <- list(
        outcome = outcome.var,
        id = id.var,
        time = time.var
    )

    group.data2 <- filter(group.data, !!time.var <= End)

    u_hv <- group.data2 |> select(!!id.var) |> distinct() |> pull()
    N <- pull(summarize(group.data2, n_distinct(!!id.var)))


    ######   Andersen-Gill model stratifying by assessment number ------


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
        group_by(..id..) |>
        arrange(..id.., ..visit_number..) |>
        mutate(
            ..prev_outcome..    := lag(..outcome.., order_by = ..time..),
            ..prev_time..       := lag(..time.., order_by =  ..time.., default = 0),
            ..delta_time..      := ..time.. - lag(..time.., order_by =  ..time.., default = 0)
        ) |>
        ungroup()
    followup_data <- data_all_with_transforms |>
        filter(..time.. > 0, !is.na(..prev_outcome..))

    ######   Intensity model  ##################################################
    intensity.formula <- update.formula(
        Surv(..prev_time.., ..time..,  !is.na(..outcome..))
        ~ ..prev_outcome.. + strata(..visit_number..)
        , intensity.covariates
    )
    intensity.model <- coxph(intensity.formula,id = ..id..,data = followup_data)

    baseline_intensity_all =
        estimate_baseline_intensity(
            intensity.model = intensity.model,
            data = followup_data,
            variables = list(prev_outcome = sym("..prev_outcome..")),
            bandwidth = intensity.bandwidth
        )
    attr(intensity.model, 'bandwidth') <- baseline_intensity_all$bandwidth
    attr(intensity.model, 'kernel') <- baseline_intensity_all$kernel



    ######   Outcome model #####################################################
    outcome.formula <- update.formula(
        ..outcome..~
            ns(..prev_outcome.., df=3) +
            scale(..time..) +
            scale(..delta_time..)
        , outcome.covariates
    )
    outcome.model <- outcome_modeler(outcome.formula, data = followup_data)

    base <- SplineBasis(knots)
    V_inverse <- solve(GramMatrix(base))

    # Compute value of the influence function: -----------------------------
    influence.terms <- purrr::map(alpha,\(a){
        compute_influence_terms(
            data_all_with_transforms,
            base = base,
            alpha = a,
            outcome.model = outcome.model,
            intensity_coef = coef(intensity.model),
            baseline_intensity_all = baseline_intensity_all$baseline_intensity,
            tol = integration.tolerance
        )
    })

    # Results
    Beta = map(influence.terms, \(IT){
        uncorrected.beta_hat <- (colSums(IT$term1) + colSums(IT$term2))/length(IT$id)
        estimate <- as.vector(V_inverse %*% uncorrected.beta_hat)
        variance <- tcrossprod(V_inverse %*% (t(IT$term1 + IT$term2) - uncorrected.beta_hat))/c(length(IT$id)^2)
        list(estimate = estimate, variance = variance)
    })


    structure(list(
        intensity.model = structure(intensity.model,
                additional.covariates = intensity.covariates),
        outcome.model = structure(outcome.model,
                additional.covariates = intensity.covariates),
        data = mf,
        variables = vars,
        End = End,
        influence = influence.terms,
        alpha = alpha,
        coefficients = map(Beta, getElement, 'estimate'),
        coefficient.variance = map(Beta, getElement, 'variance'),
        intensity.bandwidth = intensity.bandwidth,
        base=base,
        V_inverse = V_inverse
    ), class = "PCORI_within_group_model",
    call = match.call())
}



