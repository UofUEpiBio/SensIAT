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
#' @param intensity.covariates = ~1 The extra covariates for the intensity model.
#' @param outcome.covariates = ~-1 The extra covariates for the outcome model.  The default removes the intercept term.
#' @param End The end time for this data analysis, we need to set the default value as the
#'           max value of the time
#' @param control control parameters for fitting the model.
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
#' fitted.trt.sim.numeric <-
#'     fit_PCORI_within_group_model(
#'         group.data = PCORI_example_data,
#'         outcome_modeler = PCORI_sim_outcome_modeler,
#'         alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'         id.var = Subject_ID,
#'         outcome.var = Outcome,
#'         time.var = Time,
#'         End = 830,
#'         knots = c(60,60,60,60,260,460,460,460,460),
#'         control = pcori_control('numeric')
#'     )
#' time.pw <- system.time({
#' fitted.trt.sim.pw <-
#'     fit_PCORI_within_group_model(
#'         group.data = PCORI_example_data,
#'         outcome_modeler = PCORI_sim_outcome_modeler,
#'         alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'         id.var = Subject_ID,
#'         outcome.var = Outcome,
#'         time.var = Time,
#'         End = 830,
#'         knots = c(60,60,60,60,260,460,460,460,460),
#'         control = pcori_control('piecewise')
#'     )
#' })
#' fitted.trt.sim.quadv <-
#'     fit_PCORI_within_group_model(
#'         group.data = PCORI_example_data,
#'         outcome_modeler = PCORI_sim_outcome_modeler,
#'         alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'         id.var = Subject_ID,
#'         outcome.var = Outcome,
#'         time.var = Time,
#'         End = 830,
#'         knots = c(60,60,60,60,260,460,460,460,460),
#'         control = pcori_control('quadv')
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
        intensity.covariates = ~1,
        outcome.covariates = ~-1,
        End = max({{time.var}}, na.rm = TRUE) + 1,
        control = pcori_control(),
        ...
){
    id.var <- ensym(id.var)
    outcome.var <- ensym(outcome.var)
    time.var <- ensym(time.var)

    outcome_modeler <- match.fun(outcome_modeler)
    End <- rlang::enexpr(End)
    End <- rlang::eval_tidy(End, data = group.data, env =parent.frame())

    vars <- list(
        id = id.var,
        time = time.var,
        outcome = outcome.var,
        prev_outcome = rlang::sym(glue::glue("lag({outcome.var})")),
        prev_time = rlang::sym(glue::glue("lag({time.var})")),
        delta_time = rlang::sym(glue::glue("Δ({time.var})")),
        norm_time = rlang::sym(glue::glue("scale({time.var})")),
        norm_delta_time = rlang::sym(glue::glue("scale(Δ({time.var}))"))
    )

    group.data2 <- filter(group.data, !!time.var <= End)

    u_hv <- group.data2 |> select(!!id.var) |> distinct() |> pull()
    N <- pull(summarize(group.data2, n_distinct(!!id.var)))


    ######   Andersen-Gill model stratifying by assessment number

    model.data <-
        rlang::inject(!!outcome.var ~ !!id.var + !!time.var + !!rlang::f_rhs(intensity.covariates)) |>
        model.frame(data=filter(group.data, (!!time.var) <= !!End), na.action = na.pass) |>
        arrange(!!id.var, !!time.var) |>
        group_by(!!id.var) |>
        mutate(
            visit.number = seq_along(!!time.var)
        ) |>
        ungroup() |>
        complete(!!id.var, visit.number, fill = tibble::lst(!!time.var := !!End)) |>
        group_by(!!id.var) |>
        arrange(!!id.var, visit.number) |>
        mutate(
            !!vars$prev_outcome    := lag(!!outcome.var, order_by = !!time.var),
            !!vars$prev_time       := lag(!!time.var, order_by =  !!time.var, default = 0),
            !!vars$delta_time      := !!time.var - lag(!!time.var, order_by =  !!time.var, default = 0),
            across(all_of(time.var), coalesce, !!End)
        ) |>
        ungroup()

    centering.statistics <-
        summarize( ungroup(filter(model.data, !!time.var > 0, !is.na(!!outcome.var))),
            "mean({time.var})" := mean(!!time.var),
            "sd({time.var})" := sd(!!time.var),
            "mean({vars$delta_time}))" := mean(!!vars$delta_time),
            "sd({vars$delta_time})" := sd(!!vars$delta_time)
        )

    model.data <- model.data |>
        mutate(
            !!vars$norm_time := (!!vars$time - !!pull(centering.statistics, 1))/!!pull(centering.statistics, 2),
            !!vars$norm_delta_time := (!!vars$delta_time - !!pull(centering.statistics, 3))/!!pull(centering.statistics, 4)
        )
    ###########################   Model 1: Intensity model  ######################
    intensity.model <-
        rlang::inject(coxph(
            Surv(!!vars$prev_time,
                 !!time.var,
                 !is.na(!!outcome.var))~!!vars$prev_outcome+strata(visit.number),
            id = !!id.var,
            data = filter(model.data, !!time.var > 0)
        ))

    # gamma <- Int_model$coefficients # parameter lambda in lambda(t, O(t))
    # we need this gamma value




    #############  Model 2: Outcome model - single index model   ############
    outcome.formula <- rlang::inject(
        !!outcome.var~
            ns(!!vars$prev_outcome, df=3) +
            !!vars$norm_time +
            !!vars$norm_delta_time +
            !!rlang::f_rhs(outcome.covariates)
    )
    outcome.model <- outcome_modeler(outcome.formula, data = filter(model.data, !!time.var > 0, !is.na(!!outcome.var)))


    base <- SplineBasis(knots)
    V_inverse <- solve(GramMatrix(base))

    # Compute value of the influence function: -----------------------------
    influence <- purrr::map_dfr(alpha, function(alpha){
        model.data |>
            group_by(!!vars$id) |>
            group_modify(
                compute_influence_for_one_alpha_and_one_patient,
                alpha = alpha,
                variables = vars,
                intensity.model = intensity.model,
                outcome.model = outcome.model,
                base = base,
                control = control,
                centering = centering.statistics,
                ...,
                .keep=TRUE
            )
    })

    # Results
    Beta = influence |>
        group_by(alpha) |>
        summarize(
            term1 = list(term1 |> reduce(rbind) |> `rownames<-`(NULL)),
            term2 = list(term2 |> reduce(rbind) |> `rownames<-`(NULL)),
            IF = list(influence |> reduce(rbind) |> `rownames<-`(NULL))
        ) |>
        mutate(
            IF_ortho = purrr::map(IF, \(IF){IF %*% V_inverse}),
            estimate = purrr::map(IF_ortho, colMeans),
            variance = purrr::map2(IF_ortho, estimate, ~tcrossprod(t(.x) - .y)/(nrow(.x)^2))
        )


    structure(list(
        intensity_model = intensity.model,
        outcome.model = outcome.model,
        outcome.model.centering = centering.statistics,
        data = model.data,
        variables = vars,
        End = End,
        influence = influence,
        coefficients = Beta$estimate,
        coefficient_variance = Beta$variance,
        control = control,
        base=base,
        V_inverse = V_inverse
    ), class = "PCORI_within_group_model")
}


#' Control Parameters for Fitting the PCORI Within Group Model
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
#' pcori_control("quad", intensity.bandwidth = 30, tol = 1e-6)
pcori_control <-
function(
    integration.method = c('quadv', "linear", "numerical", "piecewise"),
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


#' Prediction method for `PCORI_within_group_model` objects
#'
#' @param object a `PCORI_within_group_model` object.
#' @param time The time(s) to evaluate at.
#' @param alpha The values for the sensitivity parameter to evaluate at.
#' @param spline_fn Function to to generate the values for the spline functions
#'      with knots given by `knots`.  Must accept time as the first argument and
#'      knots as the second.
#' @param spline_seq The time points used to approximate the integral in term 2
#'          of the influence function.
#' @param knots The knots used for the spline function.
#'          Includes both boundary and internal knots.
#' @param integration.method Method for integration when computing the second influence term.  See Details.
#' @param integral.resolution the number of points to use for numerical integration.
#'
#'  Evaluate the fitted model, `object`, at each combination of `time` and
#'  `alpha`.
#'
#'
#' @return A tibble/data.frame with the components: `time`, `alpha`, `mean`, `var`.
#' where `time` and `alpha` are the combinations of the respective input(s) and
#' `mean` and `var` are the estimated mean and variance of the response for the
#' given model.
#'
#'  @details
#'  For `integration.method` when computing the integral in term two of the
#'  influence function `linear` approximates the expected value as piece-wise
#'  linear and is much faster, `numerical` uses traditional numerical
#'  integration.
#'
#'
#'
#' @export
#'
#' @examples
#'
#' fitted.trt.sim <-
#'     fit_PCORI_within_group_model(
#'         group.data = PCORI_example_data,
#'         outcome_modeler = PCORI_sim_outcome_modeler,
#'         id.var = Subject_ID,
#'         outcome.var = Outcome,
#'         time.var = Time,
#'         End = 830
#'     )
#' time.pw <- system.time({
#' pred.pw <- predict(fitted.trt.sim, time = c(90, 180),
#'     alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'     intensity.bandwidth = 30,
#'     knots=c(60,60,60,60,260,460,460,460,460),
#'     integration.method = 'piecewise'
#' )
#' })
#' time.num <- system.time({
#' pred.num <- predict(fitted.trt.sim, time = c(90, 180),
#'     alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'     intensity.bandwidth = 30,
#'     knots=c(60,60,60,60,260,460,460,460,460),
#'     integration.method = 'numerical'
#' )
#' })
#' time.lin <- system.time({
#' pred.lin <- predict(
#'     fitted.trt.sim, time = c(90,180),
#'     alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
#'     knots=c(60,60,60,60,260,460,460,460,460),
#'     integration.method = 'linear'
#' )
#' })
#' pred.num$Beta_hat[[2]]
#' pred.lin$Beta_hat[[3]]
#' pred.num$term1[[2]] |> head()
#' pred.lin$term1[[3]] |> head()
#' pred.num$term2[[2]] |> head()
#' pred.lin$term2[[3]] |> head()
#'
#'
`predict.PCORI_within_group_model` <-
function(object, time, alpha,
         knots,
         intensity.bandwidth = NULL,
         integral.resolution = 1000,
         integration.method = c("linear", "numerical", "piecewise", 'quadv'),
         ...
         ){
    assert_that(
        is(object, 'PCORI_within_group_model'),
        is.numeric(time),
        is.numeric(alpha),
        is.numeric(knots),
        is.count(integral.resolution)
    )
    integration.method = match.arg(integration.method)
    a <- min(knots)
    b <- max(knots)
    base <- SplineBasis(knots)

    # Compute value of the influence function: -----------------------------
    influence <- purrr::map_dfr(alpha, function(alpha){
        object$data |>
            group_by(!!object$variables$id) |>
            group_modify(
                compute_influence_for_one_alpha_and_one_patient,
                alpha = alpha,
                object = object,
                base=base,
                integration.method = integration.method,
                numerical.ingtegration.resolution = integral.resolution,
                ...,
                .keep=TRUE
            )
    })


    B_t <- evaluate(base, time)

    # Results
    influence |>
        group_by(alpha) |>
        summarize(
            term1 = list(term1 |> reduce(rbind) |> `rownames<-`(NULL)),
            term2 = list(term2 |> reduce(rbind) |> `rownames<-`(NULL)),
            IF = list(influence |> reduce(rbind) |> `rownames<-`(NULL))
        ) |>
        mutate(
            Beta_hat = purrr::map(IF, colMeans),
            Var_beta = purrr::map2(IF, Beta_hat, ~tcrossprod(t(.x) - .y)/(nrow(.x)^2))
        ) |>
        mutate(
            # Time specific estimates.
            mean_t = purrr::map(Beta_hat, ~B_t %*% .x),
            var_t = purrr::map(Var_beta, ~B_t %*% .x %*% t(B_t))
        )
}


