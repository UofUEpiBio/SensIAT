fit_PCORI_within_group_model_old <- function(
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
        delta_time = rlang::sym(glue::glue("delta({time.var})")),
        norm_time = rlang::sym(glue::glue("scale({time.var})")),
        norm_delta_time = rlang::sym(glue::glue("scale(delta({time.var}))"))
    )

    group.data2 <- filter(group.data, !!time.var <= End)

    u_hv <- group.data2 |> select(!!id.var) |> distinct() |> pull()
    N <- pull(summarize(group.data2, n_distinct(!!id.var)))


    ######   Andersen-Gill model stratifying by assessment number

    mf <- rlang::inject(!!outcome.var ~ !!id.var + !!time.var + !!rlang::f_rhs(intensity.covariates)) |>
        model.frame(data=filter(group.data, (!!time.var) <= !!End), na.action = na.pass) |>
        arrange(!!id.var, !!time.var)

    visit.number <- NULL
    model.data <- mf |>
        group_by(!!id.var) |>
        mutate(
            visit.number = seq_along(!!time.var)
        ) |>
        ungroup() |>
        # complete(!!id.var, visit.number, fill = tibble::lst(!!time.var := !!End)) |>
        group_by(!!id.var) |>
        arrange(!!id.var, visit.number) |>
        mutate(
            !!vars$prev_outcome    := lag(!!outcome.var, order_by = !!time.var),
            !!vars$prev_time       := lag(!!time.var, order_by =  !!time.var, default = 0),
            !!vars$delta_time      := !!time.var - lag(!!time.var, order_by =  !!time.var, default = 0),
            across(all_of(time.var), \(.)coalesce(., !!End))
        ) |>
        ungroup()

    centering.statistics <-
        summarize(
            ungroup(filter(model.data, !!time.var > 0, !is.na(!!outcome.var))),
            across(c(!!vars$time, !!vars$delta_time),
                   list(mean = ~ mean(.x, na.rm = TRUE),
                        sd = ~ sd(.x, na.rm = TRUE))
            )
        ) |>
        as.numeric() |>
        matrix(ncol = 2, byrow = TRUE) |>
        `dimnames<-`(list(c("time", "delta_time"), c("mean", "sd")))


    model.data <- model.data |>
        mutate(
            !!vars$norm_time := (!!vars$time - !!(centering.statistics[1,1]))/!!(centering.statistics[1, 2]),
            !!vars$norm_delta_time := (!!vars$delta_time - !!(centering.statistics[2, 1]))/!!(centering.statistics[2, 2])
        )

    ###########################   Model 1: Intensity model  ######################
    intensity.model <-
        rlang::inject(coxph(
            Surv(!!vars$prev_time,
                 !!time.var,
                 !is.na(!!outcome.var))~!!vars$prev_outcome+strata(visit.number),
            id = !!id.var,
            data = filter(model.data, !!vars$time > 0, !is.na(!!vars$prev_outcome))
        ))

    baseline_intensity_all =
        estimate_baseline_intensity(
            intensity.model = intensity.model,
            data = na.omit(model.data),
            variables = vars,
            bandwidth = control$intensity.bandwidth
        )$baseline_intensity




    #############  Model 2: Outcome model - single index model   ############
    outcome.formula <- rlang::inject(
        !!outcome.var~
            ns(!!vars$prev_outcome, df=3) +
            scale(!!vars$time) +
            scale(!!vars$delta_time) +
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
        intensity.model = intensity.model,
        outcome.model = outcome.model,
        outcome.model.centering = centering.statistics,
        data = model.data,
        variables = vars,
        End = End,
        influence = influence,
        alpha = alpha,
        coefficients = Beta$estimate,
        coefficient.variance = Beta$variance,
        control = control,
        base=base,
        V_inverse = V_inverse
    ), class = "PCORI_within_group_model",
    call = match.call())
}
globalVariables(
    c('..id..', '..time..', '..outcome..', '..prev_outcome..', '..prev_time..',
      '..delta_time..', '..norm_time..', '..norm_delta_time..', "..visit.number..")
)
