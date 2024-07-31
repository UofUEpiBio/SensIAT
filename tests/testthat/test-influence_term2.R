test_that("Compute Influence term2 old vs. new methods", {
    model.data <- PCORI_example_data |>
        group_by(Subject_ID) |>
        arrange(Time) |>
        mutate(
            prev_time = lag(Time),
            prev_outcome = lag(Outcome),
            delta_time = Time - lag(Time),
            visit.number = seq_along(Time)
        ) |>
        filter(!is.na(Outcome))

    followup.data <- model.data |>
        filter(Time > 0)

    intensity.model <-
        rlang::inject(coxph(
            Surv(prev_time,Time,!is.na(Outcome)) ~
                prev_outcome+strata(visit.number),
            id = Subject_ID,
            data = followup.data
        ))


    outcome.model <- PCORI_sim_outcome_modeler(
        Outcome ~
            ns(prev_outcome, df=3) +
            scale(Time) +
            scale(delta_time) - 1,
        data = followup.data)


    base <- SplineBasis(c(60,60,60,60,260,460,460,460,460))


    centering.statistics <-
        summarize(
            ungroup(filter(model.data, Time > 0, !is.na(Outcome))),
            across(c(Time, delta_time),
                   list(mean = ~ mean(.x, na.rm = TRUE),
                        sd = ~ sd(.x, na.rm = TRUE))
            )
        ) |>
        as.numeric() |>
        matrix(ncol = 2, byrow = TRUE) |>
        `dimnames<-`(list(c("time", "delta_time"), c("mean", "sd")))



    df_i <- model.data |>
        filter(Subject_ID==1)



    old.method <- compute_influence_for_one_alpha_and_one_patient(
        df_i,
        alpha = -0.6,
        variables = list(
            time = "Time",
            id = "Subject_ID",
            outcome = "Outcome",
            prev_time = "prev_time",
            prev_outcome = "prev_outcome",
            delta_time = "delta_time"
        ) |> map(rlang::sym),
        intensity.model = intensity.model,
        outcome.model = outcome.model,
        base = base,
        control = pcori_control(
            integration.method = 'quadv',
            tol = 1e-6
        ),
        centering = centering.statistics
    )
    # ----

    integration_data <-
        model.frame(
            terms(outcome.model),
            data = model.data |>
                dplyr::select(-prev_outcome) |>
                dplyr::mutate(
                    prev_outcome  = Outcome,
                    delta_time    = 0
                )
            ,
            id   = Subject_ID,
            time = Time
        )

    integration_X <- model.matrix( outcome.model, data=integration_data)

    new.method <-
        compute_influence_term_2_for_individual(
            times_ind = filter(integration_data, `(id)` == 1) |> pull(`(time)`),
            X_ind = integration_X[integration_data[["(id)"]] == 1,,drop=FALSE],
            X_all = model.matrix(outcome.model),
            Y_all = model.response(model.frame(outcome.model)),
            slope = c(0,0,0,1/centering.statistics[,'sd', drop=TRUE]),
            alpha = -0.6,
            beta = coef(outcome.model),
            base = base,
            bandwidth = outcome.model$bandwidth,
            tol = 1e-6
        )

    expect_true(
        all.equal(
            old.method[[1, 'term2']][[1]],
            new.method,
            check.attributes = FALSE,
            tolerance = 1e-6
        )
    )
})
