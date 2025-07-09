test_that("Compute Influence term2 old vs. new methods", {
    model.data <- SensIAT_example_data |>
        dplyr::group_by(Subject_ID) |>
        dplyr::arrange(Time) |>
        dplyr::mutate(
            prev_time = dplyr::lag(Time),
            prev_outcome = dplyr::lag(Outcome),
            delta_time = Time - dplyr::lag(Time),
            visit.number = seq_along(Time)
        ) |>
        dplyr::filter(!is.na(Outcome))

    followup.data <- model.data |>
        filter(Time > 0)

    intensity.model <-
        rlang::inject(coxph(
            Surv(prev_time,Time,!is.na(Outcome)) ~
                prev_outcome+strata(visit.number),
            id = Subject_ID,
            data = followup.data
        ))


    outcome.model <- SensIAT_sim_outcome_modeler(
        Outcome ~
            ns(prev_outcome, df=3) +
            scale(Time) +
            scale(delta_time) - 1,
            id = Subject_ID,
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

    new.method <-
        compute_sim_influence_term_2_for_all_patients(
            outcome.model,
            integration_data = model.data |>
                filter(Subject_ID == 1) |>
                mutate(
                    prev_outcome = Outcome,
                    delta_time = 0
                ),
            alpha = -0.6,
            base = base,
            id = Subject_ID, time=Time,
            time.vars = c('delta_time')
        )


    expect_true(
        all.equal(
            old.method[[1, 'term2']][[1]],
            new.method[1,,drop=FALSE],
            check.attributes = FALSE,
            tolerance = 1e-6
        )
    )
})
