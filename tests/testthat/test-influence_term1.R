test_that("Compute Influence term1 old vs. new methods", {
    model.data <- SensIAT_example_data |>
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
        control = pcori_control(integration.method = 'quadv'),
        centering = centering.statistics
    )
    # ----

    baseline_intensity_all =
        estimate_baseline_intensity(
            intensity.model = intensity.model,
            data = followup.data,
            variables = list(
                time = "Time",
                id = "Subject_ID",
                outcome = "Outcome",
                prev_time = "prev_time",
                prev_outcome = "prev_outcome"
            ) |> map(rlang::sym),
            bandwidth = NULL
        )

    new.method <-
        compute_influence_term_1_for_all(
            X_all = model.matrix(outcome.model),
            times_all = pull(followup.data, Time),
            outcome_all = pull(followup.data, Outcome),
            prev_outcome_all = pull(followup.data, prev_outcome),
            baseline_intensity_all = baseline_intensity_all$baseline_intensity,
            alpha = -0.6,
            intensity_coef = coef(intensity.model),
            outcome_coef = coef(outcome.model),
            base = base,
            bandwidth = outcome.model$bandwidth,
            kernel = attr(outcome.model, 'kernel')
        )


    expect_equal(
        old.method[[1, 'term1']][[1]],
        new.method[pull(followup.data, Subject_ID) == 1,]
    )
})
