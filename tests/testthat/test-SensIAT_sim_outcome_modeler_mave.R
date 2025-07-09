test_that("MAVE: ise vs. mse", {
    object.mave.ise <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = SensIAT_sim_outcome_modeler_mave,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60,260,460),
            End = 830,
            intensity.args=list(bandwidth=30),
            outcome.args=list(bw.selection='ise')
        )
    object.mave.mse <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = SensIAT_sim_outcome_modeler_mave,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60,260,460),
            End = 830,
            intensity.args=list(bandwidth=30),
            outcome.args=list(bw.selection='mse')
        )

    object.mave.ise$models$outcome$bandwidth
    object.mave.mse$models$outcome$bandwidth

    expect_identical(
        coef(object.mave.ise$models$intensity),
        coef(object.mave.mse$models$intensity)
    )
})
test_that("MAVE: grid vs. optim", {
    object.mave.grid <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = SensIAT_sim_outcome_modeler_mave,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60,260,460),
            End = 830,
            intensity.args=list(bandwidth=30),
            outcome.args=list(bw.method='grid', bw.selection = 'ise')
        )
    object.mave.optim <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = SensIAT_sim_outcome_modeler_mave,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60,260,460),
            End = 830,
            intensity.args=list(bandwidth=30),
            outcome.args=list(bw.method='optim', bw.selection = 'ise')
        )
    object.mave.optimize <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = SensIAT_sim_outcome_modeler_mave,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60,260,460),
            End = 830,
            intensity.args=list(bandwidth=30),
            outcome.args=list(bw.method='optimize', bw.selection = 'ise')
        )

    object.mave.grid$models$outcome$bandwidth
    object.mave.grid$models$outcome$details$value
    object.mave.optim$models$outcome$bandwidth
    object.mave.optim$models$outcome$details$value
    object.mave.optimize$models$outcome$bandwidth
    object.mave.optimize$models$outcome$details$value


    expect_equal(
        coef(object.mave.grid$models$outcome),
        coef(object.mave.optim$models$outcome)
    )
    expect_equal(
        coef(object.mave.grid$models$outcome),
        coef(object.mave.optimize$models$outcome)
    )

    if(interactive()){
    object.mave.grid$models$outcome$details |>
        with(tibble(log_bw_seq, err_values)) |>
        mutate(bw = exp(log_bw_seq)) |>
        ggplot(aes(x=bw, y=err_values)) %+%
        geom_line() %+%
        geom_point() +
        geom_vline(
            data = tibble(
                'Selected' = c(object.mave.grid$models$outcome$bandwidth,
                               object.mave.optim$models$outcome$bandwidth,
                               exp(object.mave.optim$models$outcome$details$initial)
                               ),
                'Method' = c('grid', 'optim', 'initial')
            ) |> mutate(Method = paste(Method, " (", format(Selected, digits = 4), ")")),
            aes(xintercept = Selected, color = Method)
        )
    }
})
test_that("MAVE: reestimate.coef", {
    object.mave.wo <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = SensIAT_sim_outcome_modeler_mave,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60,260,460),
            End = 830,
            intensity.args=list(bandwidth=30),
            outcome.args=list(bw.method='grid', bw.selection = 'ise', reestimate.coef = FALSE)
        )
    object.mave.rc <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = SensIAT_sim_outcome_modeler_mave,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60,260,460),
            End = 830,
            intensity.args=list(bandwidth=30),
            outcome.args=list(bw.method='grid', bw.selection = 'ise', reestimate.coef = TRUE)
        )

    expect_identical(object.mave.wo$models$intensity$coefficients,
                     object.mave.rc$models$intensity$coefficients)
    expect_identical(object.mave.wo$models$outcome$bandwidth,
                     object.mave.rc$models$outcome$bandwidth)
    expect_true(object.mave.wo$models$outcome$details$value >= object.mave.rc$models$outcome$details.refit$fval)
})
