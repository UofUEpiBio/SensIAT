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

    object.mave.grid$models$outcome$bandwidth
    object.mave.optim$models$outcome$bandwidth

    if(!interactive()){
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
                               exp(object.mave.grid$models$outcome$details$log_bw_seq[1])/0.05
                               ),
                'Method' = c('grid', 'optim', 'initial')
            ) |> mutate(Method = paste(Method, " (", format(Selected, digits = 4), ")")),
            aes(xintercept = Selected, color = Method)
        )
    }
})
