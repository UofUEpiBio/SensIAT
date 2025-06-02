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
