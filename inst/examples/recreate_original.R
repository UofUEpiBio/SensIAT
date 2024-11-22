fixed.object <-
    fit_SensIAT_within_group_model(
        group.data = SensIAT_example_data,
        outcome_modeler = SensIAT_sim_outcome_modeler,
        alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        knots = c(60,260,460),
        End = 830,
        intensity.args = list(bandwidth = 30),
        influence.args = list(method = 'fixed', delta = 1)
    )
