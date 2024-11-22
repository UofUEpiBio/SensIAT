# A basic example using fixed intensity bandwidth.
object <-
    fit_SensIAT_within_group_model(
        group.data = SensIAT_example_data,
        outcome_modeler = SensIAT_sim_outcome_modeler,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        knots = c(60,260,460),
        End = 830,
        intensity.args=list(bandwidth=30)
    )
