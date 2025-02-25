test_that("Altering models", {
    model <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = SensIAT_sim_outcome_modeler,
            alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            End = 830,
            knots = c(60,260,460),
            outcome.args = list(
                model=~ns(..prev_outcome.., knots=c(9/6, 16/6)) + scale(..delta_time..)-1
            )
        )
    model$models$outcome |> formula() |>
    expect_equal(..outcome..~ns(..prev_outcome.., knots=c(9/6, 16/6)) + scale(..delta_time..)-1, ignore_attr=TRUE)
})
