test_that("Altering models", {
    ruinf(1)
    old.seed <- .Random.seed
    model <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
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
    expect_identical(old.seed, .Random.seed)
    model$models$outcome |> formula() |>
    expect_equal(..outcome..~ns(..prev_outcome.., knots=c(9/6, 16/6)) + scale(..delta_time..)-1, ignore_attr=TRUE)
})
test_that("including terminal rows for intensity model", {
    model.no.terminals <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
            alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            End = 830,
            knots = c(60,260,460),
            add.terminal.observations = FALSE
        )
    model.with.terminals <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
            alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            End = 830,
            knots = c(60,260,460),
            add.terminal.observations = TRUE
        )
    model.external.terminals <-
        fit_SensIAT_within_group_model(
            group.data = add_terminal_observations(SensIAT_example_data, Subject_ID, Time, end = 830),
            outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
            alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            End = 830,
            knots = c(60,260,460),
            add.terminal.observations = FALSE,
        )
    expect_error(
        fit_SensIAT_within_group_model(
            group.data = add_terminal_observations(SensIAT_example_data, Subject_ID, Time, end = 830),
            outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
            alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            End = 830,
            knots = c(60,260,460),
            add.terminal.observations = TRUE,
        ),
        "Data contains missing values, cannot add terminal observations.")


    expect_false(
        coef(model.no.terminals$models$intensity) ==
        coef(model.with.terminals$models$intensity)
    )
    expect_equal(
        coef(model.with.terminals$models$intensity),
        coef(model.external.terminals$models$intensity)
    )
})
