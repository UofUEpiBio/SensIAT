create_within_group_test_fixture <- function() {
    ids <- head(unique(SensIAT_example_data$Subject_ID), 8)
    dplyr::filter(SensIAT_example_data, Subject_ID %in% ids)
}

test_that("Altering models", {
    runif(1)
    old.seed <- .Random.seed
    fixture_data <- create_within_group_test_fixture()
    model <-
        fit_SensIAT_within_group_model(
            group.data = fixture_data,
            outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
            alpha = 0,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            End = 830,
            knots = c(60, 260, 460),
            outcome.args = list(
                model = ~ ns(..prev_outcome.., knots = c(9 / 6, 16 / 6)) + scale(..delta_time..) - 1
            )
        )
    expect_identical(old.seed, .Random.seed)
    model$models$outcome |>
        formula() |>
        expect_equal(..outcome.. ~ ns(..prev_outcome.., knots = c(9 / 6, 16 / 6)) + scale(..delta_time..) - 1, ignore_attr = TRUE)
})
test_that("including terminal rows for intensity model", {
    fixture_data <- create_within_group_test_fixture()
    model.no.terminals <-
        fit_SensIAT_within_group_model(
            group.data = fixture_data,
            outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
            alpha = 0,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            End = 830,
            knots = c(60, 260, 460),
            add.terminal.observations = FALSE
        )
    model.with.terminals <-
        fit_SensIAT_within_group_model(
            group.data = fixture_data,
            outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
            alpha = 0,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            End = 830,
            knots = c(60, 260, 460),
            add.terminal.observations = TRUE
        )
    model.external.terminals <-
        fit_SensIAT_within_group_model(
            group.data = add_terminal_observations(fixture_data, Subject_ID, Time, end = 830),
            outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
            alpha = 0,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            End = 830,
            knots = c(60, 260, 460),
            add.terminal.observations = FALSE,
        )
    expect_error(
        fit_SensIAT_within_group_model(
            group.data = add_terminal_observations(fixture_data, Subject_ID, Time, end = 830),
            outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
            alpha = 0,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            End = 830,
            knots = c(60, 260, 460),
            add.terminal.observations = TRUE,
        ),
        "Data contains missing values, cannot add terminal observations."
    )


    expect_false(
        coef(model.no.terminals$models$intensity) ==
            coef(model.with.terminals$models$intensity)
    )
    expect_equal(
        coef(model.with.terminals$models$intensity),
        coef(model.external.terminals$models$intensity)
    )
})
