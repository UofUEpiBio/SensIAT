test_that("jackknife is invariant under parallelization.", {
    skip_on_cran()
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
    rlang::check_installed("future")
    rlang::check_installed("furrr")
    future::plan(future::sequential)
    jk1 <- jackknife(model, time = c(90, 180, 270, 360, 450))

    future::plan(future::multisession, workers=2)
    jk2 <- jackknife(model, time = c(90, 180, 270, 360, 450))

    expect_identical(jk1, jk2)
})
