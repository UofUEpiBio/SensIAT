run_MAVE_expectations <- function(mave.model){
    expect_length(mave.model$models$outcome$details$bw.delta, 1)
    expect_true(
        all(is.na(mave.model$models$outcome$details$beta.delta)) ||
        (length(mave.model$models$outcome$details$beta.delta)==5)
    )

    expect_true(mave.model$models$outcome$bandwidth > 0)
    expect_true(mave.model$models$outcome$coefficients[1] > 0)

    expect_equal(crossprod(mave.model$models$outcome$coefficients)[,], 1)

    if(  !is.null(mave.model$models$outcome$details$bw.last)
      && ! is.null(mave.model$models$outcome$details$beta.last)){
        if(
            mave.model$models$outcome$details$last.optimized == 'bandwidth'){
            expect_true(mave.model$models$outcome$details$bw.last$value <=
                            mave.model$models$outcome$details$beta.last$fval)
        } else {
            expect_true(mave.model$models$outcome$details$bw.last$value >=
                            mave.model$models$outcome$details$beta.last$fval)
        }
    }
}
test_that("MAVE: ise vs. mse", {
    testthat::skip_on_cran()
    object.mave.ise <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = fit_SensIAT_single_index_norm1coef_model,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60,260,460),
            End = 830,
            intensity.args=list(bandwidth=30),
            outcome.args=list(bw.selection='ise')
        )
    run_MAVE_expectations(object.mave.ise)
    object.mave.mse <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = fit_SensIAT_single_index_norm1coef_model,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60,260,460),
            End = 830,
            intensity.args=list(bandwidth=30),
            outcome.args=list(bw.selection='mse')
        )
    run_MAVE_expectations(object.mave.mse)

    expect_identical(
        coef(object.mave.ise$models$intensity),
        coef(object.mave.mse$models$intensity)
    )
})
test_that("MAVE: grid vs. optim", {
    testthat::skip_on_cran()
    object.mave.grid <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = fit_SensIAT_single_index_norm1coef_model,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60,260,460),
            End = 830,
            intensity.args=list(bandwidth=30),
            outcome.args=list(bw.method='grid', bw.selection = 'ise')
        )
    run_MAVE_expectations(object.mave.grid)

    object.mave.optim <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = fit_SensIAT_single_index_norm1coef_model,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60,260,460),
            End = 830,
            intensity.args=list(bandwidth=30),
            outcome.args=list(bw.method='optim', bw.selection = 'ise')
        )
    run_MAVE_expectations(object.mave.optim)

    object.mave.optimize <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = fit_SensIAT_single_index_norm1coef_model,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60,260,460),
            End = 830,
            intensity.args=list(bandwidth=30),
            outcome.args=list(bw.method='optimize', bw.selection = 'ise')
        )
    run_MAVE_expectations(object.mave.optimize)

    expect_equal(
        coef(object.mave.grid$models$outcome),
        coef(object.mave.optim$models$outcome)
    )
    expect_equal(
        coef(object.mave.grid$models$outcome),
        coef(object.mave.optimize$models$outcome)
    )
})
test_that("MAVE: reestimate.coef", {
    testthat::skip_on_cran()
    object.mave.wo <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = fit_SensIAT_single_index_norm1coef_model,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60,260,460),
            End = 830,
            intensity.args=list(bandwidth=30),
            outcome.args=list(bw.method='optimize', bw.selection = 'ise', reestimate.coef = 0)
        )
    run_MAVE_expectations(object.mave.wo)

    expect_true(is.na(object.mave.wo$models$outcome$details$beta.delta))
    expect_true(is.null(object.mave.wo$models$outcome$details$beta.last))

    object.mave.rc <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = fit_SensIAT_single_index_norm1coef_model,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60,260,460),
            End = 830,
            intensity.args=list(bandwidth=30),
            outcome.args=list(bw.method='optimize', bw.selection = 'ise', reestimate.coef = 1)
        )
    run_MAVE_expectations(object.mave.wo)

    object.mave.rc.grid <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = fit_SensIAT_single_index_norm1coef_model,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60,260,460),
            End = 830,
            intensity.args=list(bandwidth=30),
            outcome.args=list(bw.method='grid', bw.selection = 'ise', reestimate.coef = 1)
        )
    run_MAVE_expectations(object.mave.rc.grid)

    object.mave.rc.optim <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = fit_SensIAT_single_index_norm1coef_model,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60,260,460),
            End = 830,
            intensity.args=list(bandwidth=30),
            outcome.args=list(bw.method='optim', bw.selection = 'ise', reestimate.coef = 1)
        )
    run_MAVE_expectations(object.mave.rc.optim)

    object.mave.rc.optimize <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = fit_SensIAT_single_index_norm1coef_model,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60,260,460),
            End = 830,
            intensity.args=list(bandwidth=30),
            outcome.args=list(bw.method='optimize', bw.selection = 'ise', reestimate.coef = 1)
        )
    run_MAVE_expectations(object.mave.rc.optimize)

    expect_identical(object.mave.wo$models$intensity$coefficients,
                     object.mave.rc.optimize$models$intensity$coefficients)
    expect_identical(object.mave.wo$models$intensity$coefficients,
                     object.mave.rc.optim$models$intensity$coefficients)
    expect_identical(object.mave.wo$models$intensity$coefficients,
                     object.mave.rc.grid$models$intensity$coefficients)

})
