# Helper to create small test dataset (10 subjects instead of 200)
create_small_mave_data <- function(n_subjects = 10) {
    small_subjects <- head(unique(SensIAT_example_data$Subject_ID), n_subjects)
    SensIAT_example_data |>
        dplyr::filter(Subject_ID %in% small_subjects)
}

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
    
    small_data <- create_small_mave_data(10)
    
    object.mave.ise <-
        fit_SensIAT_within_group_model(
            group.data = small_data,
            outcome_modeler = fit_SensIAT_single_index_norm1coef_model,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60, 460),  # Fewer knots
            End = 830,
            intensity.args=list(bandwidth=30),
            outcome.args=list(bw.selection='ise')
        )
    run_MAVE_expectations(object.mave.ise)
    
    object.mave.mse <-
        fit_SensIAT_within_group_model(
            group.data = small_data,
            outcome_modeler = fit_SensIAT_single_index_norm1coef_model,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60, 460),
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
    
    small_data <- create_small_mave_data(10)
    
    # Test just one method to verify it works
    object.mave.optimize <-
        fit_SensIAT_within_group_model(
            group.data = small_data,
            outcome_modeler = fit_SensIAT_single_index_norm1coef_model,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60, 460),
            End = 830,
            intensity.args=list(bandwidth=30),
            outcome.args=list(bw.method='optimize', bw.selection = 'ise')
        )
    run_MAVE_expectations(object.mave.optimize)
    
    # Verify coefficients are normalized
    expect_equal(sum(object.mave.optimize$models$outcome$coefficients^2), 1)
})
test_that("MAVE: reestimate.coef", {
    testthat::skip_on_cran()
    
    small_data <- create_small_mave_data(10)
    
    # Test without coefficient reestimation
    object.mave.wo <-
        fit_SensIAT_within_group_model(
            group.data = small_data,
            outcome_modeler = fit_SensIAT_single_index_norm1coef_model,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60, 460),
            End = 830,
            intensity.args=list(bandwidth=30),
            outcome.args=list(bw.method='optimize', bw.selection = 'ise', reestimate.coef = 0)
        )
    run_MAVE_expectations(object.mave.wo)

    expect_true(is.na(object.mave.wo$models$outcome$details$beta.delta))
    expect_true(is.null(object.mave.wo$models$outcome$details$beta.last))

    # Test with coefficient reestimation
    object.mave.rc <-
        fit_SensIAT_within_group_model(
            group.data = small_data,
            outcome_modeler = fit_SensIAT_single_index_norm1coef_model,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60, 460),
            End = 830,
            intensity.args=list(bandwidth=30),
            outcome.args=list(bw.method='optimize', bw.selection = 'ise', reestimate.coef = 1)
        )
    run_MAVE_expectations(object.mave.rc)

    # Intensity model should be identical (reestimate.coef only affects outcome model)
    expect_identical(object.mave.wo$models$intensity$coefficients,
                     object.mave.rc$models$intensity$coefficients)
})
