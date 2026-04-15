# Tests for fit_SensIAT_marginal_mean_model_generalized with quasi-likelihood loss

test_that("fit_SensIAT_marginal_mean_model_generalized: quasi-likelihood + identity (now integrated)", {
    # Reduced parameters to keep test under 60 seconds
    setup <- generate_test_data(link = "identity", n_subjects = 10)
    
    # Identity link should now work with unified infrastructure
    result_generalized <- suppressWarnings({
        fit_SensIAT_marginal_mean_model_generalized(
            data = setup$data,
            time = ..time..,
            id = ..id..,
            alpha = 0,
            knots = setup$knots,
            outcome.model = setup$outcome.model,
            intensity.model = setup$intensity.model,
            loss = "quasi-likelihood",
            link = "identity",
            impute_data = create_impute_fn(),
            BBsolve.control = list(maxit = 30, tol = 1e-2),
            term2_method = "fast"
        )
    })
    
    # Verify structure and basic properties
    expect_s3_class(result_generalized, "SensIAT_marginal_mean_model_generalized")
    expect_true(length(result_generalized$coefficients) > 0)
    expect_equal(length(result_generalized$coefficients), 1,
                 info = "Single alpha should return list with 1 element")
    expect_true(all(is.finite(result_generalized$coefficients[[1]])))
    expect_equal(result_generalized$alpha, 0)
})

test_that("fit_SensIAT_marginal_mean_model_generalized: quasi-likelihood + log", {
    # Log link with term2 integration is computationally expensive
    skip_if_not(identical(Sys.getenv("RUN_SLOW_TESTS"), "true"), 
                "Log link test is slow. Set RUN_SLOW_TESTS=true to run")
    
    setup <- generate_test_data(link = "log", n_subjects = 10)
    
    expect_no_error({
        suppressWarnings({
            fit_SensIAT_marginal_mean_model_generalized(
                data = setup$data,
                time = ..time..,
                id = ..id..,
                alpha = 0,
                knots = setup$knots,
                outcome.model = setup$outcome.model,
                intensity.model = setup$intensity.model,
                loss = "quasi-likelihood",
                link = "log",
                impute_data = create_impute_fn(),
                BBsolve.control = list(maxit = 30, tol = 1e-2),
                term2_method = "fast"
            )
        }, classes = "simpleWarning")
    })
})

test_that("fit_SensIAT_marginal_mean_model_generalized: quasi-likelihood + logit", {
    # Logit link is computationally expensive - skip unless explicitly enabled
    skip_if_not(identical(Sys.getenv("RUN_SLOW_TESTS"), "true"), 
                "Logit link test is slow. Set RUN_SLOW_TESTS=true to run")
    
    setup <- generate_test_data(link = "logit", n_subjects = 10)
    
    expect_no_error({
        suppressWarnings({
            fit_SensIAT_marginal_mean_model_generalized(
                data = setup$data,
                time = ..time..,
                id = ..id..,
                alpha = 0,
                knots = setup$knots,
                outcome.model = setup$outcome.model,
                intensity.model = setup$intensity.model,
                loss = "quasi-likelihood",
                link = "logit",
                impute_data = create_impute_fn(),
                BBsolve.control = list(maxit = 30, tol = 1e-2),
                term2_method = "fast"
            )
        }, classes = "simpleWarning")
    })
})
