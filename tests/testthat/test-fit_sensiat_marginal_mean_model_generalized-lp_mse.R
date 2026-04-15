# Tests for fit_SensIAT_marginal_mean_model_generalized with lp_mse loss

test_that("fit_SensIAT_marginal_mean_model_generalized: lp_mse + identity produces valid results", {
    # Reduced parameters to keep test under 60 seconds
    setup <- generate_test_data(link = "identity", n_subjects = 10)
    
    # Fit using generalized version with identity link
    result_generalized <- suppressWarnings({
        fit_SensIAT_marginal_mean_model_generalized(
            data = setup$data,
            time = ..time..,
            id = ..id..,
            alpha = 0,
            knots = setup$knots,
            outcome.model = setup$outcome.model,
            intensity.model = setup$intensity.model,
            loss = "lp_mse",
            link = "identity",
            impute_data = create_impute_fn(),
            BBsolve.control = list(maxit = 50, tol = 1e-2),
            term2_method = "fast",
            use_expected_cache = FALSE
        )
    }, classes = "simpleWarning")
    
    # Verify structure and basic properties (not exact matching with original)
    expect_true(!is.null(result_generalized$coefficients),
                info = "Coefficients should not be NULL")
    expect_equal(length(result_generalized$coefficients), 1,
                 info = "Should have one coefficient vector for single alpha")
    expect_true(all(is.finite(result_generalized$coefficients[[1]])),
                info = "All coefficient values should be finite")
    
    expect_true(!is.null(result_generalized$coefficient.variance),
                info = "Coefficient variance should not be NULL")
    expect_equal(length(result_generalized$coefficient.variance), 1,
                 info = "Should have one variance matrix for single alpha")
    expect_true(all(is.finite(result_generalized$coefficient.variance[[1]])),
                info = "All variance values should be finite")
    
    expect_true(!is.null(result_generalized$influence),
                info = "Influence should not be NULL")
    expect_equal(length(result_generalized$influence), 1,
                 info = "Should have one influence tibble for single alpha")
})

test_that("fit_SensIAT_marginal_mean_model_generalized: lp_mse + log", {
    # Log link with term2 integration is computationally expensive
    skip_if_not(identical(Sys.getenv("RUN_SLOW_TESTS"), "true"), 
                "Log link test is slow. Set RUN_SLOW_TESTS=true to run")
    
    set.seed(123)
    setup <- generate_test_data(link = "log", n_subjects = 15)

    attr(setup$intensity.model, "bandwidth") <- 7

    expect_no_error({
        suppressWarnings({
            result <- fit_SensIAT_marginal_mean_model_generalized(
                data = setup$data,
                time = ..time..,
                id = ..id..,
                alpha = 0,
                knots = setup$knots,
                outcome.model = setup$outcome.model,
                intensity.model = setup$intensity.model,
                loss = "lp_mse",
                link = "log",
                impute_data = create_impute_fn(),
                BBsolve.control = list(maxit = 50, tol = 1e-2),
                term2_method = "fast"
            )
        }, classes = "simpleWarning")
    })
})

test_that("fit_SensIAT_marginal_mean_model_generalized: lp_mse + log (original term2)", {
    skip("Manual test")
    setup <- generate_test_data(link = "log", n_subjects = 15)
    
    expect_no_error({
        fit_SensIAT_marginal_mean_model_generalized(
            data = setup$data,
            time = ..time..,
            id = ..id..,
            alpha = 0,
            knots = setup$knots,
            outcome.model = setup$outcome.model,
            intensity.model = setup$intensity.model,
            loss = "lp_mse",
            link = "log",
            impute_data = create_impute_fn(),
            BBsolve.control = list(maxit = 100, tol = 1e-3),
            term2_method = "original"
        )
    })
})

test_that("fit_SensIAT_marginal_mean_model_generalized: lp_mse + logit", {
    # Logit link is computationally expensive due to numerical integration
    # Skip unless explicitly enabled via environment variable
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
                loss = "lp_mse",
                link = "logit",
                impute_data = create_impute_fn(),
                BBsolve.control = list(maxit = 30, tol = 1e-2),
                term2_method = "fast"
            )
        }, classes = "simpleWarning")
    })
})

test_that("fit_SensIAT_marginal_mean_model_generalized: lp_mse + identity with multiple alphas", {
    # Reduced parameters to keep test under 60 seconds
    setup <- generate_test_data(link = "identity", n_subjects = 8)
    
    # Test with multiple alpha values
    alphas <- c(-0.3, 0, 0.3)
    result_generalized <- suppressWarnings({
        fit_SensIAT_marginal_mean_model_generalized(
            data = setup$data,
            time = ..time..,
            id = ..id..,
            alpha = alphas,
            knots = setup$knots,
            outcome.model = setup$outcome.model,
            intensity.model = setup$intensity.model,
            loss = "lp_mse",
            link = "identity",
            impute_data = create_impute_fn(),
            BBsolve.control = list(maxit = 30, tol = 1e-2),
            term2_method = "fast"
        )
    })
    
    # Verify structure for multiple alphas
    expect_s3_class(result_generalized, "SensIAT_marginal_mean_model_generalized")
    expect_equal(length(result_generalized$coefficients), 3)
    expect_equal(result_generalized$alpha, alphas)
    
    # Each coefficient should be finite
    for (coef in result_generalized$coefficients) {
        expect_true(all(is.finite(coef)))
    }
})
