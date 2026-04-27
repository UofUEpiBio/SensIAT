# Tests for identity link specific behavior in fit_SensIAT_marginal_mean_model_generalized

test_that("fit_SensIAT_marginal_mean_model_generalized: identity link weight function", {
    # Verify that identity link produces constant weight function
    # (independent of beta, since ds/dz = 1)
    setup <- generate_test_data(link = "identity", n_subjects = 10)
    
    result1 <- suppressWarnings({
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
            BBsolve.control = list(maxit = 20, tol = 1e-2),
            term2_method = "original"
        )
    })
    
    result2 <- suppressWarnings({
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
            BBsolve.control = list(maxit = 20, tol = 1e-2),
            term2_method = "fast"
        )
    })
    
    # Both methods should produce similar results for identity link
    # since the weight function is simple
    expect_equal(result1$coefficients, result2$coefficients,
                 tolerance = 0.1,
                 info = "Original and fast term2 should be similar for identity link")
})

test_that("fit_SensIAT_marginal_mean_model_generalized: identity link convergence", {
    setup <- generate_test_data(link = "identity", n_subjects = 12)
    
    # Fit with different convergence tolerances
    result_loose <- suppressWarnings({
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
            BBsolve.control = list(maxit = 50, tol = 1e-1),
            term2_method = "fast"
        )
    })
    
    result_tight <- suppressWarnings({
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
            BBsolve.control = list(maxit = 50, tol = 1e-5),
            term2_method = "fast"
        )
    })
    
    # Tighter tolerance should give slightly different but converged result
        expect_true(all(is.finite(result_loose$coefficients[[1]])))
        expect_true(all(is.finite(result_tight$coefficients[[1]])))
    
    # They shouldn't be too different (at least same scale)
        expect_true(all(abs(result_loose$coefficients[[1]] - result_tight$coefficients[[1]]) < 
                        max(abs(result_tight$coefficients[[1]])) * 0.5))
})

test_that("fit_SensIAT_marginal_mean_model_generalized: identity link across loss functions", {
    # Compare lp_mse vs quasi-likelihood for identity link
    setup <- generate_test_data(link = "identity", n_subjects = 15)
    
    result_lp_mse <- suppressWarnings({
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
            BBsolve.control = list(maxit = 40, tol = 1e-3),
            term2_method = "fast"
        )
    })
    
    result_quasi <- suppressWarnings({
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
            BBsolve.control = list(maxit = 40, tol = 1e-3),
            term2_method = "fast"
        )
    })
    
    # Both should have finite coefficients
    expect_true(all(is.finite(result_lp_mse$coefficients[[1]])))
    expect_true(all(is.finite(result_quasi$coefficients[[1]])))
    
    # The two losses solve different estimating equations; do not require
    # numerical agreement, only valid outputs with matching dimensions.
    expect_equal(
        length(result_lp_mse$coefficients[[1]]),
        length(result_quasi$coefficients[[1]])
    )
})

test_that("fit_SensIAT_marginal_mean_model_generalized: identity single-index matches original", {
    setup <- generate_test_data(link = "identity", n_subjects = 10)

    result_original <- suppressWarnings({
        fit_SensIAT_marginal_mean_model(
            data = setup$data,
            id = ..id..,
            alpha = 0,
            knots = setup$knots,
            outcome.model = setup$outcome.model,
            intensity.model = setup$intensity.model,
            time.vars = c("..delta_time..")
        )
    })

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
            term2_method = "fast",
            time.vars = c("..delta_time..")
        )
    })

    expect_equal(
        result_generalized$coefficients[[1]],
        result_original$coefficients[[1]],
        tolerance = 1e-6,
        info = "Generalized identity+single-index should match original coefficients"
    )

    # Normalize influence representations: original stores matrix columns,
    # generalized stores per-patient vectors as list-columns.
    orig_ids <- result_original$influence[[1]]$id
    gen_ids <- result_generalized$influence[[1]]$id
    gen_order <- match(orig_ids, gen_ids)
    expect_false(anyNA(gen_order))

    orig_term1 <- result_original$influence[[1]]$term1
    if (is.list(orig_term1)) {
        orig_term1 <- do.call(rbind, orig_term1)
    }
    gen_term1 <- result_generalized$influence[[1]]$term1
    if (is.list(gen_term1)) {
        gen_term1 <- do.call(rbind, gen_term1)
    }
    gen_term1 <- gen_term1[gen_order, , drop = FALSE]

    orig_term2 <- result_original$influence[[1]]$term2
    if (is.list(orig_term2)) {
        orig_term2 <- do.call(rbind, orig_term2)
    }
    gen_term2 <- result_generalized$influence[[1]]$term2
    if (is.list(gen_term2)) {
        gen_term2 <- do.call(rbind, gen_term2)
    }
    gen_term2 <- gen_term2[gen_order, , drop = FALSE]

    expect_equal(
        gen_term1,
        orig_term1,
        tolerance = 1e-6,
        info = "Generalized identity+single-index should match original term1"
    )

    expect_equal(
        gen_term2,
        orig_term2,
        tolerance = 1e-6,
        info = "Generalized identity+single-index should match original term2"
    )
})
