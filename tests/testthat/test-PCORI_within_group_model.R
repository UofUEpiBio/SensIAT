test_that("Altering models", {
    runif(1)
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
    model.no.terminals <-
        fit_SensIAT_within_group_model(
            group.data = SensIAT_example_data,
            outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
            alpha = c(-0.6, -0.3, 0, 0.3, 0.6),
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            End = 830,
            knots = c(60, 260, 460),
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
            knots = c(60, 260, 460),
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
            knots = c(60, 260, 460),
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

# ============================================================================
# BACKWARD COMPATIBILITY TESTS
# These tests ensure existing code continues to work after adding generalized support
# ============================================================================

test_that("default parameters maintain backward compatibility", {
    # Test 1: Default call without any new parameters
    model_default <- fit_SensIAT_within_group_model(
        group.data = SensIAT_example_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        alpha = 0,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        End = 830,
        knots = c(60, 260, 460)
    )

    # Should have default link/loss/term2_method in output
    expect_equal(model_default$link, "identity")
    expect_equal(model_default$loss, "lp_mse")
    expect_equal(model_default$term2_method, "fast")

    # Should have all standard output components
    expect_s3_class(model_default, "SensIAT_within_group_model")
    expect_true(!is.null(model_default$coefficients))
    expect_true(!is.null(model_default$coefficient.variance))
    expect_true(!is.null(model_default$influence))
    expect_true(!is.null(model_default$base))
    expect_true(!is.null(model_default$V_inverse))
})

test_that("explicit identity link produces same results as default", {
    # Model with default (no link specified)
    model_implicit <- fit_SensIAT_within_group_model(
        group.data = SensIAT_example_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        alpha = 0,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        End = 830,
        knots = c(60, 260, 460)
    )

    # Model with explicit link = "identity"
    model_explicit <- fit_SensIAT_within_group_model(
        group.data = SensIAT_example_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        alpha = 0,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        End = 830,
        knots = c(60, 260, 460),
        link = "identity"
    )

    # Coefficients should be identical
    expect_equal(model_implicit$coefficients, model_explicit$coefficients)
    expect_equal(model_implicit$coefficient.variance, model_explicit$coefficient.variance)
})

test_that("multiple alpha values work with default identity link", {
    model <- fit_SensIAT_within_group_model(
        group.data = SensIAT_example_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        alpha = c(-0.3, 0, 0.3),
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        End = 830,
        knots = c(60, 260, 460)
    )

    # Should have results for each alpha
    expect_length(model$coefficients, 3)
    expect_length(model$coefficient.variance, 3)
    expect_equal(model$alpha, c(-0.3, 0, 0.3))
})

# ============================================================================
# NEW PATTERN TESTS
# These tests verify the new generalized link/loss functionality
# ============================================================================

test_that("log link with lp_mse loss runs without error", {
    skip_on_cran()

    model_log <- fit_SensIAT_within_group_model(
        group.data = SensIAT_example_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        alpha = 0,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        End = 830,
        knots = c(60, 260, 460),
        link = "log",
        loss = "lp_mse"
    )

    expect_s3_class(model_log, "SensIAT_within_group_model")
    expect_equal(model_log$link, "log")
    expect_equal(model_log$loss, "lp_mse")
    expect_true(!is.null(model_log$coefficients))
    expect_true(!is.null(model_log$coefficient.variance))
})

test_that("logit link with lp_mse loss runs without error", {
    skip_on_cran()

    model_logit <- fit_SensIAT_within_group_model(
        group.data = SensIAT_example_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        alpha = 0,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        End = 830,
        knots = c(60, 260, 460),
        link = "logit",
        loss = "lp_mse"
    )

    expect_s3_class(model_logit, "SensIAT_within_group_model")
    expect_equal(model_logit$link, "logit")
    expect_equal(model_logit$loss, "lp_mse")
    expect_true(!is.null(model_logit$coefficients))
})

test_that("quasi-likelihood loss with log link runs without error", {
    skip_on_cran()

    model_ql <- fit_SensIAT_within_group_model(
        group.data = SensIAT_example_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        alpha = 0,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        End = 830,
        knots = c(60, 260, 460),
        link = "log",
        loss = "quasi-likelihood"
    )

    expect_s3_class(model_ql, "SensIAT_within_group_model")
    expect_equal(model_ql$link, "log")
    expect_equal(model_ql$loss, "quasi-likelihood")
})

test_that("term2_method parameter is respected for generalized models", {
    skip_on_cran()

    # Test with fixed_grid method
    model_fixed <- fit_SensIAT_within_group_model(
        group.data = SensIAT_example_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        alpha = 0,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        End = 830,
        knots = c(60, 260, 460),
        link = "log",
        term2_method = "fixed_grid"
    )

    expect_equal(model_fixed$term2_method, "fixed_grid")
    expect_true(!is.null(model_fixed$coefficients))
})

test_that("custom impute_data function is used when provided", {
    skip_on_cran()

    # Track whether custom impute_data was called
    impute_called <- FALSE

    custom_impute <- function(t, df) {
        impute_called <<- TRUE
        data_wl <- df |>
            dplyr::mutate(
                ..prev_time.. = .data$..time..,
                ..prev_outcome.. = .data$..outcome..,
                ..delta_time.. = 0
            )
        extrapolate_from_last_observation(t, data_wl, "..time..", slopes = c("..delta_time.." = 1))
    }

    model_custom <- fit_SensIAT_within_group_model(
        group.data = SensIAT_example_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        alpha = 0,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        End = 830,
        knots = c(60, 260, 460),
        link = "log",
        impute_data = custom_impute
    )

    expect_true(impute_called)
    expect_s3_class(model_custom, "SensIAT_within_group_model")
})

test_that("multiple alpha values work with log link", {
    skip_on_cran()

    model_multi <- fit_SensIAT_within_group_model(
        group.data = SensIAT_example_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        alpha = c(-0.3, 0, 0.3),
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        End = 830,
        knots = c(60, 260, 460),
        link = "log"
    )

    expect_length(model_multi$coefficients, 3)
    expect_length(model_multi$coefficient.variance, 3)
    expect_equal(model_multi$alpha, c(-0.3, 0, 0.3))
})
