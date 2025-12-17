# Helper function to prepare test data and models
prepare_test_data_and_models <- function(data = SensIAT_example_data) {
    data_with_lags <- data |>
        dplyr::group_by(Subject_ID) |>
        dplyr::mutate(
            ..prev_outcome.. = dplyr::lag(Outcome, default = NA_real_, order_by = Time),
            ..prev_time.. = dplyr::lag(Time, default = 0, order_by = Time),
            ..delta_time.. = Time - dplyr::lag(.data$Time, default = NA_real_, order_by = Time)
        )
    
    # Create the observation time intensity model
    intensity.model <-
        survival::coxph(survival::Surv(..prev_time.., Time, !is.na(Outcome)) ~ ..prev_outcome.. + strata(Visit),
            data = data_with_lags |> dplyr::filter(.data$Time > 0)
        )
    
    # Create the observed outcome model
    outcome.model <-
        fit_SensIAT_single_index_fixed_coef_model(
            Outcome ~ splines::ns(..prev_outcome.., df = 3) + ..delta_time.. - 1,
            id = Subject_ID,
            data = data_with_lags |> dplyr::filter(Time > 0)
        )
    
    list(
        data = data_with_lags,
        intensity.model = intensity.model,
        outcome.model = outcome.model
    )
}

# Helper function to create impute_data function
create_impute_fn <- function() {
    function(t, df) {
        data_wl <- df |>
            dplyr::mutate(
                ..prev_time.. = Time,
                ..prev_outcome.. = Outcome,
                ..delta_time.. = 0
            )
        extrapolate_from_last_observation(t, data_wl, "Time", slopes = c("..delta_time.." = 1))
    }
}

# Test all combinations of loss and link functions
# Note: Currently only lp_mse + log is fully implemented
# Other combinations have incomplete implementations (missing W function, etc.)

test_that("fit_SensIAT_marginal_mean_model_generalized: lp_mse + identity", {
    skip("Identity link delegates to fit_SensIAT_marginal_mean_model with API mismatch")
    setup <- prepare_test_data_and_models()
    
    expect_no_error({
        fit_SensIAT_marginal_mean_model_generalized(
            data = setup$data,
            time = setup$data$Time,
            id = setup$data$Subject_ID,
            alpha = 0,
            knots = c(60, 260, 460),
            outcome.model = setup$outcome.model,
            intensity.model = setup$intensity.model,
            loss = "lp_mse",
            link = "identity",
            impute_data = create_impute_fn(),
            BBsolve.control = list(maxit = 10, tol = 1e-4),
            term2_method = "fast"
        )
    })
})

test_that("fit_SensIAT_marginal_mean_model_generalized: lp_mse + log", {
    setup <- prepare_test_data_and_models()
    
    expect_no_error({
        fit_SensIAT_marginal_mean_model_generalized(
            data = setup$data,
            time = setup$data$Time,
            id = setup$data$Subject_ID,
            alpha = 0,
            knots = c(60, 260, 460),
            outcome.model = setup$outcome.model,
            intensity.model = setup$intensity.model,
            loss = "lp_mse",
            link = "log",
            impute_data = create_impute_fn(),
            BBsolve.control = list(maxit = 10, tol = 1e-4),
            term2_method = "fast"
        )
    })
})

test_that("fit_SensIAT_marginal_mean_model_generalized: lp_mse + log (original term2)", {
    skip("Manual test")
    setup <- prepare_test_data_and_models()
    
    expect_no_error({
        fit_SensIAT_marginal_mean_model_generalized(
            data = setup$data,
            time = setup$data$Time,
            id = setup$data$Subject_ID,
            alpha = 0,
            knots = c(60, 260, 460),
            outcome.model = setup$outcome.model,
            intensity.model = setup$intensity.model,
            loss = "lp_mse",
            link = "log",
            impute_data = create_impute_fn(),
            BBsolve.control = list(maxit = 10, tol = 1e-4),
            term2_method = "original"
        )
    })
})

test_that("fit_SensIAT_marginal_mean_model_generalized: lp_mse + logit", {
    setup <- prepare_test_data_and_models()
    
    expect_no_error({
        fit_SensIAT_marginal_mean_model_generalized(
            data = setup$data,
            time = setup$data$Time,
            id = setup$data$Subject_ID,
            alpha = 0,
            knots = c(60, 260, 460),
            outcome.model = setup$outcome.model,
            intensity.model = setup$intensity.model,
            loss = "lp_mse",
            link = "logit",
            impute_data = create_impute_fn(),
            BBsolve.control = list(maxit = 10, tol = 1e-4),
            term2_method = "fast"
        )
    })
})

test_that("fit_SensIAT_marginal_mean_model_generalized: quasi-likelihood + identity", {
    skip("Quasi-likelihood loss not fully implemented (delegates to identity which has API mismatch)")
    setup <- prepare_test_data_and_models()
    
    expect_no_error({
        fit_SensIAT_marginal_mean_model_generalized(
            data = setup$data,
            time = setup$data$Time,
            id = setup$data$Subject_ID,
            alpha = 0,
            knots = c(60, 260, 460),
            outcome.model = setup$outcome.model,
            intensity.model = setup$intensity.model,
            loss = "quasi-likelihood",
            link = "identity",
            impute_data = create_impute_fn(),
            BBsolve.control = list(maxit = 10, tol = 1e-4),
            term2_method = "fast"
        )
    })
})

test_that("fit_SensIAT_marginal_mean_model_generalized: quasi-likelihood + log", {
    skip("Quasi-likelihood loss not fully implemented (missing W function definition)")
    setup <- prepare_test_data_and_models()
    
    expect_no_error({
        fit_SensIAT_marginal_mean_model_generalized(
            data = setup$data,
            time = setup$data$Time,
            id = setup$data$Subject_ID,
            alpha = 0,
            knots = c(60, 260, 460),
            outcome.model = setup$outcome.model,
            intensity.model = setup$intensity.model,
            loss = "quasi-likelihood",
            link = "log",
            impute_data = create_impute_fn(),
            BBsolve.control = list(maxit = 10, tol = 1e-4),
            term2_method = "fast"
        )
    })
})

test_that("fit_SensIAT_marginal_mean_model_generalized: quasi-likelihood + logit", {
    skip("Quasi-likelihood loss not fully implemented (missing W function definition)")
    setup <- prepare_test_data_and_models()
    
    expect_no_error({
        fit_SensIAT_marginal_mean_model_generalized(
            data = setup$data,
            time = setup$data$Time,
            id = setup$data$Subject_ID,
            alpha = 0,
            knots = c(60, 260, 460),
            outcome.model = setup$outcome.model,
            intensity.model = setup$intensity.model,
            loss = "quasi-likelihood",
            link = "logit",
            impute_data = create_impute_fn(),
            BBsolve.control = list(maxit = 10, tol = 1e-4),
            term2_method = "fast"
        )
    })
})
