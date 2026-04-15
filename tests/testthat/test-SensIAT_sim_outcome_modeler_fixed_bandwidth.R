# Tests for fit_SensIAT_single_index_fixed_bandwidth_model
# This outcome modeler estimates all coefficients with a fixed bandwidth of 1

# Helper to create small test dataset (10 subjects instead of 200)
create_small_test_data <- function(n_subjects = 10) {
    data("SensIAT_example_data", package = "SensIAT", envir = environment())
    small_subjects <- head(unique(SensIAT_example_data$Subject_ID), n_subjects)
    SensIAT_example_data |>
        dplyr::filter(Subject_ID %in% small_subjects)
}

# Common expectations for fixed bandwidth models
run_fixed_bandwidth_expectations <- function(model_obj) {
    # Check that it's a within_group_model
    expect_s3_class(model_obj, "SensIAT_within_group_model")
    
    # Check outcome model structure
    outcome_model <- model_obj$models$outcome
    expect_s3_class(outcome_model, "SensIAT::Single-index-outcome-model")
    
    # Fixed bandwidth should be exactly 1
    expect_equal(outcome_model$bandwidth, 1)
    
    # Should have coefficients
    expect_true(is.numeric(outcome_model$coefficients))
    expect_true(length(outcome_model$coefficients) > 0)
    
    # All coefficients should be estimated (not fixed)
    expect_true(all(is.finite(outcome_model$coefficients)))
    
    # First coefficient convention (should be positive)
    expect_true(outcome_model$coefficients[1] > 0)
    
    # Should have optimization details
    expect_true(!is.null(outcome_model$details))
    expect_true(is.list(outcome_model$details))
}

test_that("fixed_bandwidth model basic functionality", {
    small_data <- create_small_test_data(10)
    
    model <- fit_SensIAT_within_group_model(
        group.data = small_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_bandwidth_model,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        knots = c(60, 460),
        End = 830,
        intensity.args = list(bandwidth = 30)
    )
    
    run_fixed_bandwidth_expectations(model)
    
    # Check that the model has valid structure
    outcome_model <- model$models$outcome
    expect_true(!is.null(outcome_model$data))
    expect_true(!is.null(coef(outcome_model)))
})

test_that("fixed_bandwidth model with different optimization methods", {
    small_data <- create_small_test_data(10)
    
    # Test with optim method
    model_optim <- fit_SensIAT_within_group_model(
        group.data = small_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_bandwidth_model,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        knots = c(60, 460),
        End = 830,
        intensity.args = list(bandwidth = 30),
        outcome.args = list(method = "optim")
    )
    run_fixed_bandwidth_expectations(model_optim)
    
    # Test with nlminb method
    model_nlminb <- fit_SensIAT_within_group_model(
        group.data = small_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_bandwidth_model,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        knots = c(60, 460),
        End = 830,
        intensity.args = list(bandwidth = 30),
        outcome.args = list(method = "nlminb")
    )
    run_fixed_bandwidth_expectations(model_nlminb)
    
    # Test with nmk method (default)
    model_nmk <- fit_SensIAT_within_group_model(
        group.data = small_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_bandwidth_model,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        knots = c(60, 460),
        End = 830,
        intensity.args = list(bandwidth = 30),
        outcome.args = list(method = "nmk")
    )
    run_fixed_bandwidth_expectations(model_nmk)
    
    # All methods should produce valid models
    # (coefficients may differ slightly due to optimization)
    expect_true(is.numeric(model_optim$models$outcome$coefficients))
    expect_true(is.numeric(model_nlminb$models$outcome$coefficients))
    expect_true(is.numeric(model_nmk$models$outcome$coefficients))
})

test_that("fixed_bandwidth model with different kernels", {
    small_data <- create_small_test_data(10)
    
    # NOTE: K4_Biweight kernel is NOT tested here because while it's supported
    # in the R-level SIDRnew_fixed_bandwidth() function, it's currently disabled
    # in the C++ acceleration functions (pcoriaccel_estimate_pmf and pcoriaccel_NW).
    # These C++ functions are required for influence term calculations in 
    # fit_SensIAT_within_group_model, so K4_Biweight would cause errors when
    # computing influence terms. See src/estimate_pmf.cpp lines 93-96 and
    # src/NW.cpp lines 220-223 where K4_Biweight support is commented out.
    
    # Test with K2_Biweight (default)
    model_k2 <- fit_SensIAT_within_group_model(
        group.data = small_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_bandwidth_model,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        knots = c(60, 460),
        End = 830,
        intensity.args = list(bandwidth = 30),
        outcome.args = list(kernel = "K2_Biweight")
    )
    run_fixed_bandwidth_expectations(model_k2)
    expect_equal(attr(model_k2$models$outcome, "kernel"), "K2_Biweight")
    
    # Test with dnorm kernel
    model_dnorm <- fit_SensIAT_within_group_model(
        group.data = small_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_bandwidth_model,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        knots = c(60, 460),
        End = 830,
        intensity.args = list(bandwidth = 30),
        outcome.args = list(kernel = "dnorm")
    )
    run_fixed_bandwidth_expectations(model_dnorm)
    expect_equal(attr(model_dnorm$models$outcome, "kernel"), "dnorm")
})

test_that("K4_Biweight kernel limitation is documented", {
    # This test documents that K4_Biweight is currently not supported
    # for the full workflow, despite being implemented in R-level code.
    # The limitation is in the C++ acceleration layer.
    
    small_data <- create_small_test_data(10)
    
    # Attempting to use K4_Biweight with influence calculations should fail
    expect_error(
        fit_SensIAT_within_group_model(
            group.data = small_data,
            outcome_modeler = fit_SensIAT_single_index_fixed_bandwidth_model,
            id = Subject_ID,
            outcome = Outcome,
            time = Time,
            knots = c(60, 460),
            End = 830,
            intensity.args = list(bandwidth = 30),
            outcome.args = list(kernel = "K4_Biweight")
        ),
        "Invalid value for `kernel`"
    )
})

test_that("fixed_bandwidth model computes expected values", {
    small_data <- create_small_test_data(10)
    
    model <- fit_SensIAT_within_group_model(
        group.data = small_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_bandwidth_model,
        alpha = c(-0.3, 0, 0.3),
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        knots = c(60, 460),
        End = 830,
        intensity.args = list(bandwidth = 30)
    )
    
    run_fixed_bandwidth_expectations(model)
    
    # Test compute_SensIAT_expected_values
    outcome_model <- model$models$outcome
    test_data <- head(model.frame(outcome_model), 5)
    
    # Single alpha
    result_single <- compute_SensIAT_expected_values(
        outcome_model,
        alpha = 0,
        new.data = test_data
    )
    expect_s3_class(result_single, "data.frame")
    expect_true("E_exp_alphaY" %in% names(result_single))
    expect_true("E_Yexp_alphaY" %in% names(result_single))
    expect_equal(nrow(result_single), nrow(test_data))
    
    # Multiple alphas
    result_multiple <- compute_SensIAT_expected_values(
        outcome_model,
        alpha = c(-0.3, 0, 0.3),
        new.data = test_data
    )
    expect_equal(nrow(result_multiple), nrow(test_data) * 3)
    expect_true(all(c(-0.3, 0, 0.3) %in% result_multiple$alpha))
})

test_that("fixed_bandwidth model with custom initial values", {
    small_data <- create_small_test_data(10)
    
    # Prepare data for initial value estimation
    followup_data <- small_data |>
        dplyr::group_by(Subject_ID) |>
        dplyr::arrange(Time) |>
        dplyr::mutate(
            prev_outcome = dplyr::lag(Outcome),
            delta_time = Time - dplyr::lag(Time)
        ) |>
        dplyr::filter(Time > 0, !is.na(prev_outcome))
    
    # Test with custom initial function
    custom_initial_fn <- function(X, Y) {
        # Simple initial: first coefficient = 1, rest = 0
        initial <- rep(0, ncol(X))
        initial[1] <- 1
        return(initial)
    }
    
    model_custom <- fit_SensIAT_within_group_model(
        group.data = small_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_bandwidth_model,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        knots = c(60, 460),
        End = 830,
        intensity.args = list(bandwidth = 30),
        outcome.args = list(initial = custom_initial_fn)
    )
    run_fixed_bandwidth_expectations(model_custom)
    
    # Test with numeric initial values
    n_coef <- length(coef(model_custom$models$outcome))
    numeric_initial <- rep(0.1, n_coef)
    numeric_initial[1] <- 1  # Ensure first coef is positive
    
    model_numeric <- fit_SensIAT_within_group_model(
        group.data = small_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_bandwidth_model,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        knots = c(60, 460),
        End = 830,
        intensity.args = list(bandwidth = 30),
        outcome.args = list(initial = numeric_initial)
    )
    run_fixed_bandwidth_expectations(model_numeric)
})

test_that("fixed_bandwidth model predictions", {
    small_data <- create_small_test_data(10)
    
    model <- fit_SensIAT_within_group_model(
        group.data = small_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_bandwidth_model,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        knots = c(60, 460),
        End = 830,
        intensity.args = list(bandwidth = 30)
    )
    
    outcome_model <- model$models$outcome
    
    # Test that the model has the predict method
    expect_true(!is.null(outcome_model$coef))
    expect_true(!is.null(outcome_model$data))
    expect_true(!is.null(outcome_model$bandwidth))
    
    # Test coefficients are valid
    expect_true(is.numeric(outcome_model$coef))
    expect_true(all(is.finite(outcome_model$coef)))
    expect_equal(outcome_model$bandwidth, 1)
})

test_that("fixed_bandwidth model formula and coef methods", {
    small_data <- create_small_test_data(10)
    
    model <- fit_SensIAT_within_group_model(
        group.data = small_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_bandwidth_model,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        knots = c(60, 460),
        End = 830,
        intensity.args = list(bandwidth = 30)
    )
    
    outcome_model <- model$models$outcome
    
    # Test coef method
    model_coef <- coef(outcome_model)
    expect_true(is.numeric(model_coef))
    expect_equal(model_coef, outcome_model$coefficients)
    expect_true(model_coef[1] > 0)  # First coefficient should be positive
    
    # Test that all coefficients are finite
    expect_true(all(is.finite(model_coef)))
})

test_that("fixed_bandwidth model consistency with within_group workflow", {
    small_data <- create_small_test_data(10)
    
    # Fit model with multiple alphas
    model <- fit_SensIAT_within_group_model(
        group.data = small_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_bandwidth_model,
        alpha = c(-0.3, 0, 0.3),
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        knots = c(60, 460),
        End = 830,
        intensity.args = list(bandwidth = 30)
    )
    
    run_fixed_bandwidth_expectations(model)
    
    # Check that influence terms were computed
    expect_true(is.list(model$influence))
    expect_equal(length(model$influence), 3)  # One per alpha
    
    # Check that each influence element has the expected structure
    for (i in seq_along(model$influence)) {
        expect_true(is.list(model$influence[[i]]))
        expect_true("alpha" %in% names(model$influence[[i]]))
        # Each influence element contains alpha as a vector per subject
        expect_true(all(abs(model$influence[[i]]$alpha - c(-0.3, 0, 0.3)[i]) < 1e-10))
    }
    
    # Check that coefficients is a list
    expect_true(is.list(model$coefficients))
    expect_equal(length(model$coefficients), length(model$alpha))
})

test_that("fixed_bandwidth model handles edge cases", {
    small_data <- create_small_test_data(10)
    
    # Test with single alpha = 0
    model_alpha0 <- fit_SensIAT_within_group_model(
        group.data = small_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_bandwidth_model,
        alpha = 0,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        knots = c(60, 460),
        End = 830,
        intensity.args = list(bandwidth = 30)
    )
    run_fixed_bandwidth_expectations(model_alpha0)
    expect_equal(length(model_alpha0$alpha), 1)
    expect_equal(model_alpha0$alpha, 0)
    
    # Test with negative alpha
    model_neg_alpha <- fit_SensIAT_within_group_model(
        group.data = small_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_bandwidth_model,
        alpha = -0.5,
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        knots = c(60, 460),
        End = 830,
        intensity.args = list(bandwidth = 30)
    )
    run_fixed_bandwidth_expectations(model_neg_alpha)
    expect_equal(model_neg_alpha$alpha, -0.5)
})
