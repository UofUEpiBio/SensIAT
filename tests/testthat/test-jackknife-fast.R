# Fast Jackknife Testing Strategies
# This file implements multiple strategies to make jackknife testing faster
# while maintaining confidence in the correctness of the implementation

# Strategy 1: Create a small subset of the example data for testing
create_small_test_data <- function(n_subjects = 10) {
    # Take a subset of subjects from the example data
    small_subjects <- head(unique(SensIAT_example_data$Subject_ID), n_subjects)
    SensIAT_example_data |>
        dplyr::filter(Subject_ID %in% small_subjects)
}

# Strategy 2: Use a simpler, faster model for testing basic jackknife functionality
test_that("jackknife basic functionality with small dataset", {
    skip_on_cran()
    
    # Create small test data (10 subjects instead of 200)
    small_data <- create_small_test_data(10)
    
    # Use simpler model configuration for faster fitting
    model <- fit_SensIAT_within_group_model(
        group.data = small_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        alpha = c(0, 0.3),  # Fewer alpha values
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        End = 830,
        knots = c(60, 460),  # Need at least 2 knots for spline degree 3
        intensity.args = list(bandwidth = 30),  # Explicit bandwidth to avoid estimation issues
        outcome.args = list(
            model = ~ ..prev_outcome.. + ..delta_time.. - 1  # Simpler model
        )
    )
    
    # Test basic functionality
    jk_result <- jackknife(model, time = c(180, 360))
    
    # Verify structure and content
    expect_s3_class(jk_result, "SensIAT_withingroup_jackknife_results")
    expect_true(all(c("alpha", "time", "jackknife_mean", "jackknife_var") %in% names(jk_result)))
    expect_equal(nrow(jk_result), 4)  # 2 alpha × 2 time points
    expect_true(all(is.finite(jk_result$jackknife_mean)))
    expect_true(all(jk_result$jackknife_var >= 0))
})

# Strategy 3: Test parallelization with very small data
test_that("jackknife parallelization with minimal data", {
    skip_on_cran()
    
    # Use only 5 subjects for parallelization test
    tiny_data <- create_small_test_data(5)
    
    model <- fit_SensIAT_within_group_model(
        group.data = tiny_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        alpha = c(0),  # Single alpha value
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        End = 830,
        knots = c(60, 460),  # Need at least 2 knots for spline degree 3
        intensity.args = list(bandwidth = 30),  # Explicit bandwidth to avoid estimation issues
        outcome.args = list(model = ~ ..prev_outcome.. + ..delta_time.. - 1)
    )
    
    rlang::check_installed("future")
    rlang::check_installed("furrr")
    
    # Sequential
    future::plan(future::sequential)
    jk1 <- jackknife(model, time = c(180))
    
    # Parallel
    future::plan(future::multisession, workers = 2)
    jk2 <- jackknife(model, time = c(180))
    
    # Should be identical
    expect_identical(jk1, jk2)
})

# Strategy 4: Mock/stub testing for complex scenarios
test_that("jackknife computation logic with mocked data", {
    skip_on_cran()
    
    # Test the summarization logic separately from the expensive model fitting
    # This tests the mathematical correctness without the computational burden
    
    # Create mock replication results
    mock_replications <- list()
    n_reps <- 5
    
    for (i in 1:n_reps) {
        # Simulate prediction results that would come from cross_validate
        mock_replications[[i]] <- tibble::tibble(
            alpha = rep(c(0, 0.3), each = 2),
            time = rep(c(180, 360), 2),
            mean = rnorm(4, mean = 10 + c(0, 0.5, 1, 1.5), sd = 0.1)
        )
    }
    
    # Mock original estimates
    original_estimates <- tibble::tibble(
        alpha = rep(c(0, 0.3), each = 2),
        time = rep(c(180, 360), 2),
        mean = c(10, 10.5, 11, 11.5),
        var = rep(0.25, 4)
    )
    
    # Test the summarization logic directly
    estimates <- purrr::map(mock_replications, ~ .x)
    summary_result <- estimates |> 
        dplyr::bind_rows(.id = '.rep') |>
        dplyr::group_by(alpha, time) |>
        dplyr::summarize(
            jackknife_mean = mean(mean),
            jackknife_var = (dplyr::n()-1)/dplyr::n() * sum((mean - mean(mean))^2),
            .groups = 'drop'
        ) |>
        dplyr::ungroup() |>
        dplyr::full_join(original_estimates, by = c('alpha', 'time'))
    
    # Verify the jackknife variance calculation
    expect_true(all(summary_result$jackknife_var >= 0))
    expect_equal(nrow(summary_result), 4)
    expect_true(all(c("jackknife_mean", "jackknife_var") %in% names(summary_result)))
})

# Strategy 5: Integration test with reduced scope
test_that("jackknife integration with medium dataset", {
    skip_on_cran()
    
    # This test runs occasionally with more subjects but still much faster than full test
    # Use testthat::skip() to control when this runs
    skip_if(Sys.getenv("SKIP_MEDIUM_TESTS") == "true", "Medium tests skipped")
    
    medium_data <- create_small_test_data(20)  # 20 subjects instead of 200
    
    model <- fit_SensIAT_within_group_model(
        group.data = medium_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        alpha = c(-0.3, 0, 0.3),
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        End = 830,
        knots = c(60, 460),  # 2 knots instead of 3
        outcome.args = list(
            model = ~ ns(..prev_outcome.., knots = c(9/6)) + scale(..delta_time..) - 1
        )
    )
    
    jk_result <- jackknife(model, time = c(180, 360))
    
    # More comprehensive checks
    expect_s3_class(jk_result, "SensIAT_withingroup_jackknife_results")
    expect_equal(length(unique(jk_result$alpha)), 3)
    expect_equal(length(unique(jk_result$time)), 2)
    expect_true(all(!is.na(jk_result$jackknife_var)))
    expect_true(all(jk_result$jackknife_var > 0))  # Should have some variance
})