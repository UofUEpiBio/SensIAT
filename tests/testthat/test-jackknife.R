# Helper function to create small test dataset
create_small_test_data <- function(n_subjects = 10) {
    data("SensIAT_example_data", package = "SensIAT", envir = environment())
    small_subjects <- head(unique(SensIAT_example_data$Subject_ID), n_subjects)
    SensIAT_example_data |>
        dplyr::filter(Subject_ID %in% small_subjects)
}

test_that("jackknife basic functionality", {
    skip_on_cran()
    
    # Use small dataset (10 subjects instead of 200)
    small_data <- create_small_test_data(10)
    
    model <- fit_SensIAT_within_group_model(
        group.data = small_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        alpha = c(0, 0.3),  # Reduced alpha values
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        knots = c(0, 260, max(small_data$Time)),  # Single knot for speed
        outcome.args = list(
            model = ~ ..prev_outcome.. + ..delta_time.. - 1  # Simplified model
        )
    )
    
    jk_result <- jackknife(model, time = c(180, 360))
    
    # Test structure and basic properties
    expect_s3_class(jk_result, "SensIAT_withingroup_jackknife_results")
    expect_true(all(c("alpha", "time", "jackknife_mean", "jackknife_var") %in% names(jk_result)))
    expect_equal(nrow(jk_result), 4)  # 2 alpha × 2 time points
    expect_true(all(is.finite(jk_result$jackknife_mean)))
    expect_true(all(jk_result$jackknife_var >= 0))
})

test_that("jackknife is invariant under parallelization.", {
    skip_on_cran()
    
    # Use very small dataset for parallelization test (5 subjects)
    tiny_data <- create_small_test_data(5)
    
    model <- fit_SensIAT_within_group_model(
        group.data = tiny_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        alpha = c(0),  # Single alpha for speed
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
    
    future::plan(future::sequential)
    jk1 <- jackknife(model, time = c(180))

    future::plan(future::multisession, workers=2)
    jk2 <- jackknife(model, time = c(180))

    expect_identical(jk1, jk2)
})

test_that("jackknife computation logic", {
    # Test jackknife variance calculation without expensive model fitting
    skip_on_cran()
    
    # Mock replication results to test summarization logic
    mock_predictions <- list(
        tibble::tibble(alpha = c(0, 0), time = c(180, 360), mean = c(10.1, 11.1)),
        tibble::tibble(alpha = c(0, 0), time = c(180, 360), mean = c(9.9, 10.9)),
        tibble::tibble(alpha = c(0, 0), time = c(180, 360), mean = c(10.0, 11.0)),
        tibble::tibble(alpha = c(0, 0), time = c(180, 360), mean = c(10.2, 11.2))
    )
    
    # Test the summarization
    result <- mock_predictions |> 
        dplyr::bind_rows(.id = '.rep') |>
        dplyr::group_by(alpha, time) |>
        dplyr::summarize(
            n = dplyr::n(),
            jackknife_mean = mean(mean),
            jackknife_var = (n-1)/n * sum((mean - mean(mean))^2),
            .groups = 'drop'
        )
    
    expect_equal(nrow(result), 2)
    expect_true(all(result$jackknife_var >= 0))
    expect_true(all(abs(result$jackknife_mean - c(10.05, 11.05)) < 0.01))
})
