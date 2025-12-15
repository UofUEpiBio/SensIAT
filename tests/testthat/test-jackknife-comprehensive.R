# Comprehensive Jackknife Tests
# These tests use the full dataset and are more computationally expensive
# Run these manually or in CI when doing major jackknife changes
#
# To run these tests:
# testthat::test_file("tests/testthat/test-jackknife-comprehensive.R")

test_that("jackknife comprehensive test with full dataset", {
    skip("Manual test - run explicitly when needed")
    skip_on_cran()

    # This is the original expensive test
    model <- fit_SensIAT_within_group_model(
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

    rlang::check_installed("future")
    rlang::check_installed("furrr")

    future::plan(future::sequential)
    jk1 <- jackknife(model, time = c(90, 180, 270, 360, 450))

    future::plan(future::multisession, workers = 2)
    jk2 <- jackknife(model, time = c(90, 180, 270, 360, 450))

    expect_identical(jk1, jk2)
})

test_that("jackknife medium scale test", {
    skip("Run when changes affect jackknife core logic")
    skip_on_cran()

    # Medium-scale test with 50 subjects
    medium_subjects <- head(unique(SensIAT_example_data$Subject_ID), 50)
    medium_data <- SensIAT_example_data |>
        dplyr::filter(Subject_ID %in% medium_subjects)

    model <- fit_SensIAT_within_group_model(
        group.data = medium_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        alpha = c(-0.3, 0, 0.3),
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        End = 830,
        knots = c(60, 460),
        outcome.args = list(
            model = ~ ns(..prev_outcome.., knots = c(9 / 6)) + scale(..delta_time..) - 1
        )
    )

    jk_result <- jackknife(model, time = c(180, 360))

    # Comprehensive validation
    expect_s3_class(jk_result, "SensIAT_withingroup_jackknife_results")
    expect_equal(length(unique(jk_result$alpha)), 3)
    expect_equal(length(unique(jk_result$time)), 2)
    expect_true(all(!is.na(jk_result$jackknife_var)))
    expect_true(all(jk_result$jackknife_var > 0))

    # Check that variance estimates seem reasonable
    expect_true(all(jk_result$jackknife_var < 100)) # Not unreasonably large
    expect_true(all(jk_result$jackknife_var > 1e-6)) # Not unreasonably small
})

test_that("jackknife performance benchmark", {
    skip("Performance benchmark - run manually")

    # Benchmark the speed improvement
    small_subjects <- head(unique(SensIAT_example_data$Subject_ID), 10)
    small_data <- SensIAT_example_data |>
        dplyr::filter(Subject_ID %in% small_subjects)

    model <- fit_SensIAT_within_group_model(
        group.data = small_data,
        outcome_modeler = fit_SensIAT_single_index_fixed_coef_model,
        alpha = c(0),
        id = Subject_ID,
        outcome = Outcome,
        time = Time,
        End = 830,
        knots = c(260),
        outcome.args = list(model = ~ ..prev_outcome.. - 1)
    )

    # Time the jackknife operation
    timing <- system.time({
        jk_result <- jackknife(model, time = c(180))
    })

    cat("Jackknife with 10 subjects took:", timing[["elapsed"]], "seconds\n")

    # Should complete in reasonable time (under 2 minutes for 10 subjects)
    expect_true(timing[["elapsed"]] < 120)
})
