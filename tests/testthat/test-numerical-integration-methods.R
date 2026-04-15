# Tests for numerical and piecewise integration methods
# These are alternative integration methods used when integration.method = "numerical" or "piecewise"

test_that("numerical integration method executes without error", {
    data("SensIAT_example_data", package = "SensIAT", envir = environment())
    
    # Prepare data (same pattern as other integration tests)
    model.data <- SensIAT_example_data |>
        dplyr::group_by(Subject_ID) |>
        dplyr::arrange(Time) |>
        dplyr::mutate(
            prev_time = dplyr::lag(Time),
            prev_outcome = dplyr::lag(Outcome),
            delta_time = Time - dplyr::lag(Time),
            visit.number = seq_along(Time)
        ) |>
        dplyr::filter(!is.na(Outcome))

    followup.data <- model.data |>
        dplyr::filter(Time > 0)

    # Fit models
    intensity.model <-
        rlang::inject(coxph(
            Surv(prev_time, Time, !is.na(Outcome)) ~
                prev_outcome + strata(visit.number),
            id = Subject_ID,
            data = followup.data
        ))

    outcome.model <- fit_SensIAT_single_index_fixed_coef_model(
        Outcome ~
            splines::ns(prev_outcome, df = 3) +
            prev_outcome +
            delta_time - 1,
        id = Subject_ID,
        data = followup.data
    )

    base <- orthogonalsplinebasis::SplineBasis(c(60, 60, 60, 60, 260, 460, 460, 460, 460))

    centering.statistics <-
        dplyr::summarize(
            dplyr::ungroup(dplyr::filter(model.data, Time > 0, !is.na(Outcome))),
            dplyr::across(
                c(Time, delta_time),
                list(
                    mean = ~ mean(.x, na.rm = TRUE),
                    sd = ~ sd(.x, na.rm = TRUE)
                )
            )
        ) |>
        as.numeric() |>
        matrix(ncol = 2, byrow = TRUE) |>
        `dimnames<-`(list(c("time", "delta_time"), c("mean", "sd")))

    df_i <- model.data |>
        dplyr::filter(Subject_ID == 1)

    variables <- list(
        time = rlang::sym("Time"),
        id = rlang::sym("Subject_ID"),
        outcome = rlang::sym("Outcome"),
        prev_time = rlang::sym("prev_time"),
        prev_outcome = rlang::sym("prev_outcome"),
        delta_time = rlang::sym("delta_time")
    )

    # Test numerical integration method with alpha = 0
    result_numerical <- compute_influence_for_one_alpha_and_one_patient(
        df_i,
        alpha = 0,
        variables = variables,
        intensity.model = intensity.model,
        outcome.model = outcome.model,
        base = base,
        control = pcori_control(
            integration.method = "numerical",
            resolution = 100  # Lower resolution for speed
        ),
        centering = centering.statistics
    )

    # Verify structure
    expect_s3_class(result_numerical, "data.frame")
    expect_true("term1" %in% names(result_numerical))
    expect_true("term2" %in% names(result_numerical))
    expect_equal(nrow(result_numerical), 1)
    
    # Verify term2 is a numeric vector with correct length
    expect_true(is.numeric(result_numerical$term2[[1]]))
    expect_equal(length(result_numerical$term2[[1]]), ncol(base))
    
    # Verify values are finite
    expect_true(all(is.finite(result_numerical$term1[[1]])))
    expect_true(all(is.finite(result_numerical$term2[[1]])))
})

test_that("piecewise integration method executes without error", {
    data("SensIAT_example_data", package = "SensIAT", envir = environment())
    
    # Use same setup as numerical test
    model.data <- SensIAT_example_data |>
        dplyr::group_by(Subject_ID) |>
        dplyr::arrange(Time) |>
        dplyr::mutate(
            prev_time = dplyr::lag(Time),
            prev_outcome = dplyr::lag(Outcome),
            delta_time = Time - dplyr::lag(Time),
            visit.number = seq_along(Time)
        ) |>
        dplyr::filter(!is.na(Outcome))

    followup.data <- model.data |>
        dplyr::filter(Time > 0)

    intensity.model <-
        rlang::inject(coxph(
            Surv(prev_time, Time, !is.na(Outcome)) ~
                prev_outcome + strata(visit.number),
            id = Subject_ID,
            data = followup.data
        ))

    outcome.model <- fit_SensIAT_single_index_fixed_coef_model(
        Outcome ~
            splines::ns(prev_outcome, df = 3) +
            prev_outcome +
            delta_time - 1,
        id = Subject_ID,
        data = followup.data
    )

    base <- orthogonalsplinebasis::SplineBasis(c(60, 60, 60, 60, 260, 460, 460, 460, 460))

    centering.statistics <-
        dplyr::summarize(
            dplyr::ungroup(dplyr::filter(model.data, Time > 0, !is.na(Outcome))),
            dplyr::across(
                c(Time, delta_time),
                list(
                    mean = ~ mean(.x, na.rm = TRUE),
                    sd = ~ sd(.x, na.rm = TRUE)
                )
            )
        ) |>
        as.numeric() |>
        matrix(ncol = 2, byrow = TRUE) |>
        `dimnames<-`(list(c("time", "delta_time"), c("mean", "sd")))

    df_i <- model.data |>
        dplyr::filter(Subject_ID == 1)

    variables <- list(
        time = rlang::sym("Time"),
        id = rlang::sym("Subject_ID"),
        outcome = rlang::sym("Outcome"),
        prev_time = rlang::sym("prev_time"),
        prev_outcome = rlang::sym("prev_outcome"),
        delta_time = rlang::sym("delta_time")
    )

    # Test piecewise integration method
    result_piecewise <- compute_influence_for_one_alpha_and_one_patient(
        df_i,
        alpha = 0,
        variables = variables,
        intensity.model = intensity.model,
        outcome.model = outcome.model,
        base = base,
        control = pcori_control(
            integration.method = "piecewise",
            resolution.within.period = 20  # Lower resolution for speed
        ),
        centering = centering.statistics
    )

    # Verify structure
    expect_s3_class(result_piecewise, "data.frame")
    expect_true("term1" %in% names(result_piecewise))
    expect_true("term2" %in% names(result_piecewise))
    expect_equal(nrow(result_piecewise), 1)
    
    # Verify term2 is a numeric vector with correct length
    expect_true(is.numeric(result_piecewise$term2[[1]]))
    expect_equal(length(result_piecewise$term2[[1]]), ncol(base))
    
    # Verify values are finite
    expect_true(all(is.finite(result_piecewise$term1[[1]])))
    expect_true(all(is.finite(result_piecewise$term2[[1]])))
})

test_that("numerical and piecewise methods produce reasonable results compared to quadv", {
    data("SensIAT_example_data", package = "SensIAT", envir = environment())
    
    # Simplified setup for comparison test
    model.data <- SensIAT_example_data |>
        dplyr::group_by(Subject_ID) |>
        dplyr::arrange(Time) |>
        dplyr::mutate(
            prev_time = dplyr::lag(Time),
            prev_outcome = dplyr::lag(Outcome),
            delta_time = Time - dplyr::lag(Time),
            visit.number = seq_along(Time)
        ) |>
        dplyr::filter(!is.na(Outcome))

    followup.data <- model.data |>
        dplyr::filter(Time > 0)

    intensity.model <-
        rlang::inject(coxph(
            Surv(prev_time, Time, !is.na(Outcome)) ~
                prev_outcome + strata(visit.number),
            id = Subject_ID,
            data = followup.data
        ))

    outcome.model <- fit_SensIAT_single_index_fixed_coef_model(
        Outcome ~
            splines::ns(prev_outcome, df = 3) +
            prev_outcome +
            delta_time - 1,
        id = Subject_ID,
        data = followup.data
    )

    base <- orthogonalsplinebasis::SplineBasis(c(60, 60, 60, 60, 260, 460, 460, 460, 460))

    centering.statistics <-
        dplyr::summarize(
            dplyr::ungroup(dplyr::filter(model.data, Time > 0, !is.na(Outcome))),
            dplyr::across(
                c(Time, delta_time),
                list(
                    mean = ~ mean(.x, na.rm = TRUE),
                    sd = ~ sd(.x, na.rm = TRUE)
                )
            )
        ) |>
        as.numeric() |>
        matrix(ncol = 2, byrow = TRUE) |>
        `dimnames<-`(list(c("time", "delta_time"), c("mean", "sd")))

    df_i <- model.data |>
        dplyr::filter(Subject_ID == 1)

    variables <- list(
        time = rlang::sym("Time"),
        id = rlang::sym("Subject_ID"),
        outcome = rlang::sym("Outcome"),
        prev_time = rlang::sym("prev_time"),
        prev_outcome = rlang::sym("prev_outcome"),
        delta_time = rlang::sym("delta_time")
    )

    # Compute with all three methods
    result_quadv <- compute_influence_for_one_alpha_and_one_patient(
        df_i, alpha = 0,
        variables = variables,
        intensity.model = intensity.model,
        outcome.model = outcome.model,
        base = base,
        control = pcori_control(integration.method = "quadv"),
        centering = centering.statistics
    )

    result_numerical <- compute_influence_for_one_alpha_and_one_patient(
        df_i, alpha = 0,
        variables = variables,
        intensity.model = intensity.model,
        outcome.model = outcome.model,
        base = base,
        control = pcori_control(integration.method = "numerical", resolution = 500),
        centering = centering.statistics
    )

    result_piecewise <- compute_influence_for_one_alpha_and_one_patient(
        df_i, alpha = 0,
        variables = variables,
        intensity.model = intensity.model,
        outcome.model = outcome.model,
        base = base,
        control = pcori_control(integration.method = "piecewise", resolution.within.period = 50),
        centering = centering.statistics
    )

    # All methods should produce similar term1 (this doesn't use integration)
    expect_equal(result_quadv$term1[[1]], result_numerical$term1[[1]], tolerance = 1e-6)
    expect_equal(result_quadv$term1[[1]], result_piecewise$term1[[1]], tolerance = 1e-6)

    # Term2 should be reasonably close (numerical methods less accurate but should be in ballpark)
    # Use relative tolerance since absolute values can vary
    term2_quadv <- result_quadv$term2[[1]]
    term2_numerical <- result_numerical$term2[[1]]
    term2_piecewise <- result_piecewise$term2[[1]]

    # Check that numerical methods are within 10% of quadv (relaxed tolerance for numerical integration)
    relative_diff_numerical <- abs(term2_numerical - term2_quadv) / (abs(term2_quadv) + 1e-10)
    relative_diff_piecewise <- abs(term2_piecewise - term2_quadv) / (abs(term2_quadv) + 1e-10)

    # Most elements should be reasonably close (relaxed threshold for piecewise method which can be less accurate)
    expect_true(mean(relative_diff_numerical < 0.1) > 0.8, 
                info = "Numerical method should match quadv for most elements")
    expect_true(mean(relative_diff_piecewise < 0.2) > 0.6,
                info = "Piecewise method should match quadv reasonably (within 20% for 60% of elements)")
})

test_that("numerical integration methods work with non-zero alpha", {
    data("SensIAT_example_data", package = "SensIAT", envir = environment())
    
    # Minimal setup
    model.data <- SensIAT_example_data |>
        dplyr::group_by(Subject_ID) |>
        dplyr::arrange(Time) |>
        dplyr::mutate(
            prev_time = dplyr::lag(Time),
            prev_outcome = dplyr::lag(Outcome),
            delta_time = Time - dplyr::lag(Time),
            visit.number = seq_along(Time)
        ) |>
        dplyr::filter(!is.na(Outcome))

    followup.data <- model.data |>
        dplyr::filter(Time > 0)

    intensity.model <-
        rlang::inject(coxph(
            Surv(prev_time, Time, !is.na(Outcome)) ~
                prev_outcome + strata(visit.number),
            id = Subject_ID,
            data = followup.data
        ))

    outcome.model <- fit_SensIAT_single_index_fixed_coef_model(
        Outcome ~
            splines::ns(prev_outcome, df = 3) +
            prev_outcome +
            delta_time - 1,
        id = Subject_ID,
        data = followup.data
    )

    base <- orthogonalsplinebasis::SplineBasis(c(60, 60, 60, 60, 260, 460, 460, 460, 460))

    centering.statistics <-
        dplyr::summarize(
            dplyr::ungroup(dplyr::filter(model.data, Time > 0, !is.na(Outcome))),
            dplyr::across(
                c(Time, delta_time),
                list(
                    mean = ~ mean(.x, na.rm = TRUE),
                    sd = ~ sd(.x, na.rm = TRUE)
                )
            )
        ) |>
        as.numeric() |>
        matrix(ncol = 2, byrow = TRUE) |>
        `dimnames<-`(list(c("time", "delta_time"), c("mean", "sd")))

    df_i <- model.data |>
        dplyr::filter(Subject_ID == 1)

    variables <- list(
        time = rlang::sym("Time"),
        id = rlang::sym("Subject_ID"),
        outcome = rlang::sym("Outcome"),
        prev_time = rlang::sym("prev_time"),
        prev_outcome = rlang::sym("prev_outcome"),
        delta_time = rlang::sym("delta_time")
    )

    # Test with alpha = 0.3
    result_numerical_pos <- compute_influence_for_one_alpha_and_one_patient(
        df_i, alpha = 0.3,
        variables = variables,
        intensity.model = intensity.model,
        outcome.model = outcome.model,
        base = base,
        control = pcori_control(integration.method = "numerical", resolution = 100),
        centering = centering.statistics
    )

    result_piecewise_pos <- compute_influence_for_one_alpha_and_one_patient(
        df_i, alpha = 0.3,
        variables = variables,
        intensity.model = intensity.model,
        outcome.model = outcome.model,
        base = base,
        control = pcori_control(integration.method = "piecewise", resolution.within.period = 20),
        centering = centering.statistics
    )

    # Basic structure checks
    expect_true(all(is.finite(result_numerical_pos$term1[[1]])))
    expect_true(all(is.finite(result_numerical_pos$term2[[1]])))
    expect_true(all(is.finite(result_piecewise_pos$term1[[1]])))
    expect_true(all(is.finite(result_piecewise_pos$term2[[1]])))
    
    # Test with negative alpha
    result_numerical_neg <- compute_influence_for_one_alpha_and_one_patient(
        df_i, alpha = -0.3,
        variables = variables,
        intensity.model = intensity.model,
        outcome.model = outcome.model,
        base = base,
        control = pcori_control(integration.method = "numerical", resolution = 100),
        centering = centering.statistics
    )

    expect_true(all(is.finite(result_numerical_neg$term1[[1]])))
    expect_true(all(is.finite(result_numerical_neg$term2[[1]])))
})
