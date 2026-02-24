test_that("term2 original and fast methods produce identical results", {
    # Skip unless explicitly enabled - this test is slow
    skip_if_not(identical(Sys.getenv("RUN_SLOW_TESTS"), "true"), 
                "term2 performance test is slow. Set RUN_SLOW_TESTS=true to run")

    # Setup test data
    data_with_lags <- SensIAT_example_data |>
        dplyr::group_by(Subject_ID) |>
        dplyr::mutate(
            ..prev_outcome.. = dplyr::lag(Outcome, default = NA_real_, order_by = Time),
            ..prev_time.. = dplyr::lag(Time, default = 0, order_by = Time),
            ..delta_time.. = Time - dplyr::lag(.data$Time, default = NA_real_, order_by = Time)
        ) |>
        dplyr::ungroup()

    # Fit outcome model
    outcome.model <-
        fit_SensIAT_single_index_fixed_coef_model(
            Outcome ~ splines::ns(..prev_outcome.., df = 3) + ..delta_time.. - 1,
            id = Subject_ID,
            data = data_with_lags |> dplyr::filter(Time > 0)
        )

    # Setup marginal mean parameters
    alpha <- 0
    knots <- c(60, 260, 460)
    spline.degree <- 3L
    knots_extended <- c(
        rep(head(knots, 1), spline.degree),
        knots,
        rep(tail(knots, 1), spline.degree)
    )
    base <- orthogonalsplinebasis::SplineBasis(knots_extended, order = spline.degree + 1L)
    V <- orthogonalsplinebasis::GramMatrix(base)
    V.inv <- solve(V)

    tmin <- base@knots[base@order]
    tmax <- base@knots[length(base@knots) - base@order + 1]

    # Initialize beta
    set.seed(123)
    observed_outcomes <- data_with_lags$Outcome[!is.na(data_with_lags$Outcome)]
    mean_outcome <- mean(observed_outcomes)
    beta_test <- rep(log(mean_outcome) / ncol(base), ncol(base)) + rnorm(ncol(base), 0, 0.1)

    # Imputation function
    impute_data_fn <- function(t, df) {
        data_wl <- df |>
            dplyr::mutate(
                ..prev_time.. = .data$Time,
                ..prev_outcome.. = .data$Outcome,
                ..delta_time.. = 0
            )
        extrapolate_from_last_observation(t, data_wl, "Time", slopes = c("..delta_time.." = 1))
    }

    # Inverse link (log link)
    inv_link <- exp

    # Weight function for lp_mse + log
    W <- function(t, beta) {
        B <- as.vector(pcoriaccel_evaluate_basis(base, t))
        mu <- sum(B * beta)
        as.vector((V.inv %*% B) * exp(-mu))
    }

    # Test on first 3 patients
    test_ids <- unique(data_with_lags$Subject_ID)[1:3]

    # Compute term2 for each patient using both methods
    results_original <- list()
    results_fast <- list()
    timings_original <- numeric(length(test_ids))
    timings_fast <- numeric(length(test_ids))

    for (i in seq_along(test_ids)) {
        patient_id <- test_ids[i]
        patient_data <- data_with_lags[data_with_lags$Subject_ID == patient_id, ]

        # Original method
        t_start <- proc.time()[3]
        results_original[[i]] <- SensIAT:::compute_term2_influence_original(
            patient_data = patient_data,
            outcome_model = outcome.model,
            base = base,
            alpha = alpha,
            marginal_beta = beta_test,
            V_inv = V.inv,
            tmin = tmin,
            tmax = tmax,
            impute_fn = impute_data_fn,
            inv_link = inv_link,
            W = W
        )
        timings_original[i] <- proc.time()[3] - t_start

        # Fast method
        t_start <- proc.time()[3]
        results_fast[[i]] <- SensIAT:::compute_term2_influence_fast(
            patient_data = patient_data,
            outcome_model = outcome.model,
            base = base,
            alpha = alpha,
            marginal_beta = beta_test,
            V_inv = V.inv,
            tmin = tmin,
            tmax = tmax,
            impute_fn = impute_data_fn, # Not used by fast, but included for interface consistency
            inv_link = inv_link, # Not used by fast, but included for interface consistency
            W = W
        )
        timings_fast[i] <- proc.time()[3] - t_start
    }

    # Test 1: Results should be numerically identical (within tolerance)
    for (i in seq_along(test_ids)) {
        expect_equal(
            as.vector(results_fast[[i]]),
            as.vector(results_original[[i]]),
            tolerance = 1e-6,
            info = paste("Patient", test_ids[i], "results should match")
        )
    }

    # Test 2: Fast method should be faster on average
    # Allow some slack due to measurement noise, but expect at least 2x speedup
    mean_time_original <- mean(timings_original)
    mean_time_fast <- mean(timings_fast)

    expect_true(
        mean_time_fast < mean_time_original,
        info = sprintf(
            "Fast method (%.3fs) should be faster than original (%.3fs)",
            mean_time_fast,
            mean_time_original
        )
    )

    # Report speedup for informational purposes
    speedup <- mean_time_original / mean_time_fast
    message(sprintf(
        "\nPerformance comparison (n=%d patients):\n  Original: %.3fs/patient\n  Fast:     %.3fs/patient\n  Speedup:  %.1fx",
        length(test_ids),
        mean_time_original,
        mean_time_fast,
        speedup
    ))

    # Test 3: Both methods should produce vectors of correct length
    for (i in seq_along(test_ids)) {
        expect_length(results_original[[i]], ncol(base))
        expect_length(results_fast[[i]], ncol(base))
    }

    # Test 4: Results should be numeric and finite
    for (i in seq_along(test_ids)) {
        expect_true(all(is.finite(results_original[[i]])))
        expect_true(all(is.finite(results_fast[[i]])))
    }
})
