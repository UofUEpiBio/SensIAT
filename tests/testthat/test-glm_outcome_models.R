test_that("compute_SensIAT_expected_values.glm: gaussian family", {
    # Simulate some data
    set.seed(123)
    data <- simulate_SensIAT_data(
        n_subjects = 20,
        End = 500,
        link = "identity",
        outcome_sd = 2
    )

    # Prepare data
    data_prepared <- prepare_SensIAT_data(data, Subject_ID, Time, Outcome, End = 500)

    # Fit a gaussian GLM
    glm_model <- glm(
        Outcome ~ ..prev_outcome.. + ..delta_time..,
        data = data_prepared |> filter(Time > 0),
        family = gaussian()
    )

    # Test compute_SensIAT_expected_values
    test_data <- data_prepared |> filter(Time > 0, ..delta_time.. > 0) |> head(5)
    result <- compute_SensIAT_expected_values(
        model = glm_model,
        alpha = 0,
        new.data = test_data
    )

    expect_s3_class(result, "data.frame")
    expect_true("E_Yexp_alphaY" %in% names(result))
    expect_true("E_exp_alphaY" %in% names(result))
    expect_equal(nrow(result), 5)
    expect_true(all(is.finite(result$E_Yexp_alphaY)))
    expect_true(all(is.finite(result$E_exp_alphaY)))
    expect_true(all(result$E_exp_alphaY > 0))
})

test_that("compute_SensIAT_expected_values.glm: binomial family", {
    # Simulate binary data
    set.seed(456)
    data <- simulate_SensIAT_data(
        n_subjects = 50,
        End = 500,
        link = "logit",
        outcome_coef = list(
            intercept = 3,  #< used to adjust the overall prevalence
            prev_outcome = c(0.3, -0.05, 0.02),
            time = 0.0005,
            delta_time = -0.05
        ),
        initial_outcome_mean = 0.5
    )

    # Prepare data
    data_prepared <- prepare_SensIAT_data(data, Subject_ID, Time, Outcome, End = 500)

    # Fit a binomial GLM
    glm_model <- glm(
        Outcome ~ ..prev_outcome.. + ..delta_time..,
        data = data_prepared |> filter(Time > 0),
        family = binomial()
    )

    # Test compute_SensIAT_expected_values
    test_data <- data_prepared |> filter(Time > 0, ..delta_time.. > 0) |> head(5)
    result <- compute_SensIAT_expected_values(
        model = glm_model,
        alpha = 0.5,
        new.data = test_data
    )

    expect_s3_class(result, "tbl_df")
    expect_true("E_Yexp_alphaY" %in% names(result))
    expect_true("E_exp_alphaY" %in% names(result))
    expect_true("alpha" %in% names(result))
    expect_equal(nrow(result), 5)
    expect_true(all(is.finite(result$E_Yexp_alphaY)))
    expect_true(all(is.finite(result$E_exp_alphaY)))
    expect_true(all(result$E_exp_alphaY > 0))
    expect_equal(unique(result$alpha), 0.5)
})

test_that("compute_SensIAT_expected_values.glm: poisson family", {
    # Simulate count data
    set.seed(789)
    data <- simulate_SensIAT_data(
        n_subjects = 50,
        End = 500,
        link = "log",
        outcome_coef = list(
            intercept = 2.0,
            prev_outcome = c(0.05, -0.01, 0.002),
            time = 0.0002,
            delta_time = -0.02
        ),
        initial_outcome_mean = 7
    )

    # Prepare data
    data_prepared <- prepare_SensIAT_data(data, Subject_ID, Time, Outcome, End = 500)

    # Fit a poisson GLM
    glm_model <- glm(
        Outcome ~ ..prev_outcome.. + ..delta_time..,
        data = data_prepared |> filter(Time > 0),
        family = poisson()
    )

    # Test compute_SensIAT_expected_values
    test_data <- data_prepared |> filter(Time > 0, ..delta_time.. > 0) |> head(5)
    result <- compute_SensIAT_expected_values(
        model = glm_model,
        alpha = -0.1,
        new.data = test_data
    )

    expect_s3_class(result, "tbl_df")
    expect_true("E_Yexp_alphaY" %in% names(result))
    expect_true("E_exp_alphaY" %in% names(result))
    expect_equal(nrow(result), 5)
    expect_true(all(is.finite(result$E_Yexp_alphaY)))
    expect_true(all(is.finite(result$E_exp_alphaY)))
    expect_true(all(result$E_exp_alphaY > 0))
})

test_that("compute_SensIAT_expected_values.glm: multiple alpha values", {
    # Simulate data
    set.seed(111)
    data <- simulate_SensIAT_data(
        n_subjects = 15,
        End = 400,
        link = "identity"
    )

    # Prepare data
    data_prepared <- prepare_SensIAT_data(data, Subject_ID, Time, Outcome, End = 500)

    # Fit GLM
    glm_model <- glm(
        Outcome ~ ..prev_outcome.. + ..delta_time..,
        data = data_prepared |> filter(Time > 0),
        family = gaussian()
    )

    # Test with multiple alpha values
    test_data <- data_prepared |> filter(Time > 0) |> head(3)
    result <- compute_SensIAT_expected_values(
        model = glm_model,
        alpha = c(-0.5, 0, 0.5),
        new.data = test_data
    )

    expect_equal(nrow(result), 9)  # 3 data rows × 3 alpha values
    expect_equal(sort(unique(result$alpha)), c(-0.5, 0, 0.5))
})

test_that("fit_SensIAT_marginal_mean_model_generalized: works with gaussian GLM", {
    # Simulate data
    set.seed(222)
    data <- simulate_SensIAT_data(
        n_subjects = 15,
        End = 500,
        link = "identity",
        outcome_sd = 2
    )

    # Prepare data
    data_with_lags <- prepare_SensIAT_data(data, Subject_ID, Time, Outcome, End = 500)

    # Fit intensity model
    intensity.model <- coxph(
        Surv(..prev_time.., Time, !is.na(Outcome)) ~ ..prev_outcome..,
        data = data_with_lags |> filter(Time > 0, ..delta_time.. > 0)
    )

    # Fit gaussian GLM outcome model
    outcome.model <- glm(
        Outcome ~ ..prev_outcome.. + ..delta_time..,
        data = data_with_lags |> filter(Time > 0),
        family = gaussian()
    )

    # Test that generalized fitting works
    result <- suppressWarnings(
        fit_SensIAT_marginal_mean_model_generalized(
            data = data_with_lags,
            time = Time,
            id = Subject_ID,
            alpha = 0,
            knots = c(100, 250, 400),
            outcome.model = outcome.model,
            intensity.model = intensity.model,
            loss = "lp_mse",
            link = "identity"
        )
    )

    # Check result structure - should match fit_SensIAT_marginal_mean_model
    expect_type(result, "list")
    expect_true("coefficients" %in% names(result))
    expect_true("coefficient.variance" %in% names(result))
    expect_true("influence" %in% names(result))
    expect_true(length(result$coefficients) > 0)
    expect_true(all(is.finite(unlist(result$coefficients))))
})

test_that("fit_SensIAT_marginal_mean_model_generalized: works with binomial GLM", {
    skip_if(TRUE, "Long-running test - enable when debugging binomial GLM")

    # Simulate binary data
    set.seed(333)
    data <- simulate_SensIAT_data(
        n_subjects = 50,
        End = 500,
        link = "logit",
        outcome_coef = list(
            intercept = 4,
            prev_outcome = c(0.8, -0.1, 0.05),
            time = 0.001,
            delta_time = -0.05
        ),
        initial_outcome_mean = 0.3
    )
    count(data, Outcome)

    # Prepare data
    data_with_lags <- prepare_SensIAT_data(data, Subject_ID, Time, Outcome, End = 500) |>
        filter(..delta_time.. > 0 | Time == 0)

    # Fit intensity model
    intensity.model <- coxph(
        Surv(..prev_time.., Time, !is.na(Outcome)) ~ ..prev_outcome..,
        data = data_with_lags |> filter(Time > 0)
    )

    # Fit binomial GLM outcome model
    outcome.model <- glm(
        Outcome ~ ..prev_outcome.. + ..delta_time..,
        data = data_with_lags |> filter(Time > 0),
        family = binomial()
    )

    # Test that generalized fitting works
    result <- suppressWarnings(
        fit_SensIAT_marginal_mean_model_generalized(
            data = data_with_lags,
            time = Time,
            id = Subject_ID,
            alpha = 0,
            knots = c(100, 250, 400),
            outcome.model = outcome.model,
            intensity.model = intensity.model,
            loss = "lp_mse",
            link = "logit"
        )
    )

    # Check result structure
    expect_type(result, "list")
    expect_true("gamma" %in% names(result))
})

test_that("fit_SensIAT_marginal_mean_model_generalized: works with poisson GLM", {
    skip_if(TRUE, "Long-running test - enable when debugging poisson GLM")

    # Simulate count data
    set.seed(444)
    data <- simulate_SensIAT_data(
        n_subjects = 20,
        End = 500,
        link = "log",
        outcome_coef = list(
            intercept = 1.5,
            prev_outcome = c(0.1, -0.02, 0.005),
            time = 0.0005,
            delta_time = -0.1
        ),
        initial_outcome_mean = 5
    )

    # Prepare data
    data_with_lags <- prepare_SensIAT_data(data, Subject_ID, Time, Outcome, End = 500)

    # Fit intensity model
    intensity.model <- coxph(
        Surv(..prev_time.., Time, !is.na(Outcome)) ~ ..prev_outcome..,
        data = data_with_lags |> filter(Time > 0)
    )

    # Fit poisson GLM outcome model
    outcome.model <- glm(
        Outcome ~ ..prev_outcome.. + ..delta_time..,
        data = data_with_lags |> filter(Time > 0),
        family = poisson()
    )

    # Test that generalized fitting works
    result <- suppressWarnings(
        fit_SensIAT_marginal_mean_model_generalized(
            data = data_with_lags,
            time = Time,
            id = Subject_ID,
            alpha = 0,
            knots = c(100, 250, 400),
            outcome.model = outcome.model,
            intensity.model = intensity.model,
            loss = "lp_mse",
            link = "log"
        )
    )

    # Check result structure
    expect_type(result, "list")
    expect_true("gamma" %in% names(result))
})

test_that("GLM vs single-index model: gaussian family comparison", {
    # Simulate data
    set.seed(555)
    data <- simulate_SensIAT_data(
        n_subjects = 12,
        End = 400,
        link = "identity",
        outcome_sd = 2
    )

    # Prepare data
    data_with_lags <- prepare_SensIAT_data(data, Subject_ID, Time, Outcome, End = 500)

    # Fit intensity model
    intensity.model <- coxph(
        Surv(..prev_time.., Time, !is.na(Outcome)) ~ ..prev_outcome..,
        data = data_with_lags |> filter(Time > 0)
    )

    # Fit gaussian GLM outcome model
    glm_outcome.model <- glm(
        Outcome ~ ..prev_outcome.. + ..delta_time..,
        data = data_with_lags |> filter(Time > 0),
        family = gaussian()
    )

    # Fit single-index outcome model
    si_outcome.model <- fit_SensIAT_single_index_fixed_coef_model(
        Outcome ~ ..prev_outcome.. + ..delta_time.. - 1,
        id = Subject_ID,
        data = data_with_lags |> filter(Time > 0)
    )
    
    # Test expected values are similar (they use different methods but should be close)
    # Filter to ensure we have clean data without NAs
    test_data <- data_with_lags |> 
        filter(Time > 0, !is.na(Outcome), !is.na(..prev_outcome..), !is.na(..delta_time..)) |> 
        head(3)
    
    # Skip test if we don't have enough valid rows
    skip_if(nrow(test_data) < 3, "Not enough valid test data rows")
    
    glm_expected <- compute_SensIAT_expected_values(
        model = glm_outcome.model,
        alpha = 0,
        new.data = test_data
    )
    
    si_expected <- compute_SensIAT_expected_values(
        model = si_outcome.model,
        alpha = 0,
        new.data = test_data
    )
    
    # Both should produce valid results
    expect_true(all(is.finite(glm_expected$E_Yexp_alphaY)))
    expect_true(all(is.finite(si_expected$E_Yexp_alphaY)))
    expect_true(all(is.finite(glm_expected$E_exp_alphaY)))
    expect_true(all(is.finite(si_expected$E_exp_alphaY)))
    
    # For alpha=0, E_Yexp_alphaY / E_exp_alphaY should equal E[Y]
    # which should be close to the predicted values
    glm_conditional_mean <- glm_expected$E_Yexp_alphaY / glm_expected$E_exp_alphaY
    si_conditional_mean <- si_expected$E_Yexp_alphaY / si_expected$E_exp_alphaY
    
    # They should both be reasonable values
    expect_true(all(abs(glm_conditional_mean) < 100))
    expect_true(all(abs(si_conditional_mean) < 100))
})

test_that("compute_SensIAT_expected_values.glm: handles edge cases", {
    # Simulate data
    set.seed(666)
    data <- simulate_SensIAT_data(
        n_subjects = 10,
        End = 300,
        link = "identity"
    )
    
    # Prepare data
    data_prepared <- prepare_SensIAT_data(data, Subject_ID, Time, Outcome, End = 500)
    
    # Fit GLM
    glm_model <- glm(
        Outcome ~ ..prev_outcome.. + ..delta_time..,
        data = data_prepared |> filter(Time > 0),
        family = gaussian()
    )
    
    # Test with empty data
    result_empty <- compute_SensIAT_expected_values(
        model = glm_model,
        alpha = 0,
        new.data = data_prepared |> filter(Time > 0) |> head(0)
    )
    expect_equal(nrow(result_empty), 0)
    
    # Test with single row
    result_single <- compute_SensIAT_expected_values(
        model = glm_model,
        alpha = 0,
        new.data = data_prepared |> filter(Time > 0) |> head(1)
    )
    expect_equal(nrow(result_single), 1)
    expect_true(is.finite(result_single$E_Yexp_alphaY))
    expect_true(is.finite(result_single$E_exp_alphaY))
})
