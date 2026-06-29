test_that("make_single_index_simulator requires outcome_model input", {
    expect_error(
        make_single_index_simulator(NULL),
        "outcome_model must be a fitted Single-index outcome model"
    )
    
    expect_error(
        make_single_index_simulator("not a model"),
        "outcome_model must be a fitted Single-index outcome model"
    )
    
    # Used from lm documentation
    ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
    trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
    group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
    weight <- c(ctl, trt)
    lm.D9 <- lm(weight ~ group)
    expect_error(
        make_single_index_simulator(lm.D9),
        "outcome_model must be a fitted Single-index outcome model"
    )
})

test_that("make_single_index_simulator returns a function with correct class", {
    # Create a minimal outcome model for testing
    set.seed(123)

    data_test <- simulate_SensIAT_data(
        n_subjects = 100,
        End = 830,
        intensity_coef = -0.05,
        outcome_coef = list(
            prev_outcome = c(0.7, -0.1, 0.05),
            time = -0.001,
            delta_time = -0.002,
            intercept = 2
        ),
        baseline_hazard = 0.005,
        outcome_sd = 1.5
    ) |>
    prepare_SensIAT_data(
        id.var = "Subject_ID",
        outcome.var = "Outcome",
        time.var = "Time",
        End = 830
    ) |> rename_with( ~ gsub("(^\\.+|\\.+$)", "", .x))

    outcome_model <- fit_SensIAT_single_index_fixed_coef_model(
        formula = Outcome ~ prev_outcome + Time + delta_time,
        data = data_test,
        id = Subject_ID,
        initial = c(0, 0, 0, 0)
    )

    sim_func <- make_single_index_simulator(outcome_model)
    
    # Check class and that it's callable
    expect_s3_class(sim_func, "SensIAT_single_index_simulator")
    expect_true(is.function(sim_func))
    
    # Check that metadata is attached
    expect_identical(attr(sim_func, "outcome_model"), outcome_model)
    
    # Test that simulator returns a single numeric value
    sampled_1 <- sim_func(newdata = tibble(prev_outcome = 5, Time = 100, delta_time = 90))
    expect_length(sampled_1, 1)
    expect_true(is.numeric(sampled_1))
    expect_true(is.finite(sampled_1))
    
    # Test multiple calls produce different values (probabilistic)
    set.seed(456)
    samples <- replicate(10, sim_func(newdata = tibble(prev_outcome = 5, Time = 100, delta_time = 90)))
    expect_true(length(unique(samples)) > 1)  # Should have variation (high probability)

    # Non-data.frame newdata
    expect_error(
        sim_func(newdata = c(1, 2, 3)),
        "newdata must be a data.frame"
    )
    
    # Empty data.frame newdata
    expect_error(
        sim_func(newdata = data.frame()),
        "newdata must be a data.frame with at least one row"
    )
    
    # Missing 'time' column in newdata
    newdata_missing <- data.frame(
        prev_outcome = 5,
        delta_time = 90
    )
    
    expect_error(
        sim_func(newdata = newdata_missing),
        "Failed to construct model matrix"
    )
    
    outcome_model_no_var <- fit_SensIAT_single_index_fixed_coef_model(
        formula = Outcome ~ prev_outcome + Time + delta_time,
        data = data_test |> mutate(Outcome = 5),  # No variation in outcome
        id = Subject_ID,
        initial = c(0, 0, 0, 0)
    )

    sim_func <- make_single_index_simulator(outcome_model)

    # Should fall back to mean when y_seq has limited variation
    sampled <- sim_func(prev_outcome = 5, time = 100, delta_time = 90)
    expect_length(sampled, 1)
    expect_true(is.numeric(sampled))
    expect_true(is.finite(sampled))
})

test_that("simulator output is consistent with model support", {
    set.seed(123)
    data_test <- data.frame(
        outcome = rnorm(20, 5, 1),
        prev_outcome = rnorm(20, 5, 1),
        time = runif(20, 0, 365),
        delta_time = rexp(20, 0.01)
    )
    
    outcome_model <- fit_SensIAT_single_index_fixed_coef_model(
        formula = outcome ~ prev_outcome + time + delta_time,
        data = data_test,
        id = NULL
    )
    
    sim_func <- make_single_index_simulator(outcome_model)
    
    # Sample many times and check they're in the support of observed outcomes
    set.seed(789)
    samples <- replicate(50, sim_func(prev_outcome = 5, time = 100, delta_time = 90))
    
    # Extract unique values from training data
    y_obs <- model.response(model.frame(outcome_model))
    y_unique <- unique(y_obs)
    
    # Each sample should be close to one of the observed values (or mean)
    expected_values <- c(y_unique, mean(y_obs, na.rm = TRUE))
    
    # Allow small numerical tolerance
    for (samp in samples) {
        expect_true(
            any(abs(samp - expected_values) < 1e-10),
            info = sprintf("Sample %f not in expected support", samp)
        )
    }
})
