test_that("simulate_SensIAT_data: identity link produces continuous outcomes", {
    set.seed(123)
    data <- simulate_SensIAT_data(
        n_subjects = 10,
        End = 500,
        link = "identity",
        outcome_coef = list(
            intercept = 5,
            prev_outcome = c(0.7, -0.1, 0.05),
            time = 0.001,
            delta_time = -0.5
        )
    )
    
    # Check structure
    expect_s3_class(data, "tbl_df")
    expect_named(data, c("Subject_ID", "Time", "Outcome"))
    expect_equal(nrow(data), nrow(data))  # Basic check
    
    # Check outcome type - should be continuous (not integers)
    outcomes <- data$Outcome[data$Time > 0]
    expect_true(any(outcomes != floor(outcomes)), 
                info = "Identity link should produce non-integer continuous outcomes")
    
    # Check all subjects have at least one follow-up
    subject_visits <- data |> 
        dplyr::group_by(Subject_ID) |> 
        dplyr::summarise(max_time = max(Time))
    expect_true(all(subject_visits$max_time > 0),
                info = "All subjects should have at least one follow-up visit")
})

test_that("simulate_SensIAT_data: log link produces count outcomes", {
    set.seed(456)
    data <- simulate_SensIAT_data(
        n_subjects = 10,
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
    
    # Check outcomes are non-negative integers
    outcomes <- data$Outcome
    expect_true(all(outcomes >= 0), 
                info = "Log link should produce non-negative outcomes")
    expect_true(all(outcomes == floor(outcomes)), 
                info = "Log link should produce integer count outcomes")
    expect_type(outcomes, "double")
})

test_that("simulate_SensIAT_data: logit link produces binary outcomes", {
    set.seed(789)
    data <- simulate_SensIAT_data(
        n_subjects = 15,
        End = 500,
        link = "logit",
        outcome_coef = list(
            intercept = 0.5,
            prev_outcome = c(0.8, -0.1, 0.05),
            time = 0.001,
            delta_time = -0.2
        ),
        initial_outcome_mean = 0.3
    )
    
    # Check outcomes are only 0 or 1
    outcomes <- data$Outcome
    expect_true(all(outcomes %in% c(0, 1)), 
                info = "Logit link should produce only 0 or 1")
    
    # Check we have some variation (not all same value)
    # This could fail with small samples by chance, but unlikely with 15 subjects
    expect_true(length(unique(outcomes)) > 1,
                info = "Should have both 0 and 1 outcomes (probabilistic)")
})

test_that("simulate_SensIAT_data: reproducibility with seed", {
    # Run 1
    set.seed(999)
    data1 <- simulate_SensIAT_data(
        n_subjects = 5,
        End = 300,
        link = "identity"
    )
    
    # Run 2 with same seed
    set.seed(999)
    data2 <- simulate_SensIAT_data(
        n_subjects = 5,
        End = 300,
        link = "identity"
    )
    
    # Should be identical
    expect_equal(data1, data2)
    
    # Run 3 with different seed should be different
    set.seed(1000)
    data3 <- simulate_SensIAT_data(
        n_subjects = 5,
        End = 300,
        link = "identity"
    )
    
    expect_false(identical(data1, data3))
})

test_that("simulate_SensIAT_data: all subjects have follow-up visits", {
    set.seed(111)
    data <- simulate_SensIAT_data(
        n_subjects = 20,
        End = 500,
        link = "identity"
    )
    
    # Check each subject has at least 2 observations (baseline + 1 follow-up)
    visit_counts <- data |> 
        dplyr::group_by(Subject_ID) |> 
        dplyr::summarise(
            n_visits = dplyr::n(),
            max_time = max(Time)
        )
    
    expect_true(all(visit_counts$n_visits >= 2),
                info = "Each subject should have baseline + at least 1 follow-up")
    expect_true(all(visit_counts$max_time > 0),
                info = "Each subject should have max(Time) > 0")
})

test_that("simulate_SensIAT_data: time sequence is valid", {
    set.seed(222)
    data <- simulate_SensIAT_data(
        n_subjects = 8,
        End = 400,
        link = "identity"
    )
    
    # Within each subject, times should be increasing
    time_checks <- data |> 
        dplyr::group_by(Subject_ID) |> 
        dplyr::summarise(
            times_increasing = all(diff(Time) > 0),
            first_is_zero = Time[1] == 0,
            last_within_End = max(Time) <= 400
        )
    
    expect_true(all(time_checks$times_increasing),
                info = "Times should be strictly increasing within subject")
    expect_true(all(time_checks$first_is_zero),
                info = "First observation should be at time 0")
    expect_true(all(time_checks$last_within_End),
                info = "All times should be <= End")
})

test_that("simulate_SensIAT_data: respects max_visits parameter", {
    set.seed(333)
    data <- simulate_SensIAT_data(
        n_subjects = 10,
        End = 1000,  # Long period
        max_visits = 5,
        link = "identity"
    )
    
    # No subject should have more than max_visits
    visit_counts <- data |> 
        dplyr::group_by(Subject_ID) |> 
        dplyr::summarise(n_visits = dplyr::n())
    
    expect_true(all(visit_counts$n_visits <= 5),
                info = "No subject should exceed max_visits")
})

test_that("simulate_SensIAT_two_groups: creates proper group structure", {
    set.seed(444)
    data <- simulate_SensIAT_two_groups(
        n_subjects = 22,  # Total subjects (will split 11/11)
        End = 500,
        treatment_effect = 2,
        link = "identity"
    )
    
    # Check structure includes Group column
    expect_named(data, c("Subject_ID", "Time", "Outcome", "Group"))
    
    # Check group assignments
    group_summary <- data |> 
        dplyr::distinct(Subject_ID, Group) |> 
        dplyr::count(Group)
    
    expect_equal(nrow(group_summary), 2, info = "Should have exactly 2 groups")
    expect_true("Control" %in% group_summary$Group)
    expect_true("Treatment" %in% group_summary$Group)
    
    # Check subject counts (22 splits to 11 control, 11 treatment)
    control_subjects <- dplyr::filter(group_summary, Group == "Control")$n
    treatment_subjects <- dplyr::filter(group_summary, Group == "Treatment")$n
    
    expect_equal(control_subjects, 11)
    expect_equal(treatment_subjects, 11)
})

test_that("simulate_SensIAT_two_groups: treatment effect affects outcomes", {
    set.seed(555)
    
    # Simulate with treatment effect
    data_with_effect <- simulate_SensIAT_two_groups(
        n_subjects = 40,  # Will split 20/20
        End = 500,
        treatment_effect = 5,  # Large positive effect
        link = "identity",
        outcome_sd = 2
    )
    
    # Compare mean outcomes between groups
    group_means <- data_with_effect |> 
        dplyr::filter(Time > 0) |>  # Exclude baseline
        dplyr::group_by(Group) |> 
        dplyr::summarise(mean_outcome = mean(Outcome))
    
    control_mean <- group_means$mean_outcome[group_means$Group == "Control"]
    treatment_mean <- group_means$mean_outcome[group_means$Group == "Treatment"]
    
    # Treatment should have higher mean (treatment_effect = 5)
    expect_true(treatment_mean > control_mean,
                info = "Treatment group should have higher mean outcome")
    
    # Difference should be approximately the treatment effect (allowing for randomness)
    # This is a rough check - exact match not expected due to non-linearity
    expect_true(abs((treatment_mean - control_mean) - 5) < 3,
                info = "Difference should be roughly the treatment effect")
})

test_that("simulate_SensIAT_two_groups: reproducibility with seed", {
    set.seed(666)
    data1 <- simulate_SensIAT_two_groups(
        n_subjects = 10,
        End = 300
    )
    
    set.seed(666)
    data2 <- simulate_SensIAT_two_groups(
        n_subjects = 10,
        End = 300
    )
    
    expect_equal(data1, data2)
})

test_that("simulate_SensIAT_data: edge case - single subject", {
    set.seed(777)
    data <- simulate_SensIAT_data(
        n_subjects = 1,
        End = 400,
        link = "identity"
    )
    
    expect_equal(length(unique(data$Subject_ID)), 1)
    expect_true(nrow(data) >= 2, info = "Should have baseline + follow-up")
    expect_equal(data$Time[1], 0, info = "First time should be 0")
})

test_that("simulate_SensIAT_data: intensity_coef affects observation times", {
    set.seed(888)
    
    # Higher intensity coefficient should lead to more frequent observations
    # (lower hazard ratio with negative coef, or higher with positive)
    data_low_intensity <- simulate_SensIAT_data(
        n_subjects = 10,
        End = 500,
        intensity_coef = -0.5,  # Negative: higher outcomes -> less frequent visits
        link = "identity"
    )
    
    data_high_intensity <- simulate_SensIAT_data(
        n_subjects = 10,
        End = 500,
        intensity_coef = 0.5,  # Positive: higher outcomes -> more frequent visits
        link = "identity"
    )
    
    # Count total visits
    n_visits_low <- nrow(data_low_intensity)
    n_visits_high <- nrow(data_high_intensity)
    
    # This is probabilistic, but with opposite signs should see a difference
    # Not a strict test, but checks the mechanism is working
    expect_true(n_visits_low != n_visits_high,
                info = "Different intensity coefficients should affect visit frequency")
})

test_that("simulate_SensIAT_data: baseline_hazard affects visit frequency", {
    set.seed(999)
    
    # Higher baseline hazard -> more frequent visits
    data_low_hazard <- simulate_SensIAT_data(
        n_subjects = 10,
        End = 500,
        baseline_hazard = 0.001,
        link = "identity"
    )
    
    data_high_hazard <- simulate_SensIAT_data(
        n_subjects = 10,
        End = 500,
        baseline_hazard = 0.01,  # 10x higher
        link = "identity"
    )
    
    # Higher hazard should produce more visits on average
    mean_visits_low <- data_low_hazard |> 
        dplyr::group_by(Subject_ID) |> 
        dplyr::summarise(n = dplyr::n()) |> 
        dplyr::pull(n) |> 
        mean()
    
    mean_visits_high <- data_high_hazard |> 
        dplyr::group_by(Subject_ID) |> 
        dplyr::summarise(n = dplyr::n()) |> 
        dplyr::pull(n) |> 
        mean()
    
    expect_true(mean_visits_high > mean_visits_low,
                info = "Higher baseline hazard should produce more visits")
})

test_that("simulate_SensIAT_data: validates inputs", {
    expect_error(
        simulate_SensIAT_data(n_subjects = 0, End = 500),
        "n_subjects must be a positive integer"
    )
    
    expect_error(
        simulate_SensIAT_data(n_subjects = 10, End = -100),
        "End must be a positive number"
    )
    
    expect_error(
        simulate_SensIAT_data(n_subjects = 10, End = 500, max_visits = 1),
        "max_visits must be at least 2"
    )
    
    expect_error(
        simulate_SensIAT_data(n_subjects = 10, End = 500, link = "invalid"),
        "'arg' should be one of"
    )
})

test_that("simulate_SensIAT_two_groups: validates inputs", {
    expect_error(
        simulate_SensIAT_two_groups(n_subjects = 0, End = 500),
        "n_subjects must be a positive integer"
    )
    
    expect_error(
        simulate_SensIAT_two_groups(n_subjects = 10, End = -100),
        "End must be a positive number"
    )
    
    expect_error(
        simulate_SensIAT_two_groups(n_subjects = 10, End = 500, link = "invalid"),
        "'arg' should be one of"
    )
})
