test_that("extrapolate_from_last_observation works with basic linear extrapolation", {
    df <- data.frame(
        time = c(0, 1, 2, 3),
        x1 = c(10, 12, 14, 16),
        x2 = c(5, 5.5, 6, 6.5),
        id = 1
    )

    result <- extrapolate_from_last_observation(
        target_time = 5,
        data = df,
        time_var = "time",
        slopes = c(x1 = 2, x2 = 0.5)
    )

    # Check extrapolated values
    # x1: 16 + (5-3)*2 = 20
    expect_equal(result$x1, 20)
    # x2: 6.5 + (5-3)*0.5 = 7.5
    expect_equal(result$x2, 7.5)
    # time should be updated
    expect_equal(result$time, 5)
    # id should be carried forward
    expect_equal(result$id, 1)

    # Check attributes
    expect_equal(attr(result, "source_period"), 4)
    expect_equal(attr(result, "source_time"), 3)
    expect_equal(attr(result, "delta_time"), 2)
})

test_that("extrapolate_from_last_observation works with partial slopes", {
    df <- data.frame(
        time = c(0, 1, 2, 3),
        x1 = c(10, 12, 14, 16),
        x2 = c(5, 5.5, 6, 6.5),
        id = 1
    )

    result <- extrapolate_from_last_observation(
        target_time = 4,
        data = df,
        time_var = "time",
        slopes = c(x1 = 2) # Only x1 has a slope
    )

    # x1 should be extrapolated: 16 + (4-3)*2 = 18
    expect_equal(result$x1, 18)
    # x2 should be carried forward unchanged
    expect_equal(result$x2, 6.5)
    # id should be carried forward
    expect_equal(result$id, 1)
    expect_equal(result$time, 4)
})

test_that("extrapolate_from_last_observation works with no slopes (carry forward)", {
    df <- data.frame(
        time = c(0, 1, 2, 3),
        x1 = c(10, 12, 14, 16),
        x2 = c(5, 5.5, 6, 6.5),
        id = 1
    )

    result <- extrapolate_from_last_observation(
        target_time = 4,
        data = df,
        time_var = "time",
        slopes = NULL
    )

    # All values should be carried forward from time = 3
    expect_equal(result$x1, 16)
    expect_equal(result$x2, 6.5)
    expect_equal(result$id, 1)
    expect_equal(result$time, 4)
})

test_that("extrapolate_from_last_observation matches influence_sim.R pattern", {
    # Simulate the influence_sim.R scenario (line 202)
    # xi <- individual_X[period, ,drop=FALSE] + (time - times[period]) * x_slope
    individual_X <- matrix(c(10, 5, 12, 5.5, 14, 6, 16, 6.5), ncol = 2, byrow = TRUE)
    times <- c(0, 1, 2, 3)
    x_slope <- c(2, 0.5)

    # Manual calculation as in influence_sim.R
    target_time <- 5
    period <- max(which(times <= target_time))
    xi_manual <- individual_X[period, , drop = FALSE] + (target_time - times[period]) * x_slope

    # Using our function
    df_sim <- data.frame(
        time = times,
        V1 = individual_X[, 1],
        V2 = individual_X[, 2]
    )
    result_sim <- extrapolate_from_last_observation(
        target_time = 5,
        data = df_sim,
        time_var = "time",
        slopes = c(V1 = x_slope[1], V2 = x_slope[2])
    )

    # Values should match exactly
    expect_equal(result_sim$V1, xi_manual[1, 1])
    expect_equal(result_sim$V2, xi_manual[1, 2])
})

test_that("extrapolate_from_last_observation_multiple works with multiple time points", {
    df <- data.frame(
        time = c(0, 1, 2, 3),
        x1 = c(10, 12, 14, 16),
        x2 = c(5, 5.5, 6, 6.5)
    )

    result <- extrapolate_from_last_observation_multiple(
        target_times = c(3.5, 4, 5),
        data = df,
        time_var = "time",
        slopes = c(x1 = 2, x2 = 0.5)
    )

    expect_equal(nrow(result), 3)

    # Check time = 3.5
    expect_equal(result$time[1], 3.5)
    expect_equal(result$x1[1], 16 + (3.5 - 3) * 2) # 17
    expect_equal(result$x2[1], 6.5 + (3.5 - 3) * 0.5) # 6.75

    # Check time = 4
    expect_equal(result$time[2], 4)
    expect_equal(result$x1[2], 16 + (4 - 3) * 2) # 18
    expect_equal(result$x2[2], 6.5 + (4 - 3) * 0.5) # 7.0

    # Check time = 5
    expect_equal(result$time[3], 5)
    expect_equal(result$x1[3], 16 + (5 - 3) * 2) # 20
    expect_equal(result$x2[3], 6.5 + (5 - 3) * 0.5) # 7.5
})

test_that("extrapolate_from_last_observation validates inputs correctly", {
    df <- data.frame(
        time = c(0, 1, 2, 3),
        x1 = c(10, 12, 14, 16)
    )

    # target_time must be scalar
    expect_error(
        extrapolate_from_last_observation(
            target_time = c(4, 5),
            data = df,
            time_var = "time",
            slopes = c(x1 = 2)
        ),
        "must be a single numeric value"
    )

    # data must be a data.frame
    expect_error(
        extrapolate_from_last_observation(
            target_time = 4,
            data = list(time = c(0, 1, 2, 3)),
            time_var = "time",
            slopes = c(x1 = 2)
        ),
        "must be a data.frame"
    )

    # time_var must be in data
    expect_error(
        extrapolate_from_last_observation(
            target_time = 4,
            data = df,
            time_var = "nonexistent",
            slopes = c(x1 = 2)
        ),
        "must be a character string naming a column"
    )

    # slopes must be named
    expect_error(
        extrapolate_from_last_observation(
            target_time = 4,
            data = df,
            time_var = "time",
            slopes = c(2) # Unnamed
        ),
        "must be a named numeric vector"
    )

    # slope variables must exist in data
    expect_error(
        extrapolate_from_last_observation(
            target_time = 4,
            data = df,
            time_var = "time",
            slopes = c(x1 = 2, nonexistent = 1)
        ),
        "not found in `data`"
    )
})

test_that("extrapolate_from_last_observation handles edge cases", {
    df <- data.frame(
        time = c(0, 1, 2, 3),
        x1 = c(10, 12, 14, 16)
    )

    # Extrapolate to exact time point (no extrapolation needed)
    result <- extrapolate_from_last_observation(
        target_time = 3,
        data = df,
        time_var = "time",
        slopes = c(x1 = 2)
    )
    expect_equal(result$x1, 16)
    expect_equal(result$time, 3)
    expect_equal(attr(result, "delta_time"), 0)

    # Target time before all observations (with strict = FALSE)
    result_early <- extrapolate_from_last_observation(
        target_time = -1,
        data = df,
        time_var = "time",
        slopes = c(x1 = 2),
        strict = FALSE
    )
    # Should use first observation
    expect_equal(result_early$x1, 10 + (-1 - 0) * 2) # 8
    expect_equal(attr(result_early, "source_period"), 1)

    # Target time before all observations (with strict = TRUE) should error
    expect_error(
        extrapolate_from_last_observation(
            target_time = -1,
            data = df,
            time_var = "time",
            slopes = c(x1 = 2),
            strict = TRUE
        ),
        "No observations found at or before"
    )
})

test_that("extrapolate_from_last_observation handles single-row data", {
    df <- data.frame(
        time = 1,
        x1 = 10,
        x2 = 5
    )

    result <- extrapolate_from_last_observation(
        target_time = 3,
        data = df,
        time_var = "time",
        slopes = c(x1 = 2, x2 = 0.5)
    )

    expect_equal(result$x1, 10 + (3 - 1) * 2) # 14
    expect_equal(result$x2, 5 + (3 - 1) * 0.5) # 6
    expect_equal(result$time, 3)
})

test_that("extrapolate_from_last_observation handles zero slopes", {
    df <- data.frame(
        time = c(0, 1, 2, 3),
        x1 = c(10, 12, 14, 16),
        x2 = c(5, 5.5, 6, 6.5)
    )

    result <- extrapolate_from_last_observation(
        target_time = 5,
        data = df,
        time_var = "time",
        slopes = c(x1 = 0, x2 = 0) # Zero slopes = carry forward
    )

    expect_equal(result$x1, 16)
    expect_equal(result$x2, 6.5)
    expect_equal(result$time, 5)
})

test_that("extrapolate_from_last_observation handles negative slopes", {
    df <- data.frame(
        time = c(0, 1, 2, 3),
        x1 = c(10, 12, 14, 16),
        x2 = c(10, 9, 8, 7)
    )

    result <- extrapolate_from_last_observation(
        target_time = 5,
        data = df,
        time_var = "time",
        slopes = c(x1 = 2, x2 = -1) # x2 decreasing
    )

    expect_equal(result$x1, 16 + (5 - 3) * 2) # 20
    expect_equal(result$x2, 7 + (5 - 3) * (-1)) # 5
    expect_equal(result$time, 5)
})
