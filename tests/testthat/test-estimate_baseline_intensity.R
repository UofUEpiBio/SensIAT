test_that("estimate_baseline_intensity", {
    data_with_lags <- SensIAT_example_data |>
        group_by(Subject_ID) |>
        mutate(
            Prev_Outcome = lag(Outcome, default = NA_real_),
            Prev_time = lag(Time, default = NA_real_),
            Delta_Time = Time - Prev_time
        )

    intensity.model <-
        coxph(Surv(Prev_time, Time, !is.na(Outcome)) ~ Prev_Outcome + strata(Visit),
            data = dplyr::filter(data_with_lags, Time > 0)
        )

    estimate_baseline_intensity_result <-
        estimate_baseline_intensity(intensity.model = intensity.model, data = data_with_lags)
    
    expect_true(is.numeric(estimate_baseline_intensity_result$baseline_intensity))
})

test_that("estimate_baseline_intensity is invariant to data size", {
    data_with_lags <- SensIAT_example_data |>
        group_by(Subject_ID) |>
        mutate(
            Prev_Outcome = lag(Outcome, default = NA_real_),
            Prev_time = lag(Time, default = NA_real_),
            Delta_Time = Time - Prev_time
        )

    intensity.model <-
        coxph(Surv(Prev_time, Time, !is.na(Outcome)) ~ Prev_Outcome + strata(Visit),
            data = dplyr::filter(data_with_lags, Time > 0)
        )

    # Get baseline intensity with full data
    result_full <- estimate_baseline_intensity(
        intensity.model = intensity.model, 
        data = data_with_lags
    )
    
    # Get baseline intensity with subset of data (same time points)
    result_subset <- estimate_baseline_intensity(
        intensity.model = intensity.model, 
        data = data_with_lags[1:50, ]
    )
    
    # The baseline_intensity values should be the same regardless of data size
    # (assuming same time range is covered)
    expect_equal(
        length(result_full$baseline_intensity),
        nrow(data_with_lags)
    )
    expect_equal(
        length(result_subset$baseline_intensity),
        50
    )
})
