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
        data = dplyr::filter(data_with_lags, Time  > 0))

    estimate_baseline_intensity_result <-
        estimate_baseline_intensity(intensity.model = intensity.model, data = data_with_lags)


})
