test_that("multiplication works", {
    data_with_lags <- SensIAT_example_data |>
        dplyr::group_by(Subject_ID) |>
        dplyr::mutate(
            ..prev_outcome.. = dplyr::lag(Outcome, default = NA_real_, order_by = Time),
            ..prev_time.. = dplyr::lag(Time, default = 0, order_by = Time),
            ..delta_time.. = Time - dplyr::lag(.data$Time, default = NA_real_, order_by = Time)
        )
    # Create the observation time intensity model
    intensity.model <-
        survival::coxph(survival::Surv(..prev_time.., Time, !is.na(Outcome)) ~ ..prev_outcome.. + strata(Visit),
            data = data_with_lags |> dplyr::filter(.data$Time > 0)
        )
    # Create the observed outcome model
    outcome.model <-
        fit_SensIAT_single_index_fixed_coef_model(
            Outcome ~ splines::ns(..prev_outcome.., df = 3) + ..delta_time.. - 1,
            id = Subject_ID,
            data = data_with_lags |> dplyr::filter(Time > 0)
        )
    fit_SensIAT_marginal_mean_model_generalized(
        data = data_with_lags,
        time = data_with_lags$Time,
        id = data_with_lags$Subject_ID,
        alpha = 0,
        knots = c(60, 260, 460),
        outcome.model = outcome.model,
        intensity.model = intensity.model,
        loss = "lp_mse",
        link = "log",
        impute_data = \(t, df){
            data_wl <- df |>
                mutate(
                    ..prev_time.. = Time,
                    ..prev_outcome.. = Outcome,
                    ..delta_time.. = 0
                )
            extrapolate_from_last_observation(t, data_wl, "Time", slopes = c("..delta_time.." = 1))
        },
        BBsolve.control = list(maxit = 10, tol = 1e-4),
        term2_method = "fast"
    )
    time <- data_with_lags$Time
    id <- data_with_lags$Subject_ID
})
