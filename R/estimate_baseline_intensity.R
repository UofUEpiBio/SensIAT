globalVariables("cumhaz")

estimate_baseline_intensity <-
    function(intensity.model,
             data = NULL,
             bandwidth = attr(intensity.model, "bandwidth"),
             kernel = attr(intensity.model, "kernel") %||% \(x) 0.75 * (1 - (x)**2) * (abs(x) < 1)) {
        
        # Filter data to remove invalid time intervals if provided
        if (!is.null(data)) {
            # Check for ..prev_time.. and ..time.. columns
            if ("..prev_time.." %in% names(data) && "..time.." %in% names(data)) {
                data <- data |> dplyr::filter(.data$..time.. > .data$..prev_time..)
            }
        }
        
        mf <- if (is.null(data)) {
            model.frame(intensity.model)
        } else {
            model.frame(terms(intensity.model), data = data, na.action = NULL)
        }

        new.time <- mf[[1]][, 2]  # Extract stop time from Surv object
        
        # Handle both stratified and non-stratified models
        strata_col_idx <- attr(terms(intensity.model), "specials")$strata
        if (!is.null(strata_col_idx) && length(strata_col_idx) > 0) {
            new.strata <- mf[[strata_col_idx]] |> as.character()
        } else {
            # No stratification - create a single dummy stratum
            new.strata <- rep("all", length(new.time))
        }


        ############  Estimated baseline intensities for each stratum ############
        data_surv <- survfit(intensity.model, newdata = construct_baseline_df(intensity.model))
        strata <- data_surv$strata
        
        # Handle case where survfit doesn't return strata (unstratified model)
        if (is.null(strata)) {
            strata <- setNames(length(data_surv$time), "all")
        }

        # the following is to find the baseline intensity function for each strata k
        cumhaz.data <-
            tibble(
                time = data_surv$time,
                cumhaz = data_surv$cumhaz,
                strata = rep(names(strata), strata)
            ) |>
            group_by(strata) |>
            mutate(hazard = cumhaz - dplyr::lag(cumhaz, default = 0, order_by = time))

        if (is.null(bandwidth)) {
            bandwidth <- dpill(cumhaz.data$time, cumhaz.data$hazard)
        }

        time.diff <- outer(cumhaz.data$time, new.time, `-`)
        strata.equal <- outer(cumhaz.data$strata, new.strata, `==`)

        list(
            baseline_intensity = colSums(kernel(time.diff / bandwidth) * strata.equal * cumhaz.data$hazard) / bandwidth,
            bandwidth = bandwidth,
            kernel = kernel
        )
    }

construct_baseline_df <- function(intensity.model) {
    df <- model.frame(intensity.model)[1, , drop = FALSE]
    mm <- model.matrix(intensity.model, data = df)
    # Exclude Surv columns from mutation
    surv_cols <- sapply(df, function(x) inherits(x, "Surv"))
    cols_to_zero <- setdiff(colnames(mm), names(df)[surv_cols])
    mutate(df, across(any_of(cols_to_zero), ~0))
}
