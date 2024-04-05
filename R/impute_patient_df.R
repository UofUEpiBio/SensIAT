impute_patient_df <- function(eval.times, df_i, variables, centering, right=TRUE){
    outcomes <- pull(df_i, !!(variables$outcome))
    orig.time <- pull(df_i, !!variables$time)

    time_mean <- centering[[1]]
    time_sd   <- centering[[2]]
    Δ_time_mean <- centering[[3]]
    Δ_time_sd   <- centering[[4]]

    period <- as.numeric(cut(eval.times, c(orig.time, Inf), right = right))
    delta_time = eval.times - orig.time[period]
    norm_time  = (eval.times - time_mean)/time_sd
    norm_delta_time = (delta_time - Δ_time_mean)/Δ_time_sd
    prev_outcome = (!!outcomes)[period]


    tibble(
        time = eval.times,
        period,
        delta_time,
        norm_time,
        norm_delta_time,
        prev_outcome,
        outcome = 0
    ) |>
        dplyr::rename(
            any_of(
                rlang::set_names(
                    names(variables),
                    sapply(variables, deparse)
                )
            )
        )
}
