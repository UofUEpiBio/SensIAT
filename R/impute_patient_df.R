impute_patient_df <- function(eval.times, df_i, object, right=TRUE){
    outcomes <- pull(df_i, !!(object$variables$outcome))
    orig.time <- pull(df_i, !!object$variables$time)

    time_mean <- object$outcome_model_centering[[1]]
    time_sd   <- object$outcome_model_centering[[2]]
    Δ_time_mean <- object$outcome_model_centering[[3]]
    Δ_time_sd   <- object$outcome_model_centering[[4]]

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
                    names(object$variables),
                    sapply(object$variables, deparse)
                )
            )
        )
}
