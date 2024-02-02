impute_patient_df <- function(eval.times, df_i, object, right=TRUE){
    outcomes <- pull(df_i, !!(object$variables$outcome))

    time_mean <- object$outcome_model_centering[[1]]
    time_sd   <- object$outcome_model_centering[[2]]
    Δ_time_mean <- object$outcome_model_centering[[3]]
    Δ_time_sd   <- object$outcome_model_centering[[4]]

    tibble(
        time = eval.times,
        period = as.numeric(cut(eval.times, c(pull(df_i, !!object$variables$time), Inf), right = right)),
    ) |>
        mutate(
            delta_time := time - df_i$time[period],
            norm_time = (time - time_mean)/time_sd,
            norm_delta_time = (delta_time - Δ_time_mean)/Δ_time_sd,
            prev_outcome = (!!outcomes)[period],
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
