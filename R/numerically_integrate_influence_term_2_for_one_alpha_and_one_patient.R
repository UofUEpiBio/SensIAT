numerically_integrate_influence_term_2_for_one_alpha_and_one_patient <-
function(
    df_i,
    object,
    alpha,
    base,
    resolution = 1000,
    a = min(base@knots),
    b = max(base@knots),
    ...
){
    expected_value <- \(data, ...){
        matrix(
            pull(pcori_conditional_means(object$outcome_model, alpha=alpha, new.data = data), 'E_Y_past'),
            nrow = nrow(data)
        )
    }
    eval.times <- seq(a, b, length.out = resolution)
    B <- evaluate(orthogonalize(base), eval.times)
    spline_df_est <-
        bind_rows(
            impute_patient_df(head(eval.times, 1), df_i, object, FALSE),
            impute_patient_df(tail(eval.times, -1), df_i, object)
        )


    Ey = pull(pcori_conditional_means(object$outcome_model, alpha=alpha, new.data = spline_df_est), 'E_Y_past') |>
        matrix(nrow = resolution, ncol = length(alpha))
    dt = diff(pull(spline_df_est, object$variables$time))

    crossprod(head(B, -1), head(Ey, -1) * dt) / 2 +
    crossprod(tail(B, -1), tail(Ey, -1) * dt) / 2
}


if(F){



}

numerically_integrate_influence_term_2_for_one_alpha_and_one_patient_piecewise <-
function(
    df_i, #< complete patient data
    ..., #< passed along
    resolution.within.period = 50
){
    df_i |>
        filter(
            ( !!min(base@knots) <= !!object$variables$time
            & !!object$variables$time <= !!max(base@knots)
            ) | ( !!min(base@knots) <= !!object$variables$prev_time
            & !!object$variables$prev_time <= !!max(base@knots)
            )
        ) |>
        transmute(
            a = pmax(!!object$variables$prev_time, !!min(base@knots)),
            b = pmin(!!object$variables$time     , !!max(base@knots))
        ) |>
        pmap(
            numerically_integrate_influence_term_2_for_one_alpha_and_one_patient,
            df_i = df_i,
            ...,
            resolution = resolution.within.period
        ) |> reduce(`+`)
}
if(F){

    numerically_integrate_influence_term_2_for_one_alpha_and_one_patient_piecewise(
        df_i,
        object = object,
        alpha = alpha,
        base = base
    )
}
