numerically_integrate_influence_term_2_for_one_alpha_and_one_patient <-
function(
    df_i,
    expected_value,
    base,
    variables,
    resolution = 1000,
    a = min(base@knots),
    b = max(base@knots),
    ...
){
    eval.times <- seq(a, b, length.out = resolution)
    B <- evaluate(base, eval.times)
    spline_df_est <-
        bind_rows(
            impute_patient_df(head(eval.times,  1), df_i, variables=variables, ..., right = FALSE),
            impute_patient_df(tail(eval.times, -1), df_i, variables=variables, ..., right = TRUE)
        )


    Ey = expected_value(spline_df_est)
    dt = diff(pull(spline_df_est, variables$time))

    crossprod(head(B, -1), head(Ey, -1) * dt) / 2 +
    crossprod(tail(B, -1), tail(Ey, -1) * dt) / 2
}


if(F){



}

numerically_integrate_influence_term_2_for_one_alpha_and_one_patient_piecewise <-
function(
    df_i, #< complete patient data
    expected_value,
    base,
    variables,
    ..., #< passed along
    resolution.within.period = 50
){
    df_i |>
        transmute(
            a = pmax(!!variables$prev_time, !!min(base@knots)),
            b = pmin(!!variables$time     , !!max(base@knots))
        ) |>
        filter(
            !!min(base@knots) <= a, a<b, b <= !!max(base@knots)
        ) |>
        pmap(
            numerically_integrate_influence_term_2_for_one_alpha_and_one_patient,
            df_i = df_i,
            expected_value = expected_value,
            base = base,
            variables = variables,
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

    library(ARCdata)
    fitted.trt.sim <-
        fit_PCORI_within_group_model(
            group.data = filter(ARC_data, Trt=='home_visits'),
            outcome_modeler = PCORI_sim_outcome_modeler,
            id.var = elig_pid,
            outcome.var = Asthma_control,
            time.var = time,
            End = 830
        )
    time.pw <- system.time({
        pred.pw <- predict(fitted.trt.sim, time = c(90, 180),
                              alpha = c(0),
                              intensity.bandwidth = 30,
                              knots=c(60,60,60,60,260,460,460,460,460),
                              integration.method = 'piecewise'
        )
    })
}
