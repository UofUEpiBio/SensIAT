compute_influence_for_one_alpha_and_one_patient <-
function(
    df_i,
    patient.df,
    alpha,
    variables,
    intensity.model,
    outcome.model,
    base,
    control,
    ...
){
    if (getOption('PCORI::do_arg_checks', TRUE))
        assert_that(
            rlang::is_atomic(alpha), is.numeric(alpha),
            # is(object, 'PCORI_within_group_model'),
            is.data.frame(df_i),
            is(base, "OrthogonalSplineBasis"),
            is.list(control)
        )

    df_i[!is.na(pull(df_i, variables$prev_outcome)), 'baseline_lambda'] <-
        estimate_baseline_intensity(intensity.model, df_i[!is.na(pull(df_i, variables$prev_outcome)), ])

    df.in.range <- df_i |>
        filter(
            !!min(base@knots) <= !!variables$time,
            !!variables$time <= !!max(base@knots)
        )

    term1 <- if(nrow(df.in.range) == 0) rep(0, ncol(base)) else
        df_i |>
        mutate(
            Exp_gamma = exp((!!coef(intensity.model))*!!variables$prev_outcome),
        ) |>
        filter(
            !!min(base@knots) <= !!variables$time,
            !!variables$time <= !!max(base@knots)
        ) |>
        pcori_conditional_means(
            outcome.model, alpha, new.data = _
        ) |>
        mutate(
            Term1_unweighted =
                (!!(variables$outcome)-E_Y_past)/
                (baseline_lambda*Exp_gamma* exp(-alpha*!!(variables$outcome))*E_exp_alphaY)
        ) |>
        summarize(
            term1 =
                list(crossprod(evaluate(base, .data$time), .data$Term1_unweighted)),
        ) |>
        pull(term1) |> unlist()

    expected_value <- \(data, ...){
        matrix(
            outcome.model |>
                pcori_conditional_means(..., alpha=alpha, new.data = data) |>
                pull('E_Y_past'),
            nrow = nrow(data)
        )}

    term2 <-
        if (control$integration.method == 'piecewise') {
            numerically_integrate_influence_term_2_for_one_alpha_and_one_patient_piecewise(
                df_i, expected_value=expected_value, base=base,
                resolution.within.period = control$resolution.within.period,
                variables = variables, ...
            ) |> as.vector()
        } else if (control$integration.method == 'numerical') {
            numerically_integrate_influence_term_2_for_one_alpha_and_one_patient(
                df_i, expected_value, base=base,
                resolution = control$resolution,
                variables = variables, ...
            ) |> as.vector()
        } else if (control$integration.method == 'linear') {
            compute_influence_term_2_linearly(
                df_i, expected_value, base=base,
                variables = variables, ...
            ) |>
                unlist() |> as.vector()
        } else if (control$integration.method == 'quadv') {
            compute_influence_term_2_quadv(
                df_i, expected_value, base=base,
                tol = control$tol,
                variables = variables, ...
            )
        }
    influence <- term1 + term2
    tibble(
        alpha,
        term1 = list(term1),
        term2 = list(term2),
        influence = list(influence)
    )
}
