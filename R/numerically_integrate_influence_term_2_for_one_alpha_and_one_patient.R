numerically_integrate_influence_term_2_for_one_alpha_and_one_patient <-
function(
    df_i,
    object,
    alpha,
    base,
    resolution = 1000,
    ...
){
    expected_value <- \(data, ...){
        matrix(
            pull(pcori_conditional_means(object$outcome_model, alpha=alpha, new.data = data), 'E_Y_past'),
            nrow = nrow(data)
        )
    }
    eval.times <- seq(min(base@knots), max(base@knots), length.out = resolution)
    B <- evaluate(orthogonalize(base), eval.times)
    spline_df_est <- impute_patient_df(eval.times, df_i, object)

    Ey = pull(pcori_conditional_means(object$outcome_model, alpha=alpha, new.data = spline_df_est), 'E_Y_past')

    summation <- crossprod(B, Ey) -
        crossprod(head(B,1), head(Ey,1))/2 -
        crossprod(tail(B,1), tail(Ey,1))/2
    bin.width <- (max(base@knots) - min(base@knots))/resolution

    return(summation * bin.width)
}


if(F){



}
